#!/usr/bin/env python3
"""
declash_gmx_molecules_v2.py

Rigid-body de-clashing for GROMACS systems that crash with "infinite force" or
LINCS warnings at step 0 due to catastrophic INTER-molecular overlaps.

Why v2
------
Compared to a residue-based de-clasher, this script moves WHOLE MOLECULES as
defined by the topology [ molecules ] section. This is critical for Amber/Lipid
multi-residue lipids (PC+PA+OL etc.), where constrained/bonded atoms can span
multiple residues. Moving residues independently can stretch those bonds and
cause LINCS failures with enormous energies.

Compared to v1, v2 uses a *sequential* (Gauss–Seidel-like) resolution of the
worst overlaps each iteration. This reduces the "vector cancellation" effect
where many small clashes overwhelm a few catastrophic ones.

Recommended workflow
--------------------
1) De-clash coordinates:
   python declash_gmx_molecules_v2.py -f step3_input.gro -p topol.top -o step3_declashed.gro \
       --fix-moltype system1 --fix-moltype UNK --cutoff 0.10 --niter 120 --max-move 0.10 \
       --report declash_report.txt -v
       
       python3 declash_gmx_molecules_v2.py -f step3_input.gro -p topol.top -o step3_declashed.gro \
       --fix-moltype system1 --fix-moltype UNK --no-heavy-only --cutoff 0.10 --niter 300 --max-move 0.05 \
       --report declash_report.txt -v

2) First EM "unbrick" with *constraints disabled*:
   constraints = none

3) Then re-enable constraints / position restraints and proceed normally.
"""

from __future__ import annotations

import argparse
import math
import re
from dataclasses import dataclass
from pathlib import Path
from typing import Dict, List, Tuple

import numpy as np


def read_text_with_includes(path: Path, max_depth: int = 20) -> str:
    """Recursively expand #include lines if the include file exists locally."""
    seen: set[Path] = set()

    def _read(p: Path, depth: int) -> str:
        if depth > max_depth:
            raise RuntimeError(f"Exceeded include depth {max_depth} while reading {p}")
        p = p.resolve()
        if p in seen:
            return ""
        seen.add(p)
        base = p.parent
        out_lines: List[str] = []
        for line in p.read_text(errors="ignore").splitlines():
            m = re.match(r'^\s*#include\s+["<]([^">]+)[">]\s*$', line)
            if not m:
                out_lines.append(line)
                continue
            inc = m.group(1)
            cand = base / inc
            if cand.exists():
                out_lines.append(_read(cand, depth + 1))
            else:
                out_lines.append(line)
        return "\n".join(out_lines)

    return _read(path, 0)


def read_gro(path: Path) -> Tuple[str, np.ndarray, np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    """Return (title, coords[n,3], resid[n], resname[n], atomname[n], box[3])."""
    lines = path.read_text(errors="ignore").splitlines()
    if len(lines) < 3:
        raise ValueError(f"{path} does not look like a .gro file.")
    title = lines[0].rstrip("\n")
    nat = int(lines[1].strip())
    if len(lines) < 3 + nat:
        raise ValueError(f"{path}: declared {nat} atoms, but file is truncated.")
    coords = np.zeros((nat, 3), dtype=np.float64)
    resid = np.zeros(nat, dtype=np.int32)
    resname = np.empty(nat, dtype=object)
    atomname = np.empty(nat, dtype=object)

    for i in range(nat):
        l = lines[2 + i]
        resid[i] = int(l[0:5])
        resname[i] = l[5:10].strip()
        atomname[i] = l[10:15].strip()
        coords[i, 0] = float(l[20:28])
        coords[i, 1] = float(l[28:36])
        coords[i, 2] = float(l[36:44])

    box = np.array(list(map(float, lines[2 + nat].split()[:3])), dtype=np.float64)
    return title, coords, resid, resname, atomname, box


def write_gro(path: Path, title: str, coords: np.ndarray, resid: np.ndarray, resname: np.ndarray, atomname: np.ndarray, box: np.ndarray):
    nat = coords.shape[0]
    with path.open("w") as f:
        f.write(f"{title}\n")
        f.write(f"{nat:5d}\n")
        for i in range(nat):
            f.write(f"{resid[i]:5d}{str(resname[i])[:5]:<5}{str(atomname[i])[:5]:>5}{(i+1)%100000:5d}"
                    f"{coords[i,0]:8.3f}{coords[i,1]:8.3f}{coords[i,2]:8.3f}\n")
        f.write(f"{box[0]:10.5f}{box[1]:10.5f}{box[2]:10.5f}\n")


def parse_moleculetype_atomcounts(top_text: str) -> Dict[str, int]:
    lines = top_text.splitlines()
    mt_idx = [i for i, l in enumerate(lines) if re.match(r"^\s*\[\s*moleculetype\s*\]\s*$", l, re.I)]
    out: Dict[str, int] = {}
    for idx in mt_idx:
        j = idx + 1
        while j < len(lines) and (lines[j].strip() == "" or lines[j].lstrip().startswith(";")):
            j += 1
        if j >= len(lines):
            continue
        name = lines[j].split()[0]

        k = j + 1
        atoms_start = None
        while k < len(lines):
            if re.match(r"^\s*\[\s*moleculetype\s*\]\s*$", lines[k], re.I):
                break
            if re.match(r"^\s*\[\s*system\s*\]\s*$", lines[k], re.I):
                break
            if re.match(r"^\s*\[\s*atoms\s*\]\s*$", lines[k], re.I):
                atoms_start = k + 1
                break
            k += 1
        if atoms_start is None:
            continue

        n = 0
        k = atoms_start
        while k < len(lines):
            if re.match(r"^\s*\[", lines[k]):
                break
            s = lines[k].strip()
            if s and (not s.startswith(";")) and (not s.startswith("#")):
                n += 1
            k += 1
        out[name] = n
    return out


def parse_molecules(top_text: str) -> List[Tuple[str, int]]:
    lines = top_text.splitlines()
    idx = None
    for i, l in enumerate(lines):
        if re.match(r"^\s*\[\s*molecules\s*\]\s*$", l, re.I):
            idx = i
            break
    if idx is None:
        raise ValueError("Topology does not contain a [ molecules ] section.")
    out: List[Tuple[str, int]] = []
    i = idx + 1
    while i < len(lines):
        l = lines[i]
        if re.match(r"^\s*\[", l):
            break
        s = l.strip()
        if s and (not s.startswith(";")) and (not s.startswith("#")):
            parts = s.split()
            if len(parts) >= 2:
                out.append((parts[0], int(float(parts[1]))))
        i += 1
    if not out:
        raise ValueError("[ molecules ] section is empty.")
    return out


@dataclass
class MolRange:
    moltype: str
    start: int
    end: int


def build_molecule_ranges(natoms: int, atomcounts: Dict[str, int], mol_list: List[Tuple[str, int]]) -> List[MolRange]:
    ranges: List[MolRange] = []
    idx = 0
    for mt, count in mol_list:
        if mt not in atomcounts:
            raise ValueError(f"Cannot find [ moleculetype ] {mt} in topology.")
        n = atomcounts[mt]
        for _ in range(count):
            ranges.append(MolRange(mt, idx, idx + n))
            idx += n
    if idx != natoms:
        raise ValueError(f"Topology expands to {idx} atoms, but .gro has {natoms}.")
    return ranges


def wrap(coords: np.ndarray, box: np.ndarray) -> np.ndarray:
    return coords - np.floor(coords / box) * box


def min_image(d: np.ndarray, box: np.ndarray) -> np.ndarray:
    return d - box * np.round(d / box)


def is_heavy(atomname: str) -> bool:
    a = atomname.strip().upper()
    if a == "EPW":
        return False
    return not (a.startswith("H") or a.startswith("D"))


def find_clashes(coords: np.ndarray, box: np.ndarray, atom_molid: np.ndarray, atomname: np.ndarray,
                 cutoff: float, max_pairs: int, heavy_only: bool) -> List[Tuple[int, int, np.ndarray, float]]:
    L = box.astype(np.float64)
    coords = coords.astype(np.float64)

    if heavy_only:
        mask = np.array([is_heavy(a) for a in atomname], dtype=bool)
    else:
        mask = np.array([a.strip().upper() != "EPW" for a in atomname], dtype=bool)

    idxs = np.nonzero(mask)[0]
    pos = wrap(coords[idxs], L)
    mol = atom_molid[idxs]

    nx = max(int(L[0] / cutoff), 1)
    ny = max(int(L[1] / cutoff), 1)
    nz = max(int(L[2] / cutoff), 1)
    sx, sy, sz = L[0] / nx, L[1] / ny, L[2] / nz

    cx = (pos[:, 0] / sx).astype(np.int32) % nx
    cy = (pos[:, 1] / sy).astype(np.int32) % ny
    cz = (pos[:, 2] / sz).astype(np.int32) % nz
    cell_id = (cx * ny + cy) * nz + cz

    order = np.argsort(cell_id)
    cell_id = cell_id[order]
    pos = pos[order]
    mol = mol[order]
    idxs = idxs[order]

    unique, start = np.unique(cell_id, return_index=True)
    end = np.r_[start[1:], len(cell_id)]
    slices = {int(u): (int(s), int(e)) for u, s, e in zip(unique, start, end)}
    neigh = [(dx, dy, dz) for dx in (-1, 0, 1) for dy in (-1, 0, 1) for dz in (-1, 0, 1)]
    c2 = cutoff * cutoff

    out: List[Tuple[int, int, np.ndarray, float]] = []
    for u, st, en in zip(unique, start, end):
        u = int(u); st = int(st); en = int(en)
        cz_u = u % nz
        tmp = u // nz
        cy_u = tmp % ny
        cx_u = tmp // ny
        for dx, dy, dz in neigh:
            u_n = (((cx_u + dx) % nx) * ny + ((cy_u + dy) % ny)) * nz + ((cz_u + dz) % nz)
            if u_n not in slices:
                continue
            st2, en2 = slices[u_n]
            if u_n < u:
                continue
            for ii in range(st, en):
                jj0 = st2 if u_n != u else ii + 1
                mi = int(mol[ii])
                for jj in range(jj0, en2):
                    mj = int(mol[jj])
                    if mi == mj:
                        continue
                    d = pos[jj] - pos[ii]
                    d = min_image(d, L)
                    r2 = float(d.dot(d))
                    if r2 < c2:
                        r = math.sqrt(r2) if r2 > 1e-16 else 0.0
                        uvec = d / r if r > 1e-6 else np.array([1.0, 0.0, 0.0], dtype=np.float64)
                        out.append((int(idxs[ii]), int(idxs[jj]), uvec, float(r)))
                        if len(out) >= max_pairs:
                            return out
    return out


def declash_sequential(coords: np.ndarray, box: np.ndarray, mol_ranges: List[MolRange],
                       atom_molid: np.ndarray, atomname: np.ndarray,
                       mobile_mol: np.ndarray, restrict_xy_mol: np.ndarray,
                       cutoff: float, niter: int, damping: float, max_move: float,
                       max_pairs: int, heavy_only: bool, seed: int, verbose: bool):
    rng = np.random.default_rng(seed if seed >= 0 else None)
    L = box.astype(np.float64)
    coords = wrap(coords.astype(np.float64), L)
    nmol = len(mol_ranges)
    mol_slices = [slice(m.start, m.end) for m in mol_ranges]

    for it in range(niter):
        clashes = find_clashes(coords, L, atom_molid, atomname, cutoff, max_pairs=max_pairs, heavy_only=heavy_only)
        if not clashes:
            if verbose:
                print(f"iter {it:3d}: clashes=0")
            return coords, {"cutoff_nm": float(cutoff), "iterations": float(it), "final_clashes": 0.0}

        clashes.sort(key=lambda x: x[3])
        worst_r = clashes[0][3]
        moved = np.zeros((nmol, 3), dtype=np.float64)
        applied = 0

        def cap(mid: int, d: np.ndarray) -> np.ndarray:
            if not mobile_mol[mid]:
                return np.zeros(3, dtype=np.float64)
            trial = moved[mid] + d
            if np.linalg.norm(trial) <= max_move:
                return d
            rem = max_move - float(np.linalg.norm(moved[mid]))
            if rem <= 0:
                return np.zeros(3, dtype=np.float64)
            dn = float(np.linalg.norm(d))
            if dn < 1e-12:
                return np.zeros(3, dtype=np.float64)
            return d * (rem / dn)

        for ai, aj, uvec, r in clashes:
            mi = int(atom_molid[ai]); mj = int(atom_molid[aj])
            if (not mobile_mol[mi]) and (not mobile_mol[mj]):
                continue
            if r < 1e-6:
                uvec = rng.normal(size=3); uvec /= (np.linalg.norm(uvec) + 1e-12)

            push = (cutoff - r) * damping
            if push <= 0:
                continue

            if mobile_mol[mi] and mobile_mol[mj]:
                di = -0.5 * push * uvec
                dj = +0.5 * push * uvec
            elif mobile_mol[mi]:
                di = -push * uvec
                dj = np.zeros(3, dtype=np.float64)
            else:
                di = np.zeros(3, dtype=np.float64)
                dj = +push * uvec

            di = cap(mi, di); dj = cap(mj, dj)
            if restrict_xy_mol[mi]:
                di[2] = 0.0
            if restrict_xy_mol[mj]:
                dj[2] = 0.0

            if np.any(di):
                coords[mol_slices[mi]] += di; moved[mi] += di; applied += 1
            if np.any(dj):
                coords[mol_slices[mj]] += dj; moved[mj] += dj; applied += 1

        coords = wrap(coords, L)

        if verbose:
            print(f"iter {it:3d}: clashes={len(clashes):6d} applied={applied:6d} worst={worst_r:8.4f} nm")

        if applied == 0:
            # escape local trap
            jitter = rng.normal(size=(nmol, 3)) * 0.02
            jitter[~mobile_mol] = 0.0
            jitter[restrict_xy_mol, 2] = 0.0
            for mid in range(nmol):
                if mobile_mol[mid]:
                    coords[mol_slices[mid]] += jitter[mid]
            coords = wrap(coords, L)

    final = find_clashes(coords, L, atom_molid, atomname, cutoff, max_pairs=500000, heavy_only=heavy_only)
    return coords, {"cutoff_nm": float(cutoff), "iterations": float(niter), "final_clashes": float(len(final))}


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("-f", "--gro", required=True, type=Path)
    ap.add_argument("-p", "--top", required=True, type=Path)
    ap.add_argument("-o", "--out", required=True, type=Path)
    ap.add_argument("--cutoff", type=float, default=0.10)
    ap.add_argument("--niter", type=int, default=120)
    ap.add_argument("--damping", type=float, default=0.8)
    ap.add_argument("--max-move", type=float, default=0.10)
    ap.add_argument("--max-pairs", type=int, default=5000)
    ap.add_argument("--scale", type=float, default=1.0)
    ap.add_argument("--no-heavy-only", action="store_true")
    ap.add_argument("--seed", type=int, default=1)
    ap.add_argument("--fix-moltype", action="append", default=[])
    ap.add_argument("--fix-first-n", type=int, default=0)
    ap.add_argument("--restrict-xy-moltype", action="append", default=[])
    ap.add_argument("--report", type=Path, default=None)
    ap.add_argument("-v", "--verbose", action="store_true")
    args = ap.parse_args()

    title, coords, resid, resname, atomname, box = read_gro(args.gro)

    if args.scale and abs(args.scale - 1.0) > 1e-8:
        sc = float(args.scale)
        center = box / 2.0
        coords = wrap(coords, box)
        coords = (coords - center) * sc + center
        box = box * sc
        coords = wrap(coords, box)

    top_text = read_text_with_includes(args.top)
    atomcounts = parse_moleculetype_atomcounts(top_text)
    mol_list = parse_molecules(top_text)
    mol_ranges = build_molecule_ranges(coords.shape[0], atomcounts, mol_list)

    atom_molid = np.empty(coords.shape[0], dtype=np.int32)
    moltype_of_molid: List[str] = []
    for mid, m in enumerate(mol_ranges):
        atom_molid[m.start:m.end] = mid
        moltype_of_molid.append(m.moltype)
    moltype_of_molid = np.array(moltype_of_molid, dtype=object)

    mobile = np.ones(len(mol_ranges), dtype=bool)
    if args.fix_first_n > 0:
        mobile[:args.fix_first_n] = False
    for mt in args.fix_moltype:
        mobile[moltype_of_molid == mt] = False

    restrict_xy = np.zeros(len(mol_ranges), dtype=bool)
    for mt in args.restrict_xy_moltype:
        restrict_xy[moltype_of_molid == mt] = True

    heavy_only = not args.no_heavy_only

    new_coords, stats = declash_sequential(coords, box, mol_ranges, atom_molid, atomname,
                                          mobile, restrict_xy, args.cutoff, args.niter,
                                          args.damping, args.max_move, args.max_pairs,
                                          heavy_only, args.seed, args.verbose)

    write_gro(args.out, title + " (declashed)", new_coords, resid, resname, atomname, box)

    if args.report is not None:
        pairs = find_clashes(new_coords, box, atom_molid, atomname, args.cutoff, max_pairs=500000, heavy_only=True)
        pairs = [(ai, aj, r) for ai, aj, _, r in pairs]
        pairs.sort(key=lambda x: x[2])
        with args.report.open("w") as f:
            f.write(f"declash stats: {stats}\n")
            f.write(f"worst inter-molecular heavy-atom contacts < {args.cutoff} nm:\n")
            for ai, aj, r in pairs[:50]:
                f.write(f"{r:10.5f} nm  atom {ai+1:8d}  atom {aj+1:8d}\n")


if __name__ == "__main__":
    main()
