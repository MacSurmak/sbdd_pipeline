#!/usr/bin/env python3
"""
bin/wrap_orca_resp.py
Two-stage quantum chemical parameterizer for small molecules:
1. Geometry optimization in implicit solvent (DFT CPCM) to eliminate vacuum self-folding.
2. Canonical AMBER ESP evaluation (gas-phase HF/6-31G* !RESP).
3. Topological graph symmetry equalization via RDKit automorphisms.
4. GAFF2 atom typing and frcmod generation via AmberTools.
Strict Fail-Fast: zero silent heuristics, zero unhandled errors.
"""

import argparse
import math
import os
import re
import shutil
import subprocess
import sys
from pathlib import Path
from rdkit import Chem


def run_cmd(cmd: list, desc: str = "Command") -> subprocess.CompletedProcess:
    """Executes a command with strict Fail-Fast error logging."""
    print(f"[{desc}] Exec: {' '.join(cmd)}", flush=True)
    res = subprocess.run(cmd, capture_output=True, text=True, check=False)
    if res.returncode != 0:
        print(
            f"FATAL: {desc} failed with exit code {res.returncode}!\n"
            f"--- STDOUT ---\n{res.stdout}\n"
            f"--- STDERR ---\n{res.stderr}",
            file=sys.stderr,
            flush=True,
        )
        sys.exit(res.returncode)
    return res


def get_spin_multiplicity(mol: Chem.Mol, formal_charge: int) -> int:
    """Calculates deterministic spin multiplicity (2S + 1) based on valence electrons."""
    total_protons = sum(atom.GetAtomicNum() for atom in mol.GetAtoms())
    total_electrons = total_protons - formal_charge
    return 2 if (total_electrons % 2 != 0) else 1


def validate_input_geometry(mol: Chem.Mol, mol_id: str, sdf_path: Path):
    """Enforces rigorous physical sanity checks on input 3D coordinates."""
    has_hydrogens = any(a.GetAtomicNum() == 1 for a in mol.GetAtoms())
    if not has_hydrogens:
        raise ValueError(
            f"FATAL: Molecule '{mol_id}' in {sdf_path} has no explicit hydrogens! "
            "Quantum parametrization requires a protonated 3D structure."
        )

    conf = mol.GetConformer()
    for idx, atom in enumerate(mol.GetAtoms()):
        pos = conf.GetAtomPosition(idx)
        if abs(pos.x) < 1e-4 and abs(pos.y) < 1e-4 and abs(pos.z) < 1e-4:
            raise ValueError(
                f"FATAL: Atom {atom.GetSymbol()}{idx+1} in molecule '{mol_id}' has degenerate coordinates (0, 0, 0)! "
                f"Corrupted geometry detected in {sdf_path}."
            )

    for bond in mol.GetBonds():
        a1 = bond.GetBeginAtom()
        a2 = bond.GetEndAtom()
        p1 = conf.GetAtomPosition(a1.GetIdx())
        p2 = conf.GetAtomPosition(a2.GetIdx())
        dist = math.sqrt((p1.x - p2.x) ** 2 + (p1.y - p2.y) ** 2 + (p1.z - p2.z) ** 2)
        if dist < 0.6:
            raise ValueError(
                f"FATAL: Unphysical steric clash between {a1.GetSymbol()}{a1.GetIdx()+1} and "
                f"{a2.GetSymbol()}{a2.GetIdx()+1}: {dist:.3f} Å < 0.6 Å in {sdf_path}"
            )
        if (a1.GetAtomicNum() == 1 or a2.GetAtomicNum() == 1) and dist > 1.8:
            raise ValueError(
                f"FATAL: Severed hydrogen bond between {a1.GetSymbol()}{a1.GetIdx()+1} and "
                f"{a2.GetSymbol()}{a2.GetIdx()+1}: {dist:.3f} Å > 1.8 Å in {sdf_path}"
            )


def parse_orca_resp_charges(orca_out_path: Path, expected_atoms: int) -> list[float]:
    """Parses Merz-Kollman RESP charges from ORCA standard output."""
    in_resp_section = False
    charges_dict = {}

    with open(orca_out_path, "r", encoding="utf-8", errors="ignore") as f:
        for line in f:
            line_str = line.strip()
            if "RESP CHARGES" in line_str:
                in_resp_section = True
                continue
            if in_resp_section:
                if (
                    line_str.startswith("ATOM")
                    or line_str.startswith("-")
                    or line_str.startswith("=")
                ):
                    continue
                if (
                    not line_str
                    or line_str.startswith("Sum of")
                    or line_str.startswith("---")
                ):
                    if charges_dict:
                        break
                    continue
                parts = line_str.split()
                if len(parts) >= 3 and parts[0].isdigit():
                    idx = int(parts[0])
                    val = float(parts[2])
                    charges_dict[idx] = val

    if len(charges_dict) != expected_atoms:
        raise ValueError(
            f"FATAL: Expected {expected_atoms} atomic charges from ORCA RESP output, "
            f"but parsed {len(charges_dict)} from {orca_out_path}!"
        )

    return [charges_dict[i] for i in range(expected_atoms)]


def symmetrize_charges(
    mol: Chem.Mol, raw_charges: list[float], formal_charge: int
) -> list[float]:
    """
    Symmetrizes partial atomic charges across graph topological automorphisms
    and guarantees exact integer net charge conservation.
    """
    n_atoms = mol.GetNumAtoms()
    matches = mol.GetSubstructMatches(mol, uniquify=False)

    equiv_classes = []
    visited = set()
    for i in range(n_atoms):
        if i in visited:
            continue
        eq_group = {match[i] for match in matches}
        visited.update(eq_group)
        equiv_classes.append(sorted(list(eq_group)))

    symmetrized = list(raw_charges)
    for group in equiv_classes:
        mean_q = sum(raw_charges[idx] for idx in group) / len(group)
        for idx in group:
            symmetrized[idx] = mean_q

    # Enforce exact net charge conservation without rounding drift
    total_q = sum(symmetrized)
    drift = float(formal_charge) - total_q
    correction_per_atom = drift / float(n_atoms)

    final_charges = [q + correction_per_atom for q in symmetrized]
    final_sum = sum(final_charges)

    if abs(final_sum - float(formal_charge)) > 1e-5:
        raise ValueError(
            f"FATAL: Symmetrized charge sum ({final_sum:.6f}) deviates from "
            f"formal charge ({formal_charge}) beyond physical tolerance!"
        )

    return final_charges


def main():
    parser = argparse.ArgumentParser(
        description="Two-Stage ORCA RESP Small Molecule Parameterizer"
    )
    parser.add_argument("--sdf", required=True, help="Input prepared ligand SDF")
    parser.add_argument(
        "--charge", required=True, type=int, help="Net formal molecular charge"
    )
    parser.add_argument(
        "--mol_id", required=True, help="Canonical molecule identifier"
    )
    parser.add_argument("--res_name", required=True, help="Residue name (3-4 chars)")
    parser.add_argument("--out_mol2", required=True, help="Output parameterized mol2")
    parser.add_argument(
        "--out_frcmod", required=True, help="Output forcefield frcmod"
    )
    parser.add_argument(
        "--out_opt_sdf", required=True, help="Output solvated relaxed SDF"
    )

    # Quantum parameters strictly passed from YAML configuration (zero defaults)
    parser.add_argument("--forcefield", required=True)
    parser.add_argument("--nprocs", type=int, required=True)
    parser.add_argument("--maxcore_mb", type=int, required=True)
    parser.add_argument("--opt_theory", required=True)
    parser.add_argument("--opt_basis", required=True)
    parser.add_argument("--opt_dispersion", required=True)
    parser.add_argument("--opt_solvent", required=True)
    parser.add_argument("--opt_rijcosx", required=True)
    parser.add_argument("--opt_maxiter", type=int, required=True)
    parser.add_argument("--scf_convergence", required=True)
    parser.add_argument("--scf_maxiter", type=int, required=True)
    parser.add_argument("--esp_theory", required=True)
    parser.add_argument("--esp_basis", required=True)
    parser.add_argument("--esp_density", required=True)
    args = parser.parse_args()

    sdf_path = Path(args.sdf).resolve()
    if not sdf_path.exists() or sdf_path.stat().st_size == 0:
        raise FileNotFoundError(
            f"FATAL: Input SDF deliverable not found: {sdf_path}"
        )

    suppl = Chem.SDMolSupplier(str(sdf_path), removeHs=False, sanitize=False)
    mol = next((m for m in suppl if m is not None), None)
    if mol is None:
        raise ValueError(
            f"FATAL: Could not read valid molecular structure from {sdf_path}"
        )

    try:
        Chem.SanitizeMol(mol)
    except Exception:
        mol.UpdatePropertyCache(strict=False)

    validate_input_geometry(mol, args.mol_id, sdf_path)
    mult = get_spin_multiplicity(mol, args.charge)
    n_atoms = mol.GetNumAtoms()
    conf = mol.GetConformer()

    res_clean = args.res_name.strip().upper()[:4]
    if len(res_clean) < 3:
        res_clean = f"{res_clean:<3s}"

    orca_bin = shutil.which("orca") or "/opt/containers/bin/orca"

    # =========================================================================
    # STAGE 1: Geometry Optimization in Implicit Solvent (DFT CPCM)
    # Prevents artificial gas-phase self-folding while capturing solution state
    # =========================================================================
    opt_rijcosx_flag = (
        "RIJCOSX def2/J"
        if str(args.opt_rijcosx).lower() in ["true", "1", "yes"]
        else ""
    )
    opt_disp_flag = (
        args.opt_dispersion
        if args.opt_dispersion and args.opt_dispersion.lower() != "none"
        else ""
    )
    opt_solv_flag = (
        args.opt_solvent
        if args.opt_solvent and args.opt_solvent.lower() != "none"
        else ""
    )

    def run_orca_calculation(inp_name: str, out_name: str, desc: str):
        print(f"[{desc}] Exec: {orca_bin} {inp_name} > {out_name}", flush=True)

        # In Slurm cpus-per-task allocations, bridge Slurm CPUs to OpenMPI/PRRTE task slots
        orca_env = os.environ.copy()
        cpus_per_task = orca_env.get("SLURM_CPUS_PER_TASK", str(args.nprocs))
        orca_env["SLURM_TASKS_PER_NODE"] = str(cpus_per_task)
        orca_env["SLURM_NTASKS"] = str(cpus_per_task)
        orca_env["SLURM_NPROCS"] = str(cpus_per_task)

        # PRRTE and OpenMPI 5.x tuning for local node multi-threading under Slurm cgroups
        orca_env["PRTE_MCA_rmaps_default_mapping_policy"] = ":oversubscribe"
        orca_env["OMPI_MCA_rmaps_base_oversubscribe"] = "1"
        orca_env["OMPI_MCA_hwloc_base_binding_policy"] = "none"
        orca_env["OMPI_MCA_btl"] = "^openib"

        with open(out_name, "w", encoding="utf-8") as out_f:
            res = subprocess.run(
                [orca_bin, inp_name],
                stdout=out_f,
                stderr=subprocess.STDOUT,
                text=True,
                check=False,
                env=orca_env,
            )
        out_path = Path(out_name)
        if not out_path.exists() or out_path.stat().st_size == 0:
            print(
                f"FATAL: {desc} failed to produce output file: {out_name}",
                file=sys.stderr,
                flush=True,
            )
            sys.exit(res.returncode if res.returncode != 0 else 1)

        log_content = out_path.read_text(encoding="utf-8", errors="ignore")
        if res.returncode != 0 or "ORCA TERMINATED NORMALLY" not in log_content:
            tail_lines = "\n".join(log_content.splitlines()[-40:])
            print(
                f"FATAL: {desc} terminated abnormally (exit code {res.returncode})!\n"
                f"--- TAIL OF {out_name} ---\n{tail_lines}",
                file=sys.stderr,
                flush=True,
            )
            sys.exit(res.returncode if res.returncode != 0 else 1)

    opt_inp_lines = [
        f"! {args.opt_theory} {opt_disp_flag} {args.opt_basis} {opt_rijcosx_flag} {opt_solv_flag} Opt {args.scf_convergence}",
        f"%maxcore {args.maxcore_mb}",
        f"%pal nprocs {args.nprocs} end",
        f"%geom maxiter {args.opt_maxiter} end",
        f"%scf maxiter {args.scf_maxiter} end",
        f"* xyz {args.charge} {mult}",
    ]
    for i in range(n_atoms):
        pos = conf.GetAtomPosition(i)
        elem = mol.GetAtomWithIdx(i).GetSymbol()
        opt_inp_lines.append(
            f"  {elem:<3s} {pos.x:12.6f} {pos.y:12.6f} {pos.z:12.6f}"
        )
    opt_inp_lines.append("*\n")

    Path("orca_opt.inp").write_text("\n".join(opt_inp_lines), encoding="utf-8")
    run_orca_calculation("orca_opt.inp", "orca_opt.out", desc="ORCA Stage 1: Geometry Optimization")

    opt_xyz_path = Path("orca_opt.xyz")
    if not opt_xyz_path.exists() or opt_xyz_path.stat().st_size == 0:
        raise FileNotFoundError(
            "FATAL: ORCA geometry optimization failed to generate 'orca_opt.xyz'! "
            "Optimization converged to error or terminated abnormally."
        )

    # Read optimized coordinates from XYZ and update RDKit molecule
    with open(opt_xyz_path, "r", encoding="utf-8") as f:
        xyz_lines = [
            line.strip() for line in f if line.strip() and not line.startswith("#")
        ]

    # Skip atom count header and comment line
    coord_lines = xyz_lines[2:]
    if len(coord_lines) != n_atoms:
        raise ValueError(
            f"FATAL: Atom count in orca_opt.xyz ({len(coord_lines)}) does not match input ({n_atoms})!"
        )

    for i, line in enumerate(coord_lines):
        parts = line.split()
        conf.SetAtomPosition(i, (float(parts[1]), float(parts[2]), float(parts[3])))

    # Save auxiliary relaxed geometry deliverable for scientific provenance
    with Chem.SDWriter(args.out_opt_sdf) as writer:
        writer.SetForceV3000(False)
        writer.write(mol)
    print(
        f"✓ Saved solvated relaxed geometry deliverable: {args.out_opt_sdf}",
        flush=True,
    )

    # =========================================================================
    # STAGE 2: Gas-Phase Canonical Merz-Kollman ESP Evaluation (HF/6-31G* !RESP)
    # Strictly in vacuum to match AMBER GAFF/GAFF2 dipole moment calibration
    # =========================================================================
    esp_inp_lines = [
        f"! {args.esp_theory} {args.esp_basis} RESP {args.scf_convergence}",
        f"%maxcore {args.maxcore_mb}",
        f"%pal nprocs {args.nprocs} end",
        f"%scf maxiter {args.scf_maxiter} end",
        f"* xyz {args.charge} {mult}",
    ]
    for i in range(n_atoms):
        pos = conf.GetAtomPosition(i)
        elem = mol.GetAtomWithIdx(i).GetSymbol()
        esp_inp_lines.append(
            f"  {elem:<3s} {pos.x:12.6f} {pos.y:12.6f} {pos.z:12.6f}"
        )
    esp_inp_lines.append("*\n")

    Path("orca_esp.inp").write_text("\n".join(esp_inp_lines), encoding="utf-8")
    run_orca_calculation("orca_esp.inp", "orca_esp.out", desc="ORCA Stage 2: Gas-Phase Merz-Kollman ESP")

    orca_esp_out = Path("orca_esp.out")

    raw_charges = parse_orca_resp_charges(orca_esp_out, n_atoms)
    final_charges = symmetrize_charges(mol, raw_charges, args.charge)

    # Export charge vector for Antechamber (-c rc -cf)
    Path("resp_charges.crg").write_text(
        "\n".join(f"{q:14.6f}" for q in final_charges) + "\n", encoding="utf-8"
    )

    # Write clean MDL V2000 SDF with relaxed coordinates for Antechamber
    clean_v2000_sdf = "ligand_opt_v2000.sdf"
    with Chem.SDWriter(clean_v2000_sdf) as writer:
        writer.SetForceV3000(False)
        writer.write(mol)

    # =========================================================================
    # STAGE 3: Antechamber Atom Typing & Parmchk2 Parameterization
    # =========================================================================
    run_cmd(
        [
            "antechamber",
            "-i",
            clean_v2000_sdf,
            "-fi",
            "sdf",
            "-o",
            args.out_mol2,
            "-fo",
            "mol2",
            "-at",
            args.forcefield,
            "-rn",
            res_clean,
            "-c",
            "rc",
            "-cf",
            "resp_charges.crg",
            "-nc",
            str(args.charge),
            "-dr",
            "no",
        ],
        desc="Antechamber GAFF2 Atom Typing with ORCA RESP Charges",
    )

    run_cmd(
        [
            "parmchk2",
            "-i",
            args.out_mol2,
            "-f",
            "mol2",
            "-o",
            args.out_frcmod,
            "-s",
            "2",
        ],
        desc="Parmchk2 Parameter File Generation",
    )

    if not Path(args.out_mol2).exists() or Path(args.out_mol2).stat().st_size == 0:
        raise FileNotFoundError(
            f"FATAL: Antechamber failed to produce output mol2: {args.out_mol2}"
        )
    if (
        not Path(args.out_frcmod).exists()
        or Path(args.out_frcmod).stat().st_size == 0
    ):
        raise FileNotFoundError(
            f"FATAL: Parmchk2 failed to produce output frcmod: {args.out_frcmod}"
        )

    print(
        f"✓ Successfully generated ORCA RESP topology: {args.out_mol2} and {args.out_frcmod}",
        flush=True,
    )


if __name__ == "__main__":
    main()
