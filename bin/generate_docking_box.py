#!/usr/bin/env python3
"""
bin/generate_docking_box.py
Calculates exact docking box boundaries and centroid from a reference ligand,
generating a ready-to-use PyMOL script (.pml) with CGO 3D wireframe boxes.
"""

import argparse
import sys
from pathlib import Path
from typing import List, Tuple


def get_coords_rdkit(path: Path) -> List[Tuple[float, float, float]]:
    from rdkit import Chem

    mols = []
    if path.suffix.lower() == ".sdf":
        suppl = Chem.SDMolSupplier(str(path), sanitize=False, removeHs=False)
        mols = [m for m in suppl if m is not None]
    elif path.suffix.lower() == ".pdb":
        m = Chem.MolFromPDBFile(str(path), sanitize=False, removeHs=False)
        if m:
            mols = [m]

    if not mols:
        return []

    conf = mols[0].GetConformer()
    return [
        (
            float(conf.GetAtomPosition(i).x),
            float(conf.GetAtomPosition(i).y),
            float(conf.GetAtomPosition(i).z),
        )
        for i in range(mols[0].GetNumAtoms())
    ]


def get_coords_schrodinger(path: Path) -> List[Tuple[float, float, float]]:
    from schrodinger import structure

    st = structure.StructureReader.read(str(path))
    return [(float(a.x), float(a.y), float(a.z)) for a in st.atom]


def parse_box_dims(val: str) -> Tuple[float, float, float]:
    parts = [float(x.strip()) for x in val.replace(",", " ").split() if x.strip()]
    if len(parts) == 3:
        return (parts[0], parts[1], parts[2])
    raise ValueError(f"FATAL: Invalid box dimensions string: {val}")


def cgo_box_lines(
    min_x: float,
    max_x: float,
    min_y: float,
    max_y: float,
    min_z: float,
    max_z: float,
    r: float,
    g: float,
    b: float,
) -> str:
    """Generates Python CGO lines representation for PyMOL."""
    return f"""[
    BEGIN, LINES,
    COLOR, {r:.2f}, {g:.2f}, {b:.2f},
    VERTEX, {min_x:.3f}, {min_y:.3f}, {min_z:.3f},  VERTEX, {max_x:.3f}, {min_y:.3f}, {min_z:.3f},
    VERTEX, {max_x:.3f}, {min_y:.3f}, {min_z:.3f},  VERTEX, {max_x:.3f}, {max_y:.3f}, {min_z:.3f},
    VERTEX, {max_x:.3f}, {max_y:.3f}, {min_z:.3f},  VERTEX, {min_x:.3f}, {max_y:.3f}, {min_z:.3f},
    VERTEX, {min_x:.3f}, {max_y:.3f}, {min_z:.3f},  VERTEX, {min_x:.3f}, {min_y:.3f}, {min_z:.3f},
    VERTEX, {min_x:.3f}, {min_y:.3f}, {max_z:.3f},  VERTEX, {max_x:.3f}, {min_y:.3f}, {max_z:.3f},
    VERTEX, {max_x:.3f}, {min_y:.3f}, {max_z:.3f},  VERTEX, {max_x:.3f}, {max_y:.3f}, {max_z:.3f},
    VERTEX, {max_x:.3f}, {max_y:.3f}, {max_z:.3f},  VERTEX, {min_x:.3f}, {max_y:.3f}, {max_z:.3f},
    VERTEX, {min_x:.3f}, {max_y:.3f}, {max_z:.3f},  VERTEX, {min_x:.3f}, {min_y:.3f}, {max_z:.3f},
    VERTEX, {min_x:.3f}, {min_y:.3f}, {min_z:.3f},  VERTEX, {min_x:.3f}, {min_y:.3f}, {max_z:.3f},
    VERTEX, {max_x:.3f}, {min_y:.3f}, {min_z:.3f},  VERTEX, {max_x:.3f}, {min_y:.3f}, {max_z:.3f},
    VERTEX, {max_x:.3f}, {max_y:.3f}, {min_z:.3f},  VERTEX, {max_x:.3f}, {max_y:.3f}, {max_z:.3f},
    VERTEX, {min_x:.3f}, {max_y:.3f}, {min_z:.3f},  VERTEX, {min_x:.3f}, {max_y:.3f}, {max_z:.3f},
    END
]"""


def main():
    parser = argparse.ArgumentParser(
        description="Generate PyMOL docking box representation script"
    )
    parser.add_argument("--reference", required=True, help="Reference ligand file")
    parser.add_argument("--receptor_name", default="receptor", help="Receptor name")
    parser.add_argument(
        "--inner_box", default="10, 10, 10", help="Glide inner box dimensions"
    )
    parser.add_argument(
        "--outer_box", default="26, 26, 26", help="Glide outer box dimensions"
    )
    parser.add_argument(
        "--autobox_add",
        type=float,
        default=4.0,
        help="GNINA autobox padding in Angstroms",
    )
    parser.add_argument("--output", required=True, help="Output .pml script path")
    args = parser.parse_args()

    ref_p = Path(args.reference).resolve()
    coords = []
    
    if ref_p.suffix.lower() in [".sdf", ".pdb"]:
        coords = get_coords_rdkit(ref_p)
    elif ref_p.suffix.lower() in [".mae", ".maegz"]:
        coords = get_coords_schrodinger(ref_p)
    else:
        raise ValueError(f"FATAL: Unsupported reference format for docking box: {ref_p.suffix}")

    if not coords:
        print(f"FATAL: No atoms found in reference ligand {ref_p}", file=sys.stderr)
        sys.exit(1)

    # 1. Centroid and bounding box
    xs = [c[0] for c in coords]
    ys = [c[1] for c in coords]
    zs = [c[2] for c in coords]

    cx = sum(xs) / len(xs)
    cy = sum(ys) / len(ys)
    cz = sum(zs) / len(zs)

    # 2. Dimensions
    inner_w, inner_h, inner_d = parse_box_dims(args.inner_box)
    outer_w, outer_h, outer_d = parse_box_dims(args.outer_box)
    auto_pad = float(args.autobox_add)

    # Glide Outer box
    go_min_x, go_max_x = cx - outer_w / 2.0, cx + outer_w / 2.0
    go_min_y, go_max_y = cy - outer_h / 2.0, cy + outer_h / 2.0
    go_min_z, go_max_z = cz - outer_d / 2.0, cz + outer_d / 2.0

    # Glide Inner box
    gi_min_x, gi_max_x = cx - inner_w / 2.0, cx + inner_w / 2.0
    gi_min_y, gi_max_y = cy - inner_h / 2.0, cy + inner_h / 2.0
    gi_min_z, gi_max_z = cz - inner_d / 2.0, cz + inner_d / 2.0

    # GNINA Autobox
    gn_min_x, gn_max_x = min(xs) - auto_pad, max(xs) + auto_pad
    gn_min_y, gn_max_y = min(ys) - auto_pad, max(ys) + auto_pad
    gn_min_z, gn_max_z = min(zs) - auto_pad, max(zs) + auto_pad

    pml_content = f"""# PyMOL Docking Grid & Autobox Visualization
# Receptor: {args.receptor_name}
# Centroid: ({cx:.3f}, {cy:.3f}, {cz:.3f})

python
from pymol.cgo import *
from pymol import cmd

# 1. Glide Outer Search Space (Cyan)
outer_cgo = {cgo_box_lines(go_min_x, go_max_x, go_min_y, go_max_y, go_min_z, go_max_z, 0.2, 0.8, 1.0)}
cmd.load_cgo(outer_cgo, "{args.receptor_name}_glide_outer_box")

# 2. Glide Inner Centroid Space (Green)
inner_cgo = {cgo_box_lines(gi_min_x, gi_max_x, gi_min_y, gi_max_y, gi_min_z, gi_max_z, 0.2, 1.0, 0.3)}
cmd.load_cgo(inner_cgo, "{args.receptor_name}_glide_inner_box")

# 3. GNINA Autobox (+{auto_pad:.1f} A padding) (Orange)
gnina_cgo = {cgo_box_lines(gn_min_x, gn_max_x, gn_min_y, gn_max_y, gn_min_z, gn_max_z, 1.0, 0.6, 0.1)}
cmd.load_cgo(gnina_cgo, "{args.receptor_name}_gnina_autobox")

# 4. Docking Center Marker
cmd.pseudoatom("{args.receptor_name}_docking_center", pos=[{cx:.3f}, {cy:.3f}, {cz:.3f}])
cmd.show("spheres", "{args.receptor_name}_docking_center")
cmd.set("sphere_scale", 0.5, "{args.receptor_name}_docking_center")
cmd.color("yellow", "{args.receptor_name}_docking_center")

cmd.zoom("{args.receptor_name}_glide_outer_box", buffer=5.0)
python end
"""

    out_path = Path(args.output)
    out_path.parent.mkdir(parents=True, exist_ok=True)
    out_path.write_text(pml_content, encoding="utf-8")
    print(f"Generated PyMOL docking box script: {out_path}")


if __name__ == "__main__":
    main()
