#!/usr/bin/env python3
"""
bin/generate_docking_box.py
Calculates exact docking box boundaries and centroid from a reference ligand,
generating a ready-to-use PyMOL script (.pml) with CGO 3D wireframe boxes.
"""

import argparse
import sys
from pathlib import Path
from typing import Tuple

FRAMEWORK_DIR = Path(__file__).resolve().parent.parent
if str(FRAMEWORK_DIR) not in sys.path:
    sys.path.insert(0, str(FRAMEWORK_DIR))

from cli.core.chem import compute_autobox_geometry


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
    auto_pad = float(args.autobox_add)
    geom = compute_autobox_geometry(ref_p, auto_pad)

    cx, cy, cz = geom["centroid"]
    gn_min_x, gn_max_x, gn_min_y, gn_max_y, gn_min_z, gn_max_z = geom["autobox_bounds"]

    inner_w, inner_h, inner_d = parse_box_dims(args.inner_box)
    outer_w, outer_h, outer_d = parse_box_dims(args.outer_box)

    # Glide Outer box
    go_min_x, go_max_x = cx - outer_w / 2.0, cx + outer_w / 2.0
    go_min_y, go_max_y = cy - outer_h / 2.0, cy + outer_h / 2.0
    go_min_z, go_max_z = cz - outer_d / 2.0, cz + outer_d / 2.0

    # Glide Inner box
    gi_min_x, gi_max_x = cx - inner_w / 2.0, cx + inner_w / 2.0
    gi_min_y, gi_max_y = cy - inner_h / 2.0, cy + inner_h / 2.0
    gi_min_z, gi_max_z = cz - inner_d / 2.0, cz + inner_d / 2.0

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

# 3. Autobox (+{auto_pad:.1f} A padding) (Orange)
autobox_cgo = {cgo_box_lines(gn_min_x, gn_max_x, gn_min_y, gn_max_y, gn_min_z, gn_max_z, 1.0, 0.6, 0.1)}
cmd.load_cgo(autobox_cgo, "{args.receptor_name}_autobox")

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
