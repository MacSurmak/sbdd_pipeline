#!/usr/bin/env python3
"""
File: bin/apply_membrane.py
Description: Injects implicit membrane slab properties for Schrödinger Prime MM-GBSA into MAE structure.
"""
import argparse
import sys
from pathlib import Path
from schrodinger import structure


def main():
    parser = argparse.ArgumentParser(description="Inject Implicit Membrane Properties")
    parser.add_argument("--input", required=True, help="Input prepared receptor MAE")
    parser.add_argument(
        "--info", required=True, help="Input membrane thickness text or json file"
    )
    parser.add_argument(
        "--output", required=True, help="Output receptor MAE with implicit membrane"
    )
    args = parser.parse_args()

    try:
        st = structure.StructureReader.read(args.input)

        info_path = Path(args.info)
        if not info_path.exists():
            raise FileNotFoundError(
                f"FATAL: Membrane info file {args.info} not found. Cannot apply implicit membrane."
            )

        content = info_path.read_text().strip()
        if not content:
            raise ValueError(f"FATAL: Membrane info file {args.info} is empty.")

        if content.startswith("{"):
            import json
            parsed = json.loads(content)
            if "thickness" not in parsed:
                raise KeyError(f"FATAL: 'thickness' key missing in membrane JSON {args.info}")
            thickness = float(parsed["thickness"])
        else:
            thickness = float(content)
            
        print(f"Applying membrane thickness: {thickness:.2f} Å", flush=True)

        half_t = thickness / 2.0
        st.property["r_psp_Memb1_x"] = 0.0
        st.property["r_psp_Memb1_y"] = 0.0
        st.property["r_psp_Memb1_z"] = -half_t

        st.property["r_psp_Memb2_x"] = 0.0
        st.property["r_psp_Memb2_y"] = 0.0
        st.property["r_psp_Memb2_z"] = half_t

        st.property["r_psp_Prime_Membrane_Thickness"] = thickness
        st.property["r_psp_Prime_Membrane_Buffer"] = 2.5

        st.write(args.output)
        print(
            f"Saved receptor with implicit membrane to {Path(args.output).name}",
            flush=True,
        )

    except Exception as e:
        print(f"Error applying implicit membrane: {e}", file=sys.stderr, flush=True)
        sys.exit(1)


if __name__ == "__main__":
    main()
