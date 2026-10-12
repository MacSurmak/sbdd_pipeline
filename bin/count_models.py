#!/usr/bin/env python3
"""
bin/count_models.py
Counts structures/NMR models in CIF, PDB, MAE, or SDF.
Outputs integer count to stdout and models_count.txt.
"""
import argparse
import sys
from pathlib import Path


def count_cif_models(file_path: Path) -> int:
    models = 0
    with open(file_path, "r", errors="ignore") as f:
        for line in f:
            if line.startswith("_pdbx_nmr_ensemble.conformers_submitted_total_number"):
                parts = line.split()
                if len(parts) >= 2:
                    try:
                        return int(parts[1])
                    except ValueError:
                        pass
            if line.startswith("_pdbx_nmr_ensemble.conformers_calculated_total_number"):
                parts = line.split()
                if len(parts) >= 2:
                    try:
                        return int(parts[1])
                    except ValueError:
                        pass
            if line.strip().startswith("_atom_site.pdbx_PDB_model_num"):
                pass
            if line.strip().startswith("data_") and not line.strip().startswith(
                "data_general"
            ):
                # some CIF files repeat data blocks
                pass
    # If not found via tags, scan model numbers or structure reader
    # Schrodinger fallback if available
    from schrodinger.structure import count_structures

    return count_structures(str(file_path))


def main():
    parser = argparse.ArgumentParser(
        description="Count structures or NMR models in macromolecular file"
    )
    parser.add_argument("--input", required=True, help="Input macromolecule file")
    parser.add_argument("--output", default="models_count.txt", help="Output text file")
    args = parser.parse_args()

    p = Path(args.input)
    count = 1

    if p.suffix.lower() in [".mae", ".maegz"]:
        from schrodinger.structure import count_structures

        count = count_structures(str(p.resolve()))
    elif p.suffix.lower() in [".cif", ".mmcif"]:
        count = count_cif_models(p)
    elif p.suffix.lower() == ".pdb":
        models = 0
        with open(p, "r", errors="ignore") as f:
            for line in f:
                if line.startswith("MODEL"):
                    models += 1
        count = max(1, models)
    else:
        raise ValueError(
            f"FATAL: Unsupported structure format for counting: {p.suffix}"
        )

    print(count)
    with open(args.output, "w") as out:
        out.write(f"{count}\n")


if __name__ == "__main__":
    main()
