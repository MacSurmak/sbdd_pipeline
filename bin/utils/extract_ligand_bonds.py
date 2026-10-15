#!/usr/bin/env python3
import sys
from schrodinger import structure
from schrodinger.structutils import assignbondorders


def main():
    if len(sys.argv) < 4:
        sys.exit(1)

    rec_file = sys.argv[1]
    lig_name = sys.argv[2].strip().upper()
    out_file = sys.argv[3]

    st = structure.StructureReader.read(rec_file)
    lig_atoms = [a.index for a in st.atom if a.pdbres.strip().upper() == lig_name]

    if not lig_atoms:
        avail = sorted(list(set(a.pdbres.strip().upper() for a in st.atom)))
        print(f"RESIDUE_NOT_FOUND: {avail}", file=sys.stderr)
        sys.exit(2)

    lig_st = st.extract(lig_atoms)
    
    # Strict Fail-Fast: Organic molecules must be a single connected graph
    if len(lig_st.molecule) > 1 and len(lig_atoms) > len(lig_st.molecule):
        print(
            f"FATAL: Residue name '{lig_name}' resolves to {len(lig_st.molecule)} disconnected molecules in the pocket! "
            f"This destroys quantum parametrization. You must assign unique residue names (e.g. ATP, AT2) to distinct molecules.",
            file=sys.stderr
        )
        sys.exit(1)

    lig_st.title = lig_name

    # Assign correct chemical bond orders and formal charges from 3D geometry
    try:
        assignbondorders.assign_st(lig_st)
    except Exception as e:
        sys.stderr.write(f"Warning: assignbondorders note: {e}\n")

    lig_st.write(out_file)


if __name__ == "__main__":
    main()
