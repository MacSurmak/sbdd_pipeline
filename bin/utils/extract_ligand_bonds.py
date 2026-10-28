#!/usr/bin/env python3
import sys
from schrodinger import structure
from schrodinger.structutils import assignbondorders


def main():
    if len(sys.argv) < 4:
        sys.exit(1)

    rec_file = sys.argv[1]
    lig_spec = sys.argv[2].strip()
    out_file = sys.argv[3]

    # Support both plain resname ("TGI") and exact chain:resnum:resname ("A:601:NA")
    if ":" in lig_spec:
        parts = lig_spec.split(":")
        req_chain = parts[0].strip()
        req_resnum = int(parts[1].strip())
        lig_name = parts[2].strip().upper()
    else:
        req_chain = None
        req_resnum = None
        lig_name = lig_spec.upper()

    st = structure.StructureReader.read(rec_file)
    lig_atoms = [
        a.index
        for a in st.atom
        if a.pdbres.strip().upper() == lig_name
        and (req_chain is None or a.chain.strip() == req_chain)
        and (req_resnum is None or a.resnum == req_resnum)
    ]

    if not lig_atoms:
        avail = sorted(
            list(
                set(
                    f"{a.chain.strip()}:{a.resnum}:{a.pdbres.strip().upper()}"
                    for a in st.atom
                )
            )
        )
        print(f"RESIDUE_NOT_FOUND: {avail}", file=sys.stderr)
        sys.exit(2)

    lig_st = st.extract(lig_atoms)

    # Strict Fail-Fast: Any extracted pose entity must be a single connected molecular graph
    if len(lig_st.molecule) > 1:
        print(
            f"FATAL: Specification '{lig_spec}' resolves to {len(lig_st.molecule)} disconnected molecules! "
            f"Each cofactor or ion instance must be extracted by its unique residue identifier.",
            file=sys.stderr,
        )
        sys.exit(1)

    lig_st.title = lig_name

    # Assign correct chemical bond orders and formal charges using CCD template when available
    try:
        assignbondorders.assign_st(lig_st, use_ccd=True)
    except Exception as e:
        sys.stderr.write(f"Warning: assignbondorders note: {e}\n")

    lig_st.write(out_file)


if __name__ == "__main__":
    main()
