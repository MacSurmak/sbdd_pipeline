#!/usr/bin/env python3
import argparse
import json
import sys
import numpy as np
from Bio.PDB import PDBParser, PDBIO, Superimposer, Structure, Model, Chain


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--crystal", required=True)
    parser.add_argument("--af2_model", required=True)
    parser.add_argument("--map", required=True)
    parser.add_argument("--out", required=True)
    args = parser.parse_args()

    with open(args.map, "r") as f:
        mapping = json.load(f)

    chain_id = mapping["chain_id"]
    parser_pdb = PDBParser(QUIET=True)

    crystal_st = parser_pdb.get_structure("crystal", args.crystal)
    af2_st = parser_pdb.get_structure("af2", args.af2_model)

    c_chain = list(crystal_st.get_chains())[0]  # Single protein chain
    a_chain = list(af2_st.get_chains())[0]

    af2_to_cryst = {
        p["af2_resSeq"]: (" ", p["crystal_seq"], p["crystal_icode"])
        for p in mapping["superposition_pairs"]
    }

    # Global Baseline
    c_atoms_global = []
    a_atoms_global = []
    for pair in mapping["superposition_pairs"]:
        c_res_id = (" ", pair["crystal_seq"], pair["crystal_icode"])
        a_res_id = (" ", pair["af2_resSeq"], " ")
        if c_res_id in c_chain and a_res_id in a_chain:
            if "CA" in c_chain[c_res_id] and "CA" in a_chain[a_res_id]:
                c_atoms_global.append(c_chain[c_res_id]["CA"])
                a_atoms_global.append(a_chain[a_res_id]["CA"])

    global_imposer = Superimposer()
    global_imposer.set_atoms(c_atoms_global, a_atoms_global)
    global_imposer.apply(af2_st.get_atoms())
    print(f"Global Baseline Core CA RMSD: {global_imposer.rms:.4f} A")

    loop_transformed_residues = {}

    for g in mapping["internal_gaps"]:
        start = g["af2_start"]
        end = g["af2_end"]
        is_terminal = g.get("is_terminal", False)

        n_stems = [
            start - offset for offset in range(1, 5) if (start - offset) in af2_to_cryst
        ]
        c_stems = [
            end + offset for offset in range(1, 5) if (end + offset) in af2_to_cryst
        ]

        has_two_sided = (len(n_stems) >= 2 and len(c_stems) >= 2) and not is_terminal

        c_stem_atoms = []
        a_stem_atoms = []

        if has_two_sided:
            for idx in n_stems + c_stems:
                c_res_id = af2_to_cryst[idx]
                a_res_id = (" ", idx, " ")
                if c_res_id in c_chain and a_res_id in a_chain:
                    for aname in ["N", "CA", "C"]:
                        if aname in c_chain[c_res_id] and aname in a_chain[a_res_id]:
                            c_stem_atoms.append(c_chain[c_res_id][aname])
                            a_stem_atoms.append(a_chain[a_res_id][aname])

        if has_two_sided and len(c_stem_atoms) >= 12:
            local_imposer = Superimposer()
            local_imposer.set_atoms(c_stem_atoms, a_stem_atoms)

            for r_idx in range(start, end + 1):
                a_res_id = (" ", r_idx, " ")
                if a_res_id in a_chain:
                    res_copy = a_chain[a_res_id].copy()
                    for atom in res_copy:
                        atom.coord = (
                            np.dot(atom.coord, local_imposer.rotran[0])
                            + local_imposer.rotran[1]
                        )
                    loop_transformed_residues[r_idx] = res_copy
            print(
                f"Loop [{start}-{end}]: Grafted with local stem RMSD {local_imposer.rms:.4f} A"
            )
        else:
            for r_idx in range(start, end + 1):
                a_res_id = (" ", r_idx, " ")
                if a_res_id in a_chain:
                    loop_transformed_residues[r_idx] = a_chain[a_res_id].copy()
            print(f"Region [{start}-{end}]: Terminal overhang preserved globally.")

    hybrid_st = Structure.Structure("hybrid")
    hybrid_model = Model.Model(0)
    hybrid_chain = Chain.Chain("A")
    hybrid_model.add(hybrid_chain)
    hybrid_st.add(hybrid_model)

    for i in range(1, mapping["fasta_length"] + 1):
        if i in loop_transformed_residues:
            res = loop_transformed_residues[i]
            res.id = (" ", i, " ")
            hybrid_chain.add(res)
        elif i in af2_to_cryst:
            c_res_id = af2_to_cryst[i]
            if c_res_id in c_chain:
                res = c_chain[c_res_id].copy()
                res.id = (" ", i, " ")
                hybrid_chain.add(res)

    io = PDBIO()
    io.set_structure(hybrid_st)
    io.save(args.out)


if __name__ == "__main__":
    main()
