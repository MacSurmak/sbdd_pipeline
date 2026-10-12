#!/usr/bin/env python3
import argparse
import json
import sys
from Bio import Align
from Bio.PDB import PDBParser, PDBIO, Select, is_aa
from Bio.PDB.Polypeptide import protein_letters_3to1


class ProteinOnlySelect(Select):
    def accept_residue(self, res):
        return 1 if is_aa(res) and res.get_id()[0] == " " else 0


def get_pdb_sequence(chain):
    seq = []
    res_ids = []
    for res in chain:
        if is_aa(res) and "CA" in res:
            resname = res.get_resname()
            char = protein_letters_3to1.get(resname, "X")
            seq.append(char)
            res_ids.append(res.get_id())
    return "".join(seq), res_ids


def read_fasta(fasta_file):
    with open(fasta_file, "r") as f:
        seq = []
        for line in f:
            if not line.startswith(">"):
                seq.append(line.strip().upper())
    return "".join(seq)


class HetatmSelect(Select):
    def accept_residue(self, res):
        from Bio.PDB import is_aa

        # The input structure is already 100% sanitized upon import.
        # Retain any non-protein residue (explicitly whitelisted ligand/cofactors).
        return 1 if not (is_aa(res) and res.get_id()[0] == " ") else 0


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--fasta", required=True)
    parser.add_argument("--crystal", required=True)
    parser.add_argument("--out_gaps", default="gaps.json")
    parser.add_argument("--out_protein", default="protein_only.pdb")
    parser.add_argument("--out_hetatm", default="hetatm_only.pdb")
    parser.add_argument("--out_clamped_fasta", default="clamped_sequence.fasta")
    args = parser.parse_args()

    cryst_path = str(args.crystal).lower()
    if cryst_path.endswith((".cif", ".mmcif")):
        from Bio.PDB import MMCIFParser
        parser_struct = MMCIFParser(QUIET=True)
    else:
        from Bio.PDB import PDBParser
        parser_struct = PDBParser(QUIET=True)

    st = parser_struct.get_structure("struct", args.crystal)

    # Auto-detect longest protein chain
    best_chain = None
    best_len = 0
    for chain in st[0]:
        seq, _ = get_pdb_sequence(chain)
        if len(seq) > best_len:
            best_len = len(seq)
            best_chain = chain

    if not best_chain:
        print("FATAL: No amino acids found in crystal structure.", file=sys.stderr)
        sys.exit(1)

    # Write clean standard PDBs
    io = PDBIO()
    io.set_structure(best_chain.get_parent())
    io.save(args.out_protein, ProteinOnlySelect())

    # Isolate HETATMs into clean standard PDB format
    io.set_structure(st[0])
    io.save(args.out_hetatm, HetatmSelect())

    fasta_seq = read_fasta(args.fasta)
    pdb_seq, pdb_res_ids = get_pdb_sequence(best_chain)

    aligner = Align.PairwiseAligner()
    aligner.mode = "global"
    aligner.match_score = 2
    aligner.mismatch_score = -1
    aligner.open_gap_score = -5
    aligner.extend_gap_score = -0.5

    alignments = aligner.align(fasta_seq, pdb_seq)
    best_alignment = alignments[0]

    blocks_fasta = best_alignment.aligned[0]
    blocks_pdb = best_alignment.aligned[1]

    raw_pairs = []
    for (f_start, f_end), (p_start, p_end) in zip(blocks_fasta, blocks_pdb):
        for offset in range(f_end - f_start):
            f_idx = f_start + offset + 1  # 1-based original FASTA index
            p_idx = p_start + offset
            pid = pdb_res_ids[p_idx]
            raw_pairs.append(
                {
                    "fasta_resSeq": int(f_idx),
                    "crystal_seq": int(pid[1]),
                    "crystal_icode": str(pid[2]),
                }
            )

    aligned_indices = sorted([p["fasta_resSeq"] for p in raw_pairs])
    first_resolved = aligned_indices[0]
    last_resolved = aligned_indices[-1]

    # Strictly clamp sequence to resolved construct (discard floppy terminal overhangs)
    clamped_seq = fasta_seq[first_resolved - 1 : last_resolved]
    clamped_offset = first_resolved - 1

    aligned_set = set(aligned_indices)
    internal_gaps = []
    in_gap = False
    gap_start = None

    # Detect strictly internal missing loops
    for i in range(first_resolved, last_resolved + 1):
        if i not in aligned_set:
            if not in_gap:
                in_gap = True
                gap_start = i
        else:
            if in_gap:
                internal_gaps.append(
                    {
                        "af2_start": int(gap_start - clamped_offset),
                        "af2_end": int((i - 1) - clamped_offset),
                        "orig_start": int(gap_start),
                        "orig_end": int(i - 1),
                        "sequence": str(fasta_seq[gap_start - 1 : i - 1]),
                    }
                )
                in_gap = False

    clamped_superposition_pairs = []
    for p in raw_pairs:
        clamped_superposition_pairs.append(
            {
                "af2_resSeq": int(p["fasta_resSeq"] - clamped_offset),
                "crystal_seq": int(p["crystal_seq"]),
                "crystal_icode": str(p["crystal_icode"]),
            }
        )

    output_data = {
        "fasta_length": len(clamped_seq),
        "chain_id": best_chain.id,
        "first_resolved": first_resolved,
        "last_resolved": last_resolved,
        "internal_gaps": internal_gaps,
        "superposition_pairs": clamped_superposition_pairs,
    }

    with open(args.out_gaps, "w", encoding="utf-8") as f:
        json.dump(output_data, f, indent=2)

    with open(args.out_clamped_fasta, "w", encoding="utf-8") as f:
        f.write(f">clamped\n{clamped_seq}\n")


if __name__ == "__main__":
    main()
