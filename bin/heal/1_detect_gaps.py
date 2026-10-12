#!/usr/bin/env python3
import argparse
import json
import sys
from Bio import Align
from Bio.PDB import PDBParser, PDBIO, Select, is_aa
from Bio.PDB.Polypeptide import protein_letters_3to1

class ProteinOnlySelect(Select):
    def accept_residue(self, res):
        return 1 if is_aa(res) and res.get_id()[0] == ' ' else 0

def get_pdb_sequence(chain):
    seq = []
    res_ids = []
    for res in chain:
        if is_aa(res) and 'CA' in res:
            resname = res.get_resname()
            char = protein_letters_3to1.get(resname, 'X')
            seq.append(char)
            res_ids.append(res.get_id())
    return "".join(seq), res_ids

def read_fasta(fasta_file):
    with open(fasta_file, 'r') as f:
        seq = []
        for line in f:
            if not line.startswith(">"):
                seq.append(line.strip().upper())
    return "".join(seq)

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--fasta", required=True)
    parser.add_argument("--crystal", required=True)
    parser.add_argument("--out_gaps", default="gaps.json")
    parser.add_argument("--out_protein", default="protein_only.pdb")
    args = parser.parse_args()

    parser_pdb = PDBParser(QUIET=True)
    st = parser_pdb.get_structure("struct", args.crystal)
    
    # Auto-detect longest chain
    best_chain = None
    best_len = 0
    for chain in st[0]:
        seq, _ = get_pdb_sequence(chain)
        if len(seq) > best_len:
            best_len = len(seq)
            best_chain = chain

    if not best_chain:
        print("FATAL: No amino acids found in PDB.", file=sys.stderr)
        sys.exit(1)

    # Write pure protein PDB to isolate from HETATM
    io = PDBIO()
    io.set_structure(best_chain.get_parent())
    io.save(args.out_protein, ProteinOnlySelect())

    fasta_seq = read_fasta(args.fasta)
    pdb_seq, pdb_res_ids = get_pdb_sequence(best_chain)

    aligner = Align.PairwiseAligner()
    aligner.mode = 'global'
    aligner.match_score = 2
    aligner.mismatch_score = -1
    aligner.open_gap_score = -5
    aligner.extend_gap_score = -0.5

    alignments = aligner.align(fasta_seq, pdb_seq)
    best_alignment = alignments[0]

    blocks_fasta = best_alignment.aligned[0]
    blocks_pdb = best_alignment.aligned[1]

    internal_gaps = []
    superposition_pairs = []

    for (f_start, f_end), (p_start, p_end) in zip(blocks_fasta, blocks_pdb):
        for offset in range(f_end - f_start):
            f_idx = f_start + offset + 1
            p_idx = p_start + offset
            pid = pdb_res_ids[p_idx]
            superposition_pairs.append({
                "af2_resSeq": int(f_idx),
                "crystal_seq": int(pid[1]),
                "crystal_icode": str(pid[2])
            })

    aligned_f_indices = set(p["af2_resSeq"] for p in superposition_pairs)
    in_gap = False
    gap_start = None

    for i in range(1, len(fasta_seq) + 1):
        if i not in aligned_f_indices:
            if not in_gap:
                in_gap = True
                gap_start = i
        else:
            if in_gap:
                is_n_term = (gap_start == 1)
                internal_gaps.append({
                    "af2_start": int(gap_start),
                    "af2_end": int(i - 1),
                    "is_terminal": is_n_term,
                    "sequence": str(fasta_seq[gap_start-1 : i-1])
                })
                in_gap = False

    if in_gap:
        internal_gaps.append({
            "af2_start": int(gap_start),
            "af2_end": int(len(fasta_seq)),
            "is_terminal": True,
            "sequence": str(fasta_seq[gap_start-1 : len(fasta_seq)])
        })

    output_data = {
        "fasta_length": len(fasta_seq),
        "chain_id": best_chain.id,
        "internal_gaps": internal_gaps,
        "superposition_pairs": superposition_pairs
    }

    with open(args.out_gaps, 'w') as f:
        json.dump(output_data, f, indent=2)

if __name__ == "__main__":
    main()
