/*
 * modules/ai/inpaint.nf
 * Monolithic GPU Worker Process for Loop Healing
 * Runs Detect Gaps -> ColabFold -> Graft -> OpenMM -> Restore HETATM in a single container invocation
 */

process HEAL_STRUCTURE {
    tag "${meta.id}"
    label 'process_gpu'

    input:
    tuple val(meta), path(raw_pdb)

    output:
    tuple val(meta), path("final_complex.pdb"), emit: healed_pdb
    path("healing_report.txt"), emit: report

    script:
    """
    echo ">target" > sequence.fasta
    echo "${meta.fasta_seq}" >> sequence.fasta

    # 1. Detect Gaps, clamp construct boundaries and isolate standard PDB components
    colabfold_python ${projectDir}/bin/heal/1_detect_gaps.py \\
        --fasta sequence.fasta \\
        --crystal ${raw_pdb} \\
        --out_gaps gaps.json \\
        --out_protein protein_only.pdb \\
        --out_hetatm hetatm_only.pdb \\
        --out_clamped_fasta clamped_sequence.fasta

    # 2. Prep ColabFold Template from clamped construct
    colabfold_python ${projectDir}/bin/heal/2_prep_cf.py \\
        --protein protein_only.pdb \\
        --fasta clamped_sequence.fasta

    # 3. Predict Loops via ColabFold (100% Offline: single_sequence)
    colabfold_batch \\
        --templates \\
        --custom-template-path custom_templates \\
        --msa-mode single_sequence \\
        --num-models 1 \\
        --num-recycle 3 \\
        query.fasta af2_output

    # 4. Surgical Graft strictly on internal loops
    colabfold_python ${projectDir}/bin/heal/3_graft_loops.py \\
        --crystal protein_only.pdb \\
        --af2_model af2_output/*rank_001*.pdb \\
        --map gaps.json \\
        --out stitched.pdb > healing_report.txt

    # 5. OpenMM Physical Seam Healing strictly on CUDA GPU
    colabfold_python ${projectDir}/bin/heal/4_heal_seams.py \\
        --stitched stitched.pdb \\
        --map gaps.json \\
        --out healed_protein.pdb >> healing_report.txt

    # 6. Clean PDB Complex Assembly
    colabfold_python ${projectDir}/bin/heal/5_restore_complex.py \\
        --healed healed_protein.pdb \\
        --hetatm hetatm_only.pdb \\
        --out final_complex.pdb
    """

    stub:
    """
    touch final_complex.pdb
    touch healing_report.txt
    """
}
