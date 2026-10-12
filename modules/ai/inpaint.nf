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

    # 1. Detect Gaps and extract pure protein
    colabfold_python ${projectDir}/bin/heal/1_detect_gaps.py \\
        --fasta sequence.fasta \\
        --crystal ${raw_pdb} \\
        --out_gaps gaps.json \\
        --out_protein protein_only.pdb

    # 2. Prep ColabFold Template (blank HETATM and correct numbering)
    colabfold_python ${projectDir}/bin/heal/2_prep_cf.py \\
        --protein protein_only.pdb \\
        --fasta sequence.fasta

    # 3. Predict Loops via ColabFold
    colabfold_batch \\
        --templates \\
        --custom-template-path custom_templates \\
        --num-models 1 \\
        --num-recycle 3 \\
        query.fasta af2_output

    # 4. Surgical Graft
    colabfold_python ${projectDir}/bin/heal/3_graft_loops.py \\
        --crystal protein_only.pdb \\
        --af2_model af2_output/*rank_001*.pdb \\
        --map gaps.json \\
        --out stitched.pdb > healing_report.txt

    # 5. OpenMM Physical Seam Healing
    colabfold_python ${projectDir}/bin/heal/4_heal_seams.py \\
        --stitched stitched.pdb \\
        --map gaps.json \\
        --out healed_protein.pdb >> healing_report.txt

    # 6. Restore Complex (Re-inject HETATM)
    colabfold_python ${projectDir}/bin/heal/5_restore_complex.py \\
        --healed healed_protein.pdb \\
        --original ${raw_pdb} \\
        --out final_complex.pdb
    """

    stub:
    """
    touch final_complex.pdb
    touch healing_report.txt
    """
}
