/*
 * modules/reporting/select_hits.nf
 * Selects top candidate hits from ranking.csv for downstream MD
 */

process SELECT_HITS {
    tag "${meta.id}"
    executor 'local'
    scratch false
    label 'process_single'

    input:
    tuple val(meta), path(ranking_csv)
    val(top_n)
    val(include_native)
    val(native_ligand)

    output:
    tuple val(meta), path("selected_hits.txt"), emit: hits_file
    path("selected_hits.txt")

    script:
    """
    python3 ${projectDir}/bin/select_hits.py \\
    script:
    def nat_arg = native_ligand ? "--native_ligand ${native_ligand}" : ""
    """
    python3 ${projectDir}/bin/select_hits.py \
        --ranking_csv ${ranking_csv} \
        --top_n ${top_n} \
        --include_native ${include_native} \
        ${nat_arg} \
        --out selected_hits.txt
    """

    stub:
    """
    echo "PPZ_1H_Ph_BnMe,1" > selected_hits.txt
    echo "${native_ligand ?: 'UNK'},1" >> selected_hits.txt
    """