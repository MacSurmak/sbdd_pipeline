/*
 * modules/membrane/immers.nf
 * PPM 3.0 Standalone Membrane Immersion Calculation
 * Self-contained worker process: zero external script file staging.
 */

process IMMERS {
    tag "${meta.id}"
    label "process_medium"

    input:
    tuple val(meta), path(pure_protein_pdb)
    val(membrane_type)
    val(membrane_topo)
    val(membrane_chains)

    output:
    tuple val(meta), path("oriented.pdb"), emit: oriented_pdb
    tuple val(meta), path("membrane.json"), emit: membrane_json
    tuple val(meta), path("membrane.txt"), emit: membrane_txt

    script:
    """
    python3 ${projectDir}/bin/wrap_membrane.py immers \\
        --pdb ${pure_protein_pdb} \\
        --membrane_type ${membrane_type} \\
        --membrane_topo ${membrane_topo} \\
        --membrane_chains ${membrane_chains} \\
        --out_id .
    mv ._oriented.pdb oriented.pdb
    mv ._membrane.json membrane.json
    mv ._membrane.txt membrane.txt
    """

    stub:
    """
    touch ${meta.id}_oriented.pdb
    echo '{"thickness": 28.0, "tilt": 0.0, "dg": 0.0}' > ${meta.id}_membrane.json
    echo '28.0' > ${meta.id}_membrane.txt
    """
}
