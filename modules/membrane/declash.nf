/*
 * modules/membrane/declash.nf
 * Relieves steric overlaps in lipid/water/protein systems
 * Self-contained worker process: zero external script file staging.
 */

process DECLASH {
    tag "${meta.id}"
    label 'process_medium'

    input:
    tuple val(meta), path(input_gro), path(topol_top)
    val(cutoff)
    val(niter)
    val(damping)
    val(max_move)
    val(jitter)
    val(tol)

    output:
    tuple val(meta), path("${meta.id}_declashed.gro"), emit: declashed_gro

    script:
    """
    python3 ${projectDir}/bin/declash_gmx.py \\
        -f ${input_gro} \\
        -p ${topol_top} \\
        -o ${meta.id}_declashed.gro \\
        --fix-moltype system1 \\
        ${meta.native_ligand ? "--fix-moltype ${meta.native_ligand} \\\\" : ""}
        --no-heavy-only \\
        --cutoff ${cutoff} \\
        --niter ${niter} \\
        --damping ${damping} \\
        --max-move ${max_move} \\
        --jitter ${jitter} \\
        --tolerance ${tol} \\
        --report declash_report.txt
    """

    stub:
    """
    touch ${meta.id}_declashed.gro
    """
}
