/*
 * modules/membrane/declash.nf
 * Relieves steric overlaps in lipid/water/protein systems
 * Self-contained worker process: zero external script file staging.
 */

process DECLASH {
    tag "${meta.id}"
    label 'process_medium'

    input:
    tuple val(meta), path(input_gro), path(topol_top), path(dry_prmtop)
    val(cutoff)
    val(niter)
    val(damping)
    val(max_move)
    val(jitter)
    val(tol)

    output:
    tuple val(meta), path("${meta.id}_declashed.gro"), emit: declashed_gro

    script:
    def fix_lig = (meta.lig_resname ?: meta.native_ligand) ? "--fix-moltype ${meta.lig_resname ?: meta.native_ligand}" : ""
    """
    python3 ${projectDir}/bin/declash_gmx.py \\
        -f ${input_gro} \\
        -p ${topol_top} \\
        --dry-top ${dry_prmtop} \\
        -o ${meta.id}_declashed.gro \\
        --fix-moltype system1 ${fix_lig} \\
        --no-heavy-only \\
        --cutoff ${cutoff} \\
        --niter ${niter} \\
        --damping ${damping} \\
        --max-move ${max_move} \\
        --jitter ${jitter} \\
        --tolerance ${tol} \\
        --report declash_report.txt
    """
}
