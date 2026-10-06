/*
 * modules/schrodinger/glide_dock.nf
 * Glide Docking Module
 */

process GLIDE_DOCK {
    tag "${meta.id}"
    label 'process_medium'

    input:
    tuple val(meta), path(grid_zip), path(ligands_mae)
    val(precision)
    val(poses_per_lig)
    val(forcefield)
    val(lig_vscale)
    val(lig_ccut)

    output:
    tuple val(meta), path("glide_dock_pv.maegz"), emit: pv_maegz
    tuple val(meta), path("glide_dock.csv"), emit: scores_csv
    tuple val(meta), path("glide_dock.sdf"), emit: poses_sdf

    script:
    """
    schrodinger run python3 ${projectDir}/bin/wrap_schrodinger.py glide_dock \\
        --grid_zip ${grid_zip} \\
        --ligands_mae ${ligands_mae} \\
        --precision ${precision} \\
        --poses_per_lig ${poses_per_lig} \\
        --forcefield ${forcefield} \\
        --lig_vscale ${lig_vscale} \\
        --lig_ccut ${lig_ccut} \\
        --out_pv glide_dock_pv.maegz \\
        --out_csv glide_dock.csv \\
        --out_sdf glide_dock.sdf
    """

    stub:
    """
    touch ${meta.id}_dock_pv.maegz
    touch ${meta.id}_dock.csv
    touch ${meta.id}_dock.sdf
    """
}

