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
    path script_file

    output:
    tuple val(meta), path("${meta.id}_dock_pv.maegz"), emit: pv_maegz
    tuple val(meta), path("${meta.id}_dock.csv"), emit: scores_csv
    tuple val(meta), path("${meta.id}_dock.sdf"), emit: poses_sdf

    script:
    """
    cat << EOF > dock.in
GRIDFILE        ${grid_zip}
LIGANDFILE      ${ligands_mae}
PRECISION       ${precision}
POSES_PER_LIG   ${poses_per_lig}
POSTDOCK        yes
EPIK_PENALTIES  yes
SAMPLE_RINGS    yes
WRITE_CSV       yes
FORCEFIELD      ${forcefield}
EOF

    schrodinger glide dock.in -OVERWRITE -WAIT
    mv dock_pv.maegz ${meta.id}_dock_pv.maegz
    mv dock.csv ${meta.id}_dock.csv

    schrodinger run ${script_file} --pv ${meta.id}_dock_pv.maegz --sdf ${meta.id}_dock.sdf
    """

    stub:
    """
    touch ${meta.id}_dock_pv.maegz
    touch ${meta.id}_dock.csv
    touch ${meta.id}_dock.sdf
    """
}

