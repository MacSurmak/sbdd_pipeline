/*
 * modules/schrodinger/ligprep.nf
 * Ligand Preparation with Epik and ForceField
 */

process LIGPREP {
    tag "${meta.id}"
    label 'process_medium'

    input:
    tuple val(meta), path(ligands_sdf)
    val(target_ph)
    val(ph_tolerance)
    val(ff_version)
    val(max_states)

    output:
    tuple val(meta), path("${meta.id}_prepared.maegz"), emit: prepared_maegz
    tuple val(meta), path("${meta.id}_prepared.sdf"), emit: prepared_sdf

    script:
    """
    cat << 'EOF' > ligprep.in
INPUT_FILE_NAME   ${ligands_sdf}
OUT_MAE           ${meta.id}_prepared.maegz
FORCE_FIELD       ${ff_version}
IONIZATION        epik
PH                ${target_ph}
PH_THRESHOLD      ${ph_tolerance}
STEREO_CHEMISTRY  retain
MAX_OUTPUT_STEREO ${max_states}
EOF

    schrodinger ligprep -inp ligprep.in -WAIT
    schrodinger structconvert ${meta.id}_prepared.maegz ${meta.id}_prepared.sdf
    """
}

