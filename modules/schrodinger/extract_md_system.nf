/*
 * modules/schrodinger/extract_md_system.nf
 * Extracts dry receptor, cofactors, and ligands for MD system setup
 */

process EXTRACT_MD_SYSTEM {
    tag "${meta.id}"
    label 'process_single'

    input:
    tuple val(meta), path(rec_file), path(lig_file), val(cofactors)

    output:
    tuple val(meta), path("receptor_dry.pdb"), path("organics.json"), path("org_files/*"), emit: extracted, optional: true
    tuple val(meta), path("receptor_dry.pdb"), path("organics.json"), emit: extracted_base

    script:
    def cof_list = cofactors instanceof List ? cofactors.join(' ') : cofactors
    def cof_arg = cof_list ? "--cofactors ${cof_list}" : ""
    def lig_arg = ""
    if (lig_file) {
        if (!meta.native_ligand) {
            error "FATAL: meta.native_ligand is required when lig_file is provided for system '${meta.id}'"
        }
        lig_arg = "--lig_file ${lig_file} --lig_name ${meta.native_ligand}"
    }
    """
    schrodinger run ${projectDir}/bin/extract_md_files.py \\
        --rec_file ${rec_file} \\
        ${lig_arg} \\
        ${cof_arg} \\
        --out_dir .
    """

    stub:
    """
    touch receptor_dry.pdb
    echo "[]" > organics.json
    """
}
