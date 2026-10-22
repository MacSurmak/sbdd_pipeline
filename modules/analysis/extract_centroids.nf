/*
 * modules/analysis/extract_centroids.nf
 * Stage 3: Direct seek extraction of 3D medoid structures from trajectories.
 * Preserves chemical topology by coordinate transfer onto SubstanceState templates.
 */

process EXTRACT_CENTROIDS {
    tag "${meta.id}"
    label 'process_single'

    input:
    tuple val(meta), path(ref_pdb), path(xtc_files, stageAs: 'trajs_?/*'), path(medoids_summary), path(template_sdf)

    output:
    tuple val(meta), path("centroid_*_receptor.pdb"), path("centroid_*_complex.pdb"), path("centroid_*_receptor.mae"), path("extracted_manifest.json"), emit: structures
    tuple val(meta), path("centroid_*_pose.sdf"), path("centroid_*_pose.mae"), emit: poses, optional: true

    script:
    def tmpl_arg = (template_sdf && template_sdf.name != 'NO_FILE') ? "--substance_template_sdf ${template_sdf}" : ""
    def lig_arg = meta.ligand ? "--ligand_resname ${meta.ligand}" : ""
    """
    python3 ${projectDir}/bin/extract_centroids.py \\
        --struct ${ref_pdb} \\
        --xtc ${xtc_files} \\
        --medoids_json ${medoids_summary} \\
        ${tmpl_arg} \\
        ${lig_arg} \\
        --outdir .

    # Standardized conversion via containerized Schrödinger wrapper in PATH
    for rec in centroid_*_receptor.pdb; do
        mae="\${rec%.pdb}.mae"
        schrodinger run structconvert "\${rec}" "\${mae}"
    done

    for sdf in centroid_*_pose.sdf; do
        if [ -f "\${sdf}" ]; then
            mae="\${sdf%.sdf}.mae"
            schrodinger run structconvert "\${sdf}" "\${mae}"
        fi
    done
    """
}
