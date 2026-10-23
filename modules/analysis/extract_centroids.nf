/*
 * modules/analysis/extract_centroids.nf
 * Stage 3: Direct seek extraction of 3D medoid structures from trajectories.
 * Preserves chemical topology by coordinate transfer onto SubstanceState templates.
 */

process EXTRACT_CENTROIDS {
    tag "${meta.id}"
    label 'process_single'

    input:
    tuple val(meta), path(ref_pdb), path(xtc_files, stageAs: 'trajs_?/*'), path(medoids_summary), path(template_sdfs), val(heterogens_json)

    output:
    tuple val(meta), path("extracted_manifest.json"), path("centroid_*_receptor.pdb"), path("centroid_*_complex.pdb"), path("centroid_*_receptor.mae"), path("centroid_*_complex.mae"), emit: structures
    tuple val(meta), path("centroid_*_pose.sdf"), path("centroid_*_pose.mae"), emit: poses, optional: true

    script:
    def lig_arg = meta.ligand ? "--ligand_resname ${meta.ligand}" : ""
    def hetero_arg = heterogens_json ? "--heterogens_spec '${heterogens_json}'" : ""
    """
    python3 ${projectDir}/bin/extract_centroids.py \\
        --struct ${ref_pdb} \\
        --xtc ${xtc_files} \\
        --medoids_json ${medoids_summary} \\
        --templates_dir templates \\
        ${hetero_arg} \\
        ${lig_arg} \\
        --outdir .

    # Batch structure conversion in a single Schrodinger runtime call
    schrodinger run python3 ${projectDir}/bin/wrap_schrodinger.py batch_convert --manifest batch_convert_tasks.json
    """
}
