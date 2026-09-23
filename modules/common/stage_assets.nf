/*
 * modules/common/stage_assets.nf
 * Stages local project assets (from Master's /home directory) to the shared NFS scratch
 * so that Slurm workers can access them seamlessly.
 */

process STAGE_ASSETS {
    tag "mdp_templates"
    executor 'local'
    scratch false
    label 'process_master'

    input:
    path mdp_dir

    output:
    path "mdp_templates", emit: mdp_templates

    script:
    """
    mkdir -p mdp_templates
    cp -Lr ${mdp_dir}/* mdp_templates/
    """

    stub:
    """
    mkdir -p mdp_templates
    touch mdp_templates/md.mdp
    """
}
