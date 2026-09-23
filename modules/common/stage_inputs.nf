/*
 * modules/common/stage_inputs.nf
 * Bridges the gap between non-shared user directories (e.g. /home, /tmp) and the Slurm NFS scratch.
 * Runs locally on the Master node to copy input files into the /shared workdir.
 */

process STAGE_INPUTS {
    tag "master_to_nfs"
    executor 'local'
    scratch false
    label 'process_master'

    input:
    val rec_path
    val lib_path
    val dock_path

    output:
    path "staged_receptor.*", emit: staged_receptor, optional: true
    path "staged_library.*",  emit: staged_library,  optional: true
    path "staged_docked.*",   emit: staged_docked,   optional: true

    script:
    def r_ext = rec_path != 'NO_FILE' ? file(rec_path).extension : ""
    def l_ext = lib_path != 'NO_FILE' ? file(lib_path).extension : ""
    def d_ext = dock_path != 'NO_FILE' ? file(dock_path).extension : ""

    def r_cmd = rec_path != 'NO_FILE' ? "cp -L \"${rec_path}\" \"staged_receptor.${r_ext}\"" : ""
    def l_cmd = lib_path != 'NO_FILE' ? "cp -L \"${lib_path}\" \"staged_library.${l_ext}\"" : ""
    def d_cmd = dock_path != 'NO_FILE' ? "cp -L \"${dock_path}\" \"staged_docked.${d_ext}\"" : ""
    """
    $r_cmd
    $l_cmd
    $d_cmd
    """

    stub:
    def r_ext = rec_path != 'NO_FILE' ? file(rec_path).extension : ""
    def l_ext = lib_path != 'NO_FILE' ? file(lib_path).extension : ""
    def d_ext = dock_path != 'NO_FILE' ? file(dock_path).extension : ""

    def r_stub = rec_path != 'NO_FILE' ? "touch \"staged_receptor.${r_ext}\"" : ""
    def l_stub = lib_path != 'NO_FILE' ? "touch \"staged_library.${l_ext}\"" : ""
    def d_stub = dock_path != 'NO_FILE' ? "touch \"staged_docked.${d_ext}\"" : ""
    """
    $r_stub
    $l_stub
    $d_stub
    """
}
