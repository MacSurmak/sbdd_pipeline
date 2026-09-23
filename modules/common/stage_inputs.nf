/*
 * modules/common/stage_inputs.nf
 * Safely stages input files into the Nextflow isolated work directory.
 * Prevents Nextflow from operating directly on user-provided files in the project folder.
 */

process STAGE_INPUTS {
    tag "staging"
    executor 'local'
    scratch false
    label 'process_master'

    input:
    path receptor_file
    path library_file
    path ref_file

    output:
    path "staged_receptor.*", emit: staged_receptor, optional: true
    path "staged_library.*",  emit: staged_library,  optional: true
    path "staged_ref.*",      emit: staged_ref,      optional: true

    script:
    def r_cmd = receptor_file.name != 'NO_FILE' ? "cp -L \"${receptor_file}\" \"staged_receptor.${receptor_file.extension}\"" : ""
    def l_cmd = library_file.name != 'NO_FILE' ? "cp -L \"${library_file}\" \"staged_library.${library_file.extension}\"" : ""
    def f_cmd = ref_file.name != 'NO_FILE' ? "cp -L \"${ref_file}\" \"staged_ref.${ref_file.extension}\"" : ""
    """
    $r_cmd
    $l_cmd
    $f_cmd
    """

    stub:
    def r_stub = receptor_file.name != 'NO_FILE' ? "touch \"staged_receptor.${receptor_file.extension}\"" : ""
    def l_stub = library_file.name != 'NO_FILE' ? "touch \"staged_library.${library_file.extension}\"" : ""
    def f_stub = ref_file.name != 'NO_FILE' ? "touch \"staged_ref.${ref_file.extension}\"" : ""
    """
    $r_stub
    $l_stub
    $f_stub
    """
}
