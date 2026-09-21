/*
 * modules/common/stage_inputs.nf
 * Stages initial input files from Master local filesystem into /shared scratch storage
 */

process STAGE_INPUTS {
    executor 'local'
    scratch false
    label 'process_single'

    input:
    path rec_file, stageAs: 'in_rec/*'
    path lib_file, stageAs: 'in_lib/*'
    path ref_file, stageAs: 'in_ref/*'

    output:
    path("receptor.*"), emit: staged_receptor
    path("library.sdf"), emit: staged_library, optional: true
    path("ref_native.*"), emit: staged_ref, optional: true

    script:
    def lib_copy = (lib_file && lib_file.name != 'NO_FILE' && lib_file.name != '.nobackup') ? "cp -L ${lib_file} library.sdf" : ""
    def ref_copy = (ref_file && ref_file.name != 'NO_FILE' && ref_file.name != '.nobackup' && ref_file.name != rec_file.name) ? "cp -L ${ref_file} ref_native.${ref_file.extension}" : ""
    """
    cp -L ${rec_file} receptor.${rec_file.extension}
    ${lib_copy}
    ${ref_copy}
    """

    stub:
    """
    touch receptor.${rec_file.extension}
    touch library.sdf
    touch ref_native.sdf
    """
}

