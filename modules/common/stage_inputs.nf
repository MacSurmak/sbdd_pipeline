/*
 * modules/common/stage_inputs.nf
 * Stages initial input files from Master local filesystem into /shared scratch storage
 */

process STAGE_INPUTS {
    executor 'local'
    scratch false
    label 'process_single'

    input:
    val rec_file
    val lib_file
    val ref_file

    output:
    path("receptor.*"), emit: staged_receptor
    path("library.sdf"), emit: staged_library, optional: true
    path("ref_native.*"), emit: staged_ref, optional: true

    script:
    """
    ext="\$(basename "${rec_file}" | awk -F. '{if (NF>1) print \$NF; else print "pdb"}')"
    cp -L "${rec_file}" "receptor.\${ext}"

    if [ -n "${lib_file}" ] && [ "${lib_file}" != "NO_FILE" ] && [ -f "${lib_file}" ]; then
        cp -L "${lib_file}" library.sdf
    fi

    if [ -n "${ref_file}" ] && [ "${ref_file}" != "NO_FILE" ] && [ "${ref_file}" != "${rec_file}" ] && [ -f "${ref_file}" ]; then
        ref_ext="\$(basename "${ref_file}" | awk -F. '{if (NF>1) print \$NF; else print "sdf"}')"
        cp -L "${ref_file}" "ref_native.\${ref_ext}"
    fi
    """

    stub:
    """
    touch receptor.pdb
    touch library.sdf
    touch ref_native.sdf
    """
}
