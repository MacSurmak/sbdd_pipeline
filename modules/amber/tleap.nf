/*
 * modules/amber/tleap.nf
 * Amber tLEaP System Topology and Coordinate Builder
 */

process TLEAP {
    tag "${meta.id}"
    label 'process_single'

    input:
    tuple val(meta), path(complex_pdb), path(mol2s), path(frcmods)
    val(water_leaprc)
    val(water_box)
    val(box_type)
    val(is_membrane)
    val(box_padding)
    val(ss_cutoff)

    output:
    tuple val(meta), path("complex.prmtop"), path("complex.inpcrd"), emit: amber_system
    tuple val(meta), path("complex_tleap.pdb"), emit: tleap_pdb

    script:
    def mem_arg = (is_membrane.toString() == 'true') ? "--membrane" : ""
    """
    python3 ${projectDir}/bin/wrap_amber.py tleap \\
        --complex_pdb ${complex_pdb} \\
        --water_leaprc ${water_leaprc} \\
        --water_box ${water_box} \\
        --box_type ${box_type} \\
        --box_padding ${box_padding} \\
        --ss_cutoff ${ss_cutoff} \\
        --mol2s ${mol2s} \\
        --frcmods ${frcmods} \\
        ${mem_arg}
    """

    stub:
    """
    touch complex.prmtop
    touch complex.inpcrd
    touch complex_tleap.pdb
    """
}
