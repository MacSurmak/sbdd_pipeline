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
    val(water_names)
    val(salt_conc)
    val(protein_ff)
    val(gaff_ff)
    val(lipid_ff)

    output:
    tuple val(meta), path("complex.prmtop"), path("complex.inpcrd"), emit: amber_system
    tuple val(meta), path("complex_tleap.pdb"), emit: tleap_pdb

    script:
    def mem_arg = (is_membrane.toString() == 'true') ? "--membrane" : ""
    def mol2_arg = (mol2s && mol2s.size() > 0) ? "--mol2s ${mol2s}" : ""
    def frcmod_arg = (frcmods && frcmods.size() > 0) ? "--frcmods ${frcmods}" : ""
    """
    python3 ${projectDir}/bin/wrap_amber.py tleap \\
        --complex_pdb ${complex_pdb} \\
        --water_leaprc ${water_leaprc} \\
        --water_box ${water_box} \\
        --box_type ${box_type} \\
        --box_padding ${box_padding} \\
        --ss_cutoff ${ss_cutoff} \\
        --water_names "${water_names}" \\
        --salt_conc ${salt_conc} \\
        --protein_ff ${protein_ff} \\
        --gaff_ff ${gaff_ff} \\
        --lipid_ff ${lipid_ff} \\
        ${mol2_arg} \\
        ${frcmod_arg} \\
        ${mem_arg}
    """
}

process TLEAP_DRY {
    tag "${meta.id}"
    label 'process_single'

    input:
    tuple val(meta), path(complex_dry_pdb), path(mol2s), path(frcmods)
    val(water_leaprc)
    val(ss_cutoff)
    val(protein_ff)
    val(gaff_ff)
    val(lipid_ff)

    output:
    tuple val(meta), path("complex_dry.prmtop"), path("complex_dry.inpcrd"), emit: dry_amber

    script:
    def mol2_arg = (mol2s && mol2s.size() > 0) ? "--mol2s ${mol2s}" : ""
    def frcmod_arg = (frcmods && frcmods.size() > 0) ? "--frcmods ${frcmods}" : ""
    """
    python3 ${projectDir}/bin/wrap_amber.py tleap_dry \\
        --complex_pdb ${complex_dry_pdb} \\
        --water_leaprc ${water_leaprc} \\
        --ss_cutoff ${ss_cutoff} \\
        --protein_ff ${protein_ff} \\
        --gaff_ff ${gaff_ff} \\
        --lipid_ff ${lipid_ff} \\
        ${mol2_arg} \\
        ${frcmod_arg}
    """
}
