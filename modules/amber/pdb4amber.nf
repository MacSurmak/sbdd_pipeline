/*
 * modules/amber/pdb4amber.nf
 * Official AmberTools utility to standardize PDB files for tLEaP
 * Resolves terminal hydrogens, protonation states, and extracts disulfide bonds.
 */

process PDB4AMBER {
    tag "${meta.id}"
    label 'process_single'

    input:
    tuple val(meta), path(receptor_pdb)

    output:
    tuple val(meta), path("rec_amber.pdb"), path("rec_amber_sslink"), emit: amber_rec

    script:
    """
    # Let pdb4amber properly handle all naming conventions, terminal caps, and SS-bonds natively
    pdb4amber -i ${receptor_pdb} -o rec_amber.pdb
    """

    stub:
    """
    touch rec_amber.pdb
    touch rec_amber_sslink
    """
}
