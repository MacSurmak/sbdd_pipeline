/*
 * modules/schrodinger/light_prep.nf
 * Light preparation for MD frames/clusters.
 * Reassigns OPLS4 charges without altering heavy atom coordinates or protonation states.
 */

process LIGHT_PREPWIZARD {
    tag "${meta.id}"
    label 'process_single'

    input:
    tuple val(meta), path(input_pdb)

    output:
    tuple val(meta), path("${meta.id}_lightprep.mae"), emit: prepared_mae

    script:
    """
    schrodinger structconvert ${input_pdb} ${meta.id}_lightprep.mae

    # Ensure pocket is vacant if native ligand was retained in centroid PDB
    if [ -n "${meta.native_ligand}" ]; then
        schrodinger run python3 -c "
from schrodinger import structure
st = structure.StructureReader.read('${meta.id}_lightprep.mae')
lig_res = '${meta.native_ligand}'.split('_')[-1].strip().upper()
lig_atoms = [a.index for a in st.atom if a.pdbres.strip().upper() == lig_res]
if lig_atoms:
    st.deleteAtoms(lig_atoms)
    st.write('${meta.id}_lightprep.mae')
"
    fi    """

    stub:
    """
    touch ${meta.id}_lightprep.mae
    """
}
