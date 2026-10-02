/*
 * modules/schrodinger/gridgen.nf
 * Receptor Grid Generation Module for Glide Docking
 */

process GLIDE_GRIDGEN {
    tag "${meta.id}"
    label 'process_single'

    input:
    tuple val(meta), path(receptor_mae), path(ref_lig_mae)
    val(inner_box)
    val(outer_box)
    val(forcefield)

    output:
    tuple val(meta), path("${meta.id}_grid.zip"), emit: grid_zip

    script:
    """
    # 1. Ensure reference ligand is in native MAE format
    if [[ "${ref_lig_mae}" == *.mae || "${ref_lig_mae}" == *.maegz ]]; then
        REF_FILE="${ref_lig_mae}"
    else
        schrodinger structconvert "${ref_lig_mae}" ref_ligand.mae
        REF_FILE="ref_ligand.mae"
    fi

    # 2. Compute exact centroid coordinates directly from staged file
    GRID_CENTER=\$(schrodinger run python3 -c "
import sys
from schrodinger import structure

ref_path = sys.argv[1]
st = structure.StructureReader.read(ref_path)
if st.atom_total == 0:
    sys.exit(1)
cx = sum(a.x for a in st.atom) / st.atom_total
cy = sum(a.y for a in st.atom) / st.atom_total
cz = sum(a.z for a in st.atom) / st.atom_total
print(f'{cx:.4f}, {cy:.4f}, {cz:.4f}')
" "\$REF_FILE")

    if [ -z "\$GRID_CENTER" ]; then
        echo "FATAL: Could not determine GRID_CENTER from \$REF_FILE"
        exit 1
    fi

    # 3. Build input and generate grid
    cat << EOF > gridgen.in
GRIDFILE ${meta.id}_grid.zip
RECEP_FILE ${receptor_mae}
GRID_CENTER \${GRID_CENTER}
REF_LIGAND_FILE \${REF_FILE}
INNERBOX ${inner_box}
OUTERBOX ${outer_box}
FORCEFIELD ${forcefield}
EOF

    schrodinger glide gridgen.in -OVERWRITE -WAIT -HOST localhost
    """

    stub:
    """
    touch ${meta.id}_grid.zip
    """
}
