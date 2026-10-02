/*
 * modules/schrodinger/gridgen.nf
 * Receptor Grid Generation Module for Glide Docking
 */

process GLIDE_GRIDGEN {
    tag "${meta.id}"
    label 'process_single'

    input:
    tuple val(meta), path(receptor_mae, stageAs: 'rec_dir/*'), path(ref_lig_mae, stageAs: 'ref_dir/*')
    val(inner_box)
    val(outer_box)
    val(forcefield)

    output:
    tuple val(meta), path("${meta.id}_grid.zip"), emit: grid_zip

    script:
    """
    REF_BASE=\$(basename "${ref_lig_mae}")
    REF_MAE="ref_ligand.mae"
    if [[ "\$REF_BASE" == *.mae || "\$REF_BASE" == *.maegz ]]; then
        REF_MAE="${ref_lig_mae}"
    elif [[ "\$REF_BASE" != "NO_FILE" && -f "${ref_lig_mae}" ]]; then
        schrodinger structconvert "${ref_lig_mae}" ref_ligand.mae
        REF_MAE="ref_ligand.mae"
    else
        REF_MAE="${ref_lig_mae}"
    fi

    # Compute exact centroid coordinates of reference pose for Glide GRID_CENTER (Fail-Fast)
    GRID_CENTER=$(schrodinger run python3 -c "
import sys
from schrodinger import structure
try:
    st = structure.StructureReader.read('${REF_MAE}')
    if st.atom_total > 0:
        cx = sum(a.x for a in st.atom) / st.atom_total
        cy = sum(a.y for a in st.atom) / st.atom_total
        cz = sum(a.z for a in st.atom) / st.atom_total
        print(f'{cx:.4f}, {cy:.4f}, {cz:.4f}')
    else:
        print('FATAL: Reference ligand structure is empty.', file=sys.stderr)
        sys.exit(1)
except Exception as e:
    print(f'FATAL: Failed to calculate grid center from ${REF_MAE}: {e}', file=sys.stderr)
    sys.exit(1)
")
    export SCHRODINGER_JOBDB="/tmp/schrodinger_jobdb_\${SLURM_JOB_ID:-\$\$}"
    mkdir -p "\$SCHRODINGER_JOBDB"

    cat << EOF > gridgen.in
GRIDFILE ${meta.id}_grid.zip
RECEP_FILE ${receptor_mae}
GRID_CENTER \${GRID_CENTER}
REF_LIGAND_FILE \${REF_MAE}
INNERBOX ${inner_box}
OUTERBOX ${outer_box}
FORCEFIELD ${forcefield}
EOF

    schrodinger glide gridgen.in -OVERWRITE -WAIT -HOST localhost    """

    stub:
    """
    touch ${meta.id}_grid.zip
    """
}

