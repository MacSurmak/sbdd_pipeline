/*
 * modules/gromacs/trjconv.nf
 * PBC correction, rotational/translational fitting, and water/ion stripping
 */

process TRJCONV {
    tag "${meta.id}_rep${meta.rep}"
    label 'process_low'

    input:
    tuple val(meta), path(concat_xtc), path(tpr), path(index_ndx)

    output:
    tuple val(meta), path("trajectory_clean.xtc"), path("complex_ref.pdb"), emit: clean_traj

    script:
    """
    # 1. Inject [ COMPLEX ] (Protein + UNK strictly, no lipids) into index.ndx if not present
    python3 -c "
import sys
from pathlib import Path
sys.path.insert(0, '${projectDir}/bin')
from run_production_md import ensure_complex_ndx_group
ensure_complex_ndx_group(Path('${index_ndx}'))
"

    # 2. PBC correction steps: whole -> center -> extract pure COMPLEX
    echo "0" | gmx trjconv -s ${tpr} -f ${concat_xtc} -o pbc_whole.xtc -pbc whole

    echo -e "Protein\\nCOMPLEX" | gmx trjconv -s ${tpr} -f pbc_whole.xtc -n ${index_ndx} -o trajectory_clean.xtc -center -pbc mol -ur compact

    echo "COMPLEX" | gmx trjconv -s ${tpr} -f trajectory_clean.xtc -n ${index_ndx} -o complex_ref.pdb -dump 0

    rm -f pbc_whole.xtc
    """

    stub:
    """
    touch trajectory_clean.xtc
    touch complex_ref.pdb
    """
}

