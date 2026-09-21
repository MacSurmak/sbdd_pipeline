/*
 * modules/gromacs/mmpbsa.nf
 * Free Energy Calculation and Per-Residue Decomposition via gmx_MMPBSA
 */

process GMX_MMPBSA {
    tag "${meta.id}_rep${meta.rep}"
    label 'process_high'

    input:
    tuple val(meta), path(clean_xtc), path(ref_pdb), path(topol_top), path(index_ndx)
    val(cores)
    val(start_ns)
    val(step_frames)
    val(temp)
    val(decomp)

    output:
    tuple val(meta), path("FINAL_RESULTS_MMPBSA.dat"), path("FINAL_RESULTS_MMPBSA.csv"), emit: results
    tuple val(meta), path("FINAL_DECOMP_MMPBSA.dat"), emit: decomp_results, optional: true

    script:
    def start_frame = ((start_ns * 1000) / 10) as int // assumes 10 ps recording interval
    def decomp_block = decomp ? """
&decomp
idecomp=1,
dec_verbose=1,
print_res="within 6.0"
/
""" : ""
    """
    cat << EOF > mmpbsa.in
&general
startframe=${start_frame},
interval=${step_frames},
temperature=${temp},
/
&gb
igb=5,
saltcon=0.150,
/
${decomp_block}
EOF

    # Execute gmx_MMPBSA using available cores
    mpirun -np ${cores} gmx_MMPBSA MPI \
        -O \
        -i mmpbsa.in \
        -cs ${ref_pdb} \
        -ci ${index_ndx} \
        -cg Protein UNK \
        -ct ${clean_xtc} \
        -cp ${topol_top} \
        -o FINAL_RESULTS_MMPBSA.dat \
        -eo FINAL_RESULTS_MMPBSA.csv \
        -do FINAL_DECOMP_MMPBSA.dat || true
    """

    stub:
    """
    touch FINAL_RESULTS_MMPBSA.dat
    touch FINAL_RESULTS_MMPBSA.csv
    touch FINAL_DECOMP_MMPBSA.dat
    """
}

