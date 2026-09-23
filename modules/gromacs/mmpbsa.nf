/*
 * modules/gromacs/mmpbsa.nf
 * Free Energy Calculation (Membrane PB) and Per-Residue Decomposition (GB) via gmx_MMPBSA
 * Accurately reproducing original physical parameters and thermodynamic models.
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
    val(salt_conc)
    val(is_membrane)
    val(membrane_thickness)
    val(decomp)

    output:
    tuple val(meta), path("thermodynamics_mmpbsa_pb.csv"), path("thermodynamics_mmpbsa_pb.dat"), emit: pb_results
    tuple val(meta), path("decomposition_mmpbsa_gb.csv"), path("thermodynamics_mmpbsa_gb.dat"), emit: gb_decomp_results, optional: true

    script:
    def mthick_val = membrane_thickness ? (membrane_thickness as float) : 28.0
    def pb_block = is_membrane ? """
&pb
memopt=1,
mthick=${mthick_val},
mctrdz=0.0,
emem=7.0,
indi=4.0,
istrng=${salt_conc},
radiopt=0,
poretype=1,
fillratio=1.25,
ipb=1,
inp=2,
sasopt=0,
solvopt=2,
bcopt=10,
nfocus=1,
linit=1000,
eneopt=1,
cutfd=7.0,
cutnb=99.0,
maxarcdot=15000
/
""" : """
&gb
igb=2,
saltcon=${salt_conc}
/
"""
    def mpirun_cmd = cores > 1 ? "mpirun --oversubscribe -np ${cores} gmx_MMPBSA MPI" : "gmx_MMPBSA"
    """
    # 1. Total Binding Free Energy via Membrane-aware Poisson-Boltzmann (PB)
    start_frame=\$(get_mda_startframe.py --pdb ${ref_pdb} --xtc ${clean_xtc} --start_ns ${start_ns})

    cat << EOF > mmpbsa_pb.in
&general
sys_name="Protein_Ligand",
startframe=\${start_frame},
interval=${step_frames},
temperature=${temp}
/
${pb_block}
EOF

    ${mpirun_cmd} \\
        -O \\
        -i mmpbsa_pb.in \\
        -cs ${ref_pdb} \\
        -ci ${index_ndx} \\
        -cg Protein UNK \\
        -ct ${clean_xtc} \\
        -cp ${topol_top} \\
        -nogui \\
        -o thermodynamics_mmpbsa_pb.dat \\
        -eo thermodynamics_mmpbsa_pb.csv

    # 2. Per-Residue Decomposition via Generalized Born (GB)
    if [ "${decomp}" = "true" ]; then
        cat << EOF > mmpbsa_gb.in
&general
sys_name="Protein_Ligand",
startframe=\${start_frame},
interval=${step_frames},
temperature=${temp}
/
&gb
igb=2,
saltcon=${salt_conc}
/
&decomp
idecomp=2,
dec_verbose=0
/
EOF

        ${mpirun_cmd} \\
            -O \\
            -prefix _GMXMMPBSA_GB_ \\
            -i mmpbsa_gb.in \\
            -cs ${ref_pdb} \\
            -ci ${index_ndx} \\
            -cg Protein UNK \\
            -ct ${clean_xtc} \\
            -cp ${topol_top} \\
            -nogui \\
            -o thermodynamics_mmpbsa_gb.dat \\
            -deo decomposition_mmpbsa_gb.csv
    fi
    """

    stub:
    """
    touch thermodynamics_mmpbsa_pb.csv
    touch thermodynamics_mmpbsa_pb.dat
    touch decomposition_mmpbsa_gb.csv
    touch thermodynamics_mmpbsa_gb.dat
    """
}

