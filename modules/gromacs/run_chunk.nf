/*
 * modules/gromacs/run_chunk.nf
 * Executes production MD chunks on Slurm GPU
 */

process RUN_MD_INIT {
    tag "${meta.id}_rep${meta.rep}"
    label 'process_gpu'

    input:
    tuple val(meta), path(in_gro), path(in_cpt), path(topol_top), path(index_ndx)
    val(is_membrane)
    val(chunk_ns)
    val(temperature)

    output:
    tuple val(meta), path("md_1.gro"), path("md_1.cpt"), path("md_1.tpr"), emit: next_chunk_input
    tuple val(meta), path("md_1.xtc"), emit: xtc
    tuple val(meta), path("md_1.edr"), emit: edr
    tuple val(meta), path("md_1.log"), emit: log

    script:
    def nsteps = ((chunk_ns * 1000) / 0.002) as int
    def pcoupl_type = is_membrane ? "semiisotropic" : "isotropic"
    def ref_p = is_membrane ? "1.0 1.0" : "1.0"
    def comp = is_membrane ? "4.5e-5 4.5e-5" : "4.5e-5"
    """
    echo "=== Worker: \$(hostname) - Executing Chunk 1 for ${meta.id} (Rep ${meta.rep}) ==="
    nvidia-smi --query-gpu=name,memory.total --format=csv,noheader || true

    cat << EOF > md.mdp
; Production MD template
integrator               = md
dt                       = 0.002
nsteps                   = ${nsteps}

nstxout-compressed       = 5000
nstcalcenergy            = 100
nstenergy                = 1000
nstlog                   = 1000

cutoff-scheme            = Verlet
nstlist                  = 20
coulombtype              = PME
rcoulomb                 = 1.0
vdwtype                  = Cut-off
vdw-modifier             = None
DispCorr                 = EnerPres
rvdw                     = 1.0
rlist                    = 1.0
pbc                      = xyz

tcoupl                   = v-rescale
tc-grps                  = PROT_LIG SOL_ION
tau-t                    = 1.0 1.0
ref-t                    = ${temperature} ${temperature}

pcoupl                   = Parrinello-Rahman
pcoupltype               = ${pcoupl_type}
tau-p                    = 5.0
ref-p                    = ${ref_p}
compressibility          = ${comp}

constraints              = h-bonds
constraint-algorithm     = LINCS
continuation             = yes
lincs-order              = 4
lincs-iter               = 1

nstcomm                  = 100
comm-mode                = linear
comm-grps                = PROT_LIG SOL_ION

gen-vel                  = yes
gen-temp                 = ${temperature}
gen-seed                 = ${meta.rep * 12345 + 101}
EOF

    gmx grompp -f md.mdp -o md_1.tpr -c ${in_gro} -t ${in_cpt} -p ${topol_top} -n ${index_ndx} -maxwarn 2
    gmx mdrun -deffnm md_1 -nb gpu -pme gpu -update gpu -ntomp 8 -noappend
    """

    stub:
    """
    touch md_1.gro
    touch md_1.cpt
    touch md_1.tpr
    touch md_1.xtc
    touch md_1.edr
    touch md_1.log
    """
}

process RUN_MD_CHUNK {
    tag "${meta.id}_rep${meta.rep}_chunk${chunk_idx}"
    label 'process_gpu'

    input:
    tuple val(meta), path(in_cpt), path(in_tpr), val(chunk_idx)
    val(chunk_ns)

    output:
    tuple val(meta), path("md_${chunk_idx}.gro"), path("md_${chunk_idx}.cpt"), path("md_${chunk_idx}.tpr"), emit: next_chunk_input
    tuple val(meta), path("md_${chunk_idx}.xtc"), emit: xtc
    tuple val(meta), path("md_${chunk_idx}.edr"), emit: edr
    tuple val(meta), path("md_${chunk_idx}.log"), emit: log

    script:
    def extend_ps = (chunk_ns * 1000) as int
    """
    echo "=== Worker: \$(hostname) - Extending Chunk ${chunk_idx} for ${meta.id} (Rep ${meta.rep}) ==="
    nvidia-smi --query-gpu=name,memory.total --format=csv,noheader || true

    gmx convert-tpr -s ${in_tpr} -extend ${extend_ps} -o md_${chunk_idx}.tpr
    gmx mdrun -deffnm md_${chunk_idx} -cpi ${in_cpt} -nb gpu -pme gpu -update gpu -ntomp 8 -noappend
    """

    stub:
    """
    touch md_${chunk_idx}.gro
    touch md_${chunk_idx}.cpt
    touch md_${chunk_idx}.tpr
    touch md_${chunk_idx}.xtc
    touch md_${chunk_idx}.edr
    touch md_${chunk_idx}.log
    """
}
