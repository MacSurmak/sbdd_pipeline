/*
 * modules/gromacs/ndx_posre.nf
 * Generates GROMACS index.ndx, thermostat groups PROT_LIG/SOL_ION, and multi-level position restraints
 */

process NDX_POSRE {
    tag "${meta.id}"
    label 'process_low'

    input:
    tuple val(meta), path(input_gro), path(topol_top)

    output:
    tuple val(meta), path("index.ndx"), path("topol.top"), path("posre*.itp"), emit: restraints

    script:
    """
    echo "q" | gmx make_ndx -f ${input_gro} -o index.ndx

    python3 - << 'EOF'
    import sys
    from pathlib import Path
    sys.path.insert(0, "${projectDir}/bin")
    import gromacs_ndx_fixer as g
    g.add_thermostat_groups_by_name(Path("index.ndx"))
    EOF

    # Generate multi-level restraints for Protein and Ligand
    for fc in "1000:" "250:_250" "50:_50"; do
        val=\${fc%%:*}
        sfx=\${fc##*:}
        echo "Protein" | gmx genrestr -f ${input_gro} -n index.ndx -o "posre\${sfx}.itp" -fc \$val \$val \$val
        echo "UNK"     | gmx genrestr -f ${input_gro} -n index.ndx -o "posre_lig\${sfx}.itp" -fc \$val \$val \$val
    done

    python3 - << 'EOF'
    import sys
    from pathlib import Path
    sys.path.insert(0, "${projectDir}/bin")
    import gromacs_ndx_fixer as g
    for sfx in ["", "_250", "_50"]:
        g.fix_posre_indices(Path(f"posre{sfx}.itp"))
        g.fix_posre_indices(Path(f"posre_lig{sfx}.itp"))
    g.inject_ligand_posre(Path("${topol_top}"))
    EOF
    """

    stub:
    """
    touch index.ndx
    touch topol.top
    touch posre.itp
    touch posre_250.itp
    touch posre_50.itp
    """
}

