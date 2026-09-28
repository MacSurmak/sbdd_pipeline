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
    // For omnivorous setup, we dynamically extract all non-protein organics
    """
    echo "q" | gmx make_ndx -f ${input_gro} -o index.ndx
    python3 ${projectDir}/bin/gromacs_ndx_fixer.py --ndx index.ndx --action add_groups

    # 1. Protein Restraints (Isolated to guarantee 1-to-N local index)
    echo "Protein" | gmx trjconv -f ${input_gro} -n index.ndx -o prot.gro
    echo "System" | gmx genrestr -f prot.gro -o posre_PROT.itp -fc 9999 9999 9999
    sed -i 's/9999/POSRES_PROT_FC/g' posre_PROT.itp

    # 2. Organic Restraints
    ORG_LIST=""
    grep -E "^\\[" index.ndx | tr -d '[] ' > groups.txt
    for grp in \$(cat groups.txt); do
        if [[ "\$grp" != "System" && "\$grp" != "Protein" && "\$grp" != "SOL" && "\$grp" != "Water" && "\$grp" != "PROT_LIG" && "\$grp" != "SOL_ION" && "\$grp" != "COMPLEX" ]]; then
            # Ignore single atom ions
            NATOMS=\$(echo "\$grp" | gmx make_ndx -f ${input_gro} -n index.ndx -o /dev/null 2>&1 | grep -E "^\s*[0-9]+ \$grp" | awk '{print \$3}')
            if [ -n "\$NATOMS" ] && [ "\$NATOMS" -gt 1 ]; then
                echo "\$grp" | gmx trjconv -f ${input_gro} -n index.ndx -o \${grp}.gro
                echo "System" | gmx genrestr -f \${grp}.gro -o posre_\${grp}.itp -fc 9999 9999 9999
                sed -i "s/9999/POSRES_\${grp}_FC/g" posre_\${grp}.itp
                ORG_LIST="\$ORG_LIST,\$grp"
            fi
        fi
    done

    python3 ${projectDir}/bin/gromacs_ndx_fixer.py --top ${topol_top} --action inject_posre --orgs "\$ORG_LIST"
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
