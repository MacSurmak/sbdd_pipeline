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
    python3 ${projectDir}/bin/gromacs_ndx_fixer.py --ndx index.ndx --action add_groups

    # 1. Protein Restraints (Isolated to guarantee 1-to-N local index)
    if grep -q "\\[ Protein \\]" index.ndx; then
        echo "Protein" | gmx trjconv -s ${input_gro} -f ${input_gro} -n index.ndx -o prot.gro
        echo "0" | gmx genrestr -f prot.gro -o posre_PROT.itp
        sed -i -E 's/1000\\s+1000\\s+1000/POSRES_PROT_FC POSRES_PROT_FC POSRES_PROT_FC/g' posre_PROT.itp
    fi

    # 2. Organic Restraints from Topology
    ORG_LIST=""
    awk '/^\\[ molecules \\]/{flag=1; next} flag && NF>0 && !/^;/{print \$1}' ${topol_top} | sort | uniq > mol_types.txt

    for grp in \$(cat mol_types.txt); do
        if [[ "\$grp" == "system1" || "\$grp" == "Protein" || "\$grp" == "SOL" || "\$grp" == "Water" ]]; then continue; fi
        if echo "\$grp" | grep -qE "^(K\\+|NA\\+|CL-|Cl-|Na\\+|ZN|MG|CA|ZN2|MG2|CA2)\$"; then continue; fi

        if grep -q "^\\[ \$grp \\]" index.ndx; then
            echo "\$grp" | gmx trjconv -s ${input_gro} -f ${input_gro} -n index.ndx -o \${grp}.gro
            echo "0" | gmx genrestr -f \${grp}.gro -o posre_\${grp}.itp
            sed -i -E "s/1000\\s+1000\\s+1000/POSRES_\${grp}_FC POSRES_\${grp}_FC POSRES_\${grp}_FC/g" posre_\${grp}.itp
            ORG_LIST="\$ORG_LIST,\$grp"
        fi
    done

    python3 ${projectDir}/bin/gromacs_ndx_fixer.py --ndx index.ndx --top ${topol_top} --action inject_posre --orgs "\$ORG_LIST"
    """

    stub:
    """
    touch index.ndx
    touch topol.top
    touch posre_PROT.itp
    """
}