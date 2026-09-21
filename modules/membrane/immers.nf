/*
 * modules/membrane/immers.nf
 * PPM 3.0 Standalone Membrane Immersion Calculation
 * Pure external scientific tool (No Schroedinger dependency)
 */

process IMMERS {
    tag "${meta.id}"
    label "process_medium"
    label "process_master"

    input:
    tuple val(meta), path(protein_pdb)
    path(parse_script)
    val(membrane_type)
    val(membrane_topo)
    val(membrane_chains)

    output:
    tuple val(meta), path("${meta.id}_oriented.pdb"), emit: oriented_pdb
    tuple val(meta), path("${meta.id}_membrane.json"), emit: membrane_json
    tuple val(meta), path("${meta.id}_membrane.txt"), emit: membrane_txt

    script:
    """
    if [ ! -x "/usr/local/bin/immers" ] || [ ! -f "/usr/local/share/immers/res.lib" ]; then
        echo "ERROR: immers binary or parameter library missing on host" >&2
        exit 1
    fi

    ln -sf /usr/local/share/immers/res.lib res.lib

    cat <<EOF > immers.inp
2
no
${protein_pdb}
1
${membrane_type}
planar
${membrane_topo}
${membrane_chains}
EOF

    /usr/local/bin/immers < immers.inp > immers.log

    pdb_base="\$(basename ${protein_pdb} .pdb)"
    if [ -f "\${pdb_base}out.pdb" ]; then
        mv "\${pdb_base}out.pdb" "${meta.id}_oriented.pdb"
    elif [ -f "protout.pdb" ]; then
        mv "protout.pdb" "${meta.id}_oriented.pdb"
    else
        echo "ERROR: immers did not produce expected oriented PDB." >&2
        cat immers.log >&2
        exit 1
    fi

    python3 ${parse_script} \
        --datapar datapar1 \
        --out_json ${meta.id}_membrane.json \
        --out_txt ${meta.id}_membrane.txt
    """

    stub:
    """
    touch ${meta.id}_oriented.pdb
    echo '{"thickness": 28.0, "tilt": 0.0, "dg": 0.0}' > ${meta.id}_membrane.json
    echo '28.0' > ${meta.id}_membrane.txt
    """
}
