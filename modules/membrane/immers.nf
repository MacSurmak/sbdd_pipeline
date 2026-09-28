/*
 * modules/membrane/immers.nf
 * PPM 3.0 Standalone Membrane Immersion Calculation
 * Self-contained worker process: zero external script file staging.
 */

process IMMERS {
    tag "${meta.id}"
    label "process_medium"

    input:
    tuple val(meta), path(pure_protein_pdb)
    val(membrane_type)
    val(membrane_topo)
    val(membrane_chains)

    output:
    tuple val(meta), path("${meta.id}_oriented.pdb"), emit: oriented_pdb
    tuple val(meta), path("${meta.id}_membrane.json"), emit: membrane_json
    tuple val(meta), path("${meta.id}_membrane.txt"), emit: membrane_txt

    script:
    """
    ln -sf /usr/local/share/immers/res.lib res.lib 2>/dev/null || true

    cat <<EOF > immers.inp
2
no
${pure_protein_pdb}
1
${membrane_type}
planar
${membrane_topo}
${membrane_chains}
EOF

    IMMERS_BIN=\$(which immers 2>/dev/null || echo "/opt/containers/bin/immers")
    if [ ! -x "\$IMMERS_BIN" ]; then
        IMMERS_BIN="/usr/local/bin/immers"
    fi

    if [ ! -x "\$IMMERS_BIN" ]; then
        echo "FATAL: immers binary not found in PATH, /opt/containers/bin, or /usr/local/bin!" >&2
        exit 1
    fi

    "\$IMMERS_BIN" < immers.inp > immers.log

    pdb_base="\$(basename ${pure_protein_pdb} .pdb)"
    if [ -f "\${pdb_base}out.pdb" ]; then
        mv "\${pdb_base}out.pdb" "${meta.id}_oriented.pdb"
    elif [ -f "protout.pdb" ]; then
        mv "protout.pdb" "${meta.id}_oriented.pdb"
    else
        echo "FATAL: immers did not produce expected oriented PDB." >&2
        cat immers.log >&2
        exit 1
    fi

    python3 ${projectDir}/bin/parse_immers.py \\
        --datapar datapar1 \\
        --out_json ${meta.id}_membrane.json \\
        --out_txt ${meta.id}_membrane.txt
    """

    stub:
    """
    touch ${meta.id}_oriented.pdb
    echo '{"thickness": 28.0, "tilt": 0.0, "dg": 0.0}' > ${meta.id}_membrane.json
    echo '28.0' > ${meta.id}_membrane.txt
    """
}
