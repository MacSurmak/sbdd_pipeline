/*
 * modules/membrane/immers.nf
 * PPM 3.0 Standalone Membrane Immersion Calculation
 * Self-contained worker process: zero external script file staging.
 */

process IMMERS {
    tag "${meta.id}"
    label "process_medium"

    input:
    tuple val(meta), path(protein_pdb)
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

    cat << 'EOF' > parse_immers.py
import argparse, json, sys
from pathlib import Path

parser = argparse.ArgumentParser()
parser.add_argument("--datapar", default="datapar1")
parser.add_argument("--out_json", required=True)
parser.add_argument("--out_txt", required=True)
args = parser.parse_args()

data = {"thickness": 28.0, "tilt": 0.0, "dg": 0.0}
datapar_path = Path(args.datapar)
if datapar_path.exists():
    try:
        with open(datapar_path, "r") as f:
            parts = [p.strip() for p in f.readline().split(";")]
            data["thickness"] = float(parts[1])
            data["tilt"] = float(parts[3])
            data["dg"] = float(parts[5])
    except Exception:
        pass

with open(args.out_json, "w") as f_json:
    json.dump(data, f_json, indent=2)

with open(args.out_txt, "w") as f_txt:
    f_txt.write(str(data["thickness"]))
EOF

    python3 parse_immers.py \
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
