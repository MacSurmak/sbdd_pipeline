/*
 * modules/membrane/apply_membrane.nf
 * Implicit Membrane Property Injection for Prime
 * Self-contained worker process: zero external script file staging.
 */

process APPLY_MEMBRANE {
    tag "${meta.id}"
    label "process_single"

    input:
    tuple val(meta), path(prepared_mae), path(membrane_info)

    output:
    tuple val(meta), path("${meta.id}_membrane.mae"), emit: membrane_mae

    script:
    """
    cat << 'EOF' > apply_membrane.py
import argparse, sys, json
from pathlib import Path
from schrodinger import structure

parser = argparse.ArgumentParser()
parser.add_argument("--input", required=True)
parser.add_argument("--info", required=True)
parser.add_argument("--output", required=True)
args = parser.parse_args()

st = structure.StructureReader.read(args.input)
thickness = 28.0

info_path = Path(args.info)
if info_path.exists():
    content = info_path.read_text().strip()
    if content.startswith("{"):
        thickness = float(json.loads(content).get("thickness", 28.0))
    else:
        thickness = float(content)

half_t = thickness / 2.0
st.property['r_psp_Memb1_x'] = 0.0
st.property['r_psp_Memb1_y'] = 0.0
st.property['r_psp_Memb1_z'] = -half_t

st.property['r_psp_Memb2_x'] = 0.0
st.property['r_psp_Memb2_y'] = 0.0
st.property['r_psp_Memb2_z'] = half_t

st.property['r_psp_Prime_Membrane_Thickness'] = thickness
st.property['r_psp_Prime_Membrane_Buffer'] = 2.5

st.write(args.output)
EOF

    schrodinger run python3 apply_membrane.py \\
        --input ${prepared_mae} \\
        --info ${membrane_info} \\
        --output ${meta.id}_membrane.mae
    """

    stub:
    """
    touch ${meta.id}_membrane.mae
    """
}
