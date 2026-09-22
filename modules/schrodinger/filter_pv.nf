/*
 * modules/schrodinger/filter_pv.nf
 * Filters Pose Viewer file based on funnel whitelist CSV for a specific frame,
 * splitting into chunked PV files for massive distributed Prime MM-GBSA execution.
 * Self-contained worker process: zero external script file staging.
 */

process FILTER_PV {
    tag "${meta.id}"
    label 'process_single'

    input:
    tuple val(meta), path(pv_maegz)
    path(whitelist_csv)
    val(batch_size)

    output:
    tuple val(meta), path("${meta.id}_filtered_pv.maegz"), emit: filtered_pv
    tuple val(meta), path("${meta.id}_filtered_c*.maegz"), emit: filtered_pv_chunks

    script:
    def frame_id = meta.frame ?: 1
    """
    cat << 'EOF' > filter_pv.py
import argparse, sys, csv
from schrodinger import structure

parser = argparse.ArgumentParser()
parser.add_argument("--input", required=True)
parser.add_argument("--whitelist", required=True)
parser.add_argument("--frame", required=True)
parser.add_argument("--max_poses", type=int, default=1)
parser.add_argument("--out_prefix", default="")
parser.add_argument("--batch_size", type=int, default=25)
parser.add_argument("--output", default="")
args = parser.parse_args()

allowed_titles = set()
with open(args.whitelist, "r") as f:
    for row in csv.DictReader(f):
        if str(row["Frame"]) == str(args.frame):
            allowed_titles.add(row["Title"])

receptor = None
whitelisted_poses = []
pose_counts = {}

with structure.StructureReader(args.input) as reader:
    for i, st in enumerate(reader):
        if i == 0:
            receptor = st.copy()
            continue
        title = (st.property.get("s_m_title") or st.title or "").strip()
        if title in allowed_titles and pose_counts.get(title, 0) < args.max_poses:
            whitelisted_poses.append(st.copy())
            pose_counts[title] = pose_counts.get(title, 0) + 1

if args.output:
    with structure.StructureWriter(args.output) as writer:
        if receptor:
            writer.append(receptor)
        for p in whitelisted_poses:
            writer.append(p)

if args.out_prefix:
    if not whitelisted_poses:
        c0_file = f"{args.out_prefix}_c0.maegz"
        with structure.StructureWriter(c0_file) as writer:
            if receptor:
                writer.append(receptor)
    else:
        chunk_idx = 0
        for i in range(0, len(whitelisted_poses), args.batch_size):
            chunk = whitelisted_poses[i : i + args.batch_size]
            chunk_file = f"{args.out_prefix}_c{chunk_idx}.maegz"
            with structure.StructureWriter(chunk_file) as writer:
                if receptor:
                    writer.append(receptor)
                for p in chunk:
                    writer.append(p)
            chunk_idx += 1
EOF

    schrodinger run python3 filter_pv.py \\
        --input ${pv_maegz} \\
        --output ${meta.id}_filtered_pv.maegz \\
        --out_prefix ${meta.id}_filtered \\
        --whitelist ${whitelist_csv} \\
        --frame ${frame_id} \\
        --batch_size ${batch_size} \\
        --max_poses 1
    """

    stub:
    """
    touch ${meta.id}_filtered_pv.maegz
    touch ${meta.id}_filtered_c0.maegz
    """
}
