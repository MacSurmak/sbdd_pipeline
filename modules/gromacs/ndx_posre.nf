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
    cat << 'EOF' > gromacs_ndx_fixer.py
import os, sys
from pathlib import Path
from typing import Dict, List, Set

def parse_ndx(ndx_path: Path) -> Dict[str, List[int]]:
    groups: Dict[str, List[int]] = {}
    current_group = None
    with open(ndx_path, "r") as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
            if line.startswith("[") and line.endswith("]"):
                current_group = line[1:-1].strip()
                groups[current_group] = []
            elif current_group is not None:
                groups[current_group].extend([int(x) for x in line.split()])
    return groups

def write_ndx_group(f, name: str, atom_ids: List[int]):
    f.write(f"\\n[ {name} ]\\n")
    for i in range(0, len(atom_ids), 15):
        f.write(" ".join(f"{x:5d}" for x in atom_ids[i:i+15]) + "\\n")

def add_thermostat_groups_by_name(ndx_path: Path):
    groups = parse_ndx(ndx_path)
    if "System" not in groups:
        raise KeyError("Group [ System ] not found in index.ndx")

    water_and_ion_atoms: Set[int] = set()
    ion_names = {"K+", "Cl-", "NA+", "CL-", "NA", "CL", "K", "POT", "CLA", "SOD", "MG", "ZN", "CA", "CAL"}
    if "SOL" in groups:
        water_and_ion_atoms.update(groups["SOL"])
    elif "Water" in groups:
        water_and_ion_atoms.update(groups["Water"])
    else:
        raise KeyError("Neither [ SOL ] nor [ Water ] found in index.ndx")

    for ion_name in ion_names:
        if ion_name in groups:
            water_and_ion_atoms.update(groups[ion_name])

    system_atoms = set(groups["System"])
    solute_membrane_atoms = system_atoms - water_and_ion_atoms

    if not solute_membrane_atoms:
        raise ValueError("PROT_LIG group is empty! Check index.ndx names.")
    if not water_and_ion_atoms:
        raise ValueError("SOL_ION group is empty! Check index.ndx names.")
    if len(solute_membrane_atoms) + len(water_and_ion_atoms) != len(system_atoms):
        raise ValueError("Partition mismatch: PROT_LIG + SOL_ION != System!")

    with open(ndx_path, "a") as f:
        write_ndx_group(f, "PROT_LIG", sorted(list(solute_membrane_atoms)))
        write_ndx_group(f, "SOL_ION", sorted(list(water_and_ion_atoms)))

def inject_ligand_posre(top_file: Path):
    if not top_file.exists(): return
    lines = top_file.read_text().splitlines()
    out_lines = []
    in_unk = False
    injected = False
    posre_block = (
        "\\n#ifdef POSRES_LIG\\n"
        "#include \\"posre_lig.itp\\"\\n"
        "#endif\\n"
        "#ifdef POSRES_LIG_250\\n"
        "#include \\"posre_lig_250.itp\\"\\n"
        "#endif\\n"
        "#ifdef POSRES_LIG_50\\n"
        "#include \\"posre_lig_50.itp\\"\\n"
        "#endif\\n"
    )
    for line in lines:
        stripped = line.strip()
        if in_unk and not injected and (stripped.startswith("[ moleculetype ]") or stripped.startswith("[ system ]")):
            out_lines.append(posre_block)
            injected = True
            in_unk = False
        if stripped.startswith("[ moleculetype ]"):
            in_unk = False
        if not stripped.startswith(";") and stripped.startswith("UNK") and len(stripped.split()) >= 2:
            in_unk = True
        out_lines.append(line)
    if in_unk and not injected:
        out_lines.append(posre_block)

    out_lines_prot = []
    in_prot = False
    injected_prot = False
    posre_block_prot = (
        "\\n#ifdef POSRES\\n"
        "#include \\"posre.itp\\"\\n"
        "#endif\\n"
        "#ifdef POSRES_250\\n"
        "#include \\"posre_250.itp\\"\\n"
        "#endif\\n"
        "#ifdef POSRES_50\\n"
        "#include \\"posre_50.itp\\"\\n"
        "#endif\\n"
    )
    for line in out_lines:
        stripped = line.strip()
        if in_prot and not injected_prot and (stripped.startswith("[ moleculetype ]") or stripped.startswith("[ system ]")):
            out_lines_prot.append(posre_block_prot)
            injected_prot = True
            in_prot = False
        if stripped.startswith("[ moleculetype ]"):
            in_prot = False
        if not stripped.startswith(";") and (stripped.startswith("system1") or stripped.startswith("Protein")) and len(stripped.split()) >= 2:
            in_prot = True
        out_lines_prot.append(line)
    top_file.write_text("\\n".join(out_lines_prot) + "\\n")

def fix_posre_indices(itp_file: Path):
    if not itp_file.exists(): return
    lines = itp_file.read_text().splitlines()
    out_lines = []
    atom_idx = 1
    for line in lines:
        stripped = line.strip()
        if not stripped or stripped.startswith(";") or stripped.startswith("["):
            out_lines.append(line)
            continue
        parts = stripped.split()
        if len(parts) >= 5:
            out_lines.append(f"{atom_idx:5d}{int(parts[1]):5d}{float(parts[2]):10.1f}{float(parts[3]):10.1f}{float(parts[4]):10.1f}")
            atom_idx += 1
        else:
            out_lines.append(line)
    itp_file.write_text("\\n".join(out_lines) + "\\n")
EOF

    echo "q" | gmx make_ndx -f ${input_gro} -o index.ndx

    python3 - << 'EOF'
    import sys
    from pathlib import Path
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
    from pathlib import Path
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

