# FILE: src/utils/pymol_export.py
"""
File: src/utils/pymol_export.py
Description: Generates PyMOL sessions with 3.5 A pocket residues and Test #4 H-bonds.
"""
import re
import shutil
import pandas as pd
from pathlib import Path
from loguru import logger
from rdkit import Chem


def _find_3d_mol(cwd: Path, title: str, frame_id: str, is_native: bool):
    """Loads exact 3D RDKit Mol object keeping coordinates and hydrogens."""
    if is_native:
        candidates = [
            cwd / "01_prep_receptor" / "export" / f"ref_native_f{frame_id}.sdf",
            cwd / "04a_redocking" / "export" / f"glide_dock_f{frame_id}_poses.sdf",
            cwd / "04a_redocking" / "export" / f"gnina_dock_f{frame_id}_poses.sdf",
            cwd / "02_prep_ligands" / "export" / "native_prepared.sdf"
        ]
    else:
        dock_export = cwd / "04b_screening" / "export"
        candidates = [
            dock_export / f"glide_dock_f{frame_id}_poses.sdf",
            dock_export / f"gnina_dock_f{frame_id}_poses.sdf",
            dock_export / f"gnina_minimize_f{frame_id}_poses.sdf",
            dock_export / f"gnina_score_only_f{frame_id}_poses.sdf",
            cwd / "02_prep_ligands" / "export" / "library_prepared.sdf"
        ]

    clean_target = re.sub(r'(_model_\d+|_dock_\d+|_\d+)$', '', title).strip()

    for cand in candidates:
        if not cand.exists():
            continue
        try:
            suppl = Chem.SDMolSupplier(str(cand), removeHs=False, sanitize=False)
            for mol in suppl:
                if mol is None:
                    continue
                name = mol.GetProp("_Name") if mol.HasProp("_Name") else (
                    mol.GetProp("s_m_title") if mol.HasProp("s_m_title") else ""
                )
                clean_name = re.sub(r'(_model_\d+|_dock_\d+|_\d+)$', '', name).strip()
                if name == title or clean_name == clean_target or name == clean_target:
                    return mol
            if is_native and "ref_native" in cand.name:
                suppl = Chem.SDMolSupplier(str(cand), removeHs=False, sanitize=False)
                for mol in suppl:
                    if mol is not None:
                        return mol
        except Exception:
            continue
    return None


def generate_pymol_session(df: pd.DataFrame, cwd: Path, out_dir: Path):
    if df.empty:
        return

    pymol_dir = out_dir / "pymol_session"
    pymol_dir.mkdir(parents=True, exist_ok=True)
    pml_file = pymol_dir / "view_hits.pml"

    hit_items = []
    first_focus = None

    with open(pml_file, "w", encoding="utf-8") as pml:
        # 1. Global Viewport & Rendering Settings
        pml.write("# PyMOL SBDD Hit Session\n")
        pml.write("set bg_rgb, white\n")
        pml.write("set ray_opaque_background, on\n")
        pml.write("set defer_builds_mode, 3\n")
        
        # Fancy Helices
        pml.write("set cartoon_fancy_helices, 1\n")
        pml.write("set cartoon_highlight_color, grey90\n")
        pml.write("set cartoon_flat_sheets, 1\n")
        pml.write("set cartoon_smooth_loops, 1\n")
        pml.write("set cartoon_side_chain_helper, on\n")
        pml.write("set stick_radius, 0.22\n")
        pml.write("set valence, 1\n")
        pml.write("set auto_show_spheres, off\n")

        # Labels (Strictly in front)
        pml.write("set float_labels, on\n")
        pml.write("set label_position, [0.0, 0.0, 5.0]\n")
        pml.write("set label_color, black\n")
        pml.write("set label_outline_color, none\n")
        pml.write("set label_size, 18\n")
        pml.write("set label_font_id, 7\n")

        # H-Bond Dashes Styling
        pml.write("set dash_color, marine\n")
        pml.write("set dash_gap, 0.20\n")
        pml.write("set dash_width, 3.0\n")
        pml.write("set dash_radius, 0.06\n\n")

        # 2. Iterate hits
        for _, row in df.iterrows():
            raw_title = str(row["Title"])
            frame_id = str(row["Frame"]).replace("F", "")
            is_native = bool(row.get("Is_Native", False))
            is_bench = bool(row.get("Is_Benchmark", False))

            safe_code = re.sub(r'[^a-zA-Z0-9_]', '_', raw_title)

            if is_native:
                code_name = f"REF_Native_F{frame_id}"
            elif is_bench:
                code_name = f"CTRL_{safe_code}_F{frame_id}"
            else:
                code_name = f"{safe_code}_F{frame_id}"

            grp_name = code_name
            rec_obj = f"rec_{code_name}"
            lig_obj = f"lig_{code_name}"
            ref_obj = f"native_{code_name}"
            hb_obj = f"hb_{code_name}"

            if first_focus is None and not is_native:
                first_focus = lig_obj
            elif first_focus is None and is_native:
                first_focus = ref_obj

            # Copy receptor PDB
            src_rec = cwd / "01_prep_receptor" / "export" / f"receptor_prepared_f{frame_id}.pdb"
            dst_rec = pymol_dir / f"{rec_obj}.pdb"
            if src_rec.exists() and not dst_rec.exists():
                shutil.copy(src_rec, dst_rec)

            # Copy native reference SDF
            src_ref = cwd / "01_prep_receptor" / "export" / f"ref_native_f{frame_id}.sdf"
            dst_ref = pymol_dir / f"{ref_obj}.sdf"
            if src_ref.exists() and not dst_ref.exists():
                shutil.copy(src_ref, dst_ref)

            # Write docked ligand pose SDF
            dst_lig = pymol_dir / f"{lig_obj}.sdf"
            mol = _find_3d_mol(cwd, raw_title, frame_id, is_native)
            if mol is not None:
                writer = Chem.SDWriter(str(dst_lig))
                writer.write(mol)
                writer.close()

            # Load objects into PyMOL
            pml.write(f"# --- Load {code_name} ---\n")
            loaded_objs = []

            if dst_rec.exists():
                pml.write(f"load {dst_rec.name}, {rec_obj}\n")
                loaded_objs.append(rec_obj)

            if dst_ref.exists():
                pml.write(f"load {dst_ref.name}, {ref_obj}\n")
                loaded_objs.append(ref_obj)

            if dst_lig.exists():
                pml.write(f"load {dst_lig.name}, {lig_obj}\n")
                loaded_objs.append(lig_obj)

            active_lig = lig_obj if dst_lig.exists() else ref_obj

            if loaded_objs:
                hit_items.append((grp_name, rec_obj, active_lig, ref_obj, hb_obj, loaded_objs))

        # 3. Base Clean Styling (Global)
        pml.write("\n# --- Base Styling ---\n")
        pml.write("hide everything\n")
        pml.write("show cartoon, rec_*\n")
        pml.write("color gray85, rec_*\n")
        pml.write("color gray85, rec_* and elem c\n\n")

        pml.write("show sticks, lig_*\n")
        pml.write("color cyan, lig_* and elem c\n\n")

        pml.write("show sticks, native_*\n")
        pml.write("color forest, native_* and elem c\n\n")

        # 4. Pocket 3.5 A Residues + Method 4 H-Bonds
        pml.write("# --- Pocket Contacts (3.5 A) & Method 4 H-Bonds ---\n")
        for grp_name, rec_obj, active_lig, ref_obj, hb_obj, loaded_objs in hit_items:
            sel_pocket = f"p_{grp_name}"
            
            # Simple direct 3.5 A shell selection
            pml.write(f"select {sel_pocket}, (byres ({rec_obj} within 3.5 of {active_lig}))\n")
            pml.write(f"show sticks, ({sel_pocket} and not name c+n+o)\n")
            pml.write(f"set stick_radius, 0.16, {sel_pocket}\n")
            pml.write(f"color gray70, {sel_pocket} and elem c\n")
            pml.write(f'label ({sel_pocket} and (name CB or (name CA and resn GLY))), "%s%s" % (resn, resi)\n')

            # Method 4 H-Bonds
            pml.write(f"dist {hb_obj}, ({active_lig} and elem n+o+s+f+cl+br), ({rec_obj} and elem n+o+s), 3.5\n")
            pml.write(f"show dashes, {hb_obj}\n")
            pml.write(f"hide labels, {hb_obj}\n")

            # Group everything into compound name
            group_members = [o for o in loaded_objs] + [hb_obj]
            pml.write(f"group {grp_name}, {' '.join(group_members)}\n\n")

        pml.write("delete p_*\n")
        pml.write("deselect\n\n")

        # 5. Viewport Focus
        pml.write("# --- Initial Camera View ---\n")
        pml.write("disable all\n")
        if hit_items:
            pml.write(f"enable {hit_items[0][0]}\n")
            target = first_focus if first_focus else hit_items[0][0]
            pml.write(f"zoom {target}, 5\n")

    logger.success(f"PyMOL Session cleanly exported to: {pml_file.absolute()}")