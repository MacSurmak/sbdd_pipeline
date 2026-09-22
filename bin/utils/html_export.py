"""
File: src/utils/html_export.py
Description: Generates standalone HTML report with 2D structures using RDKit and templates.
"""
import os
import pandas as pd
import numpy as np
from pathlib import Path
from loguru import logger
from rdkit import Chem, RDLogger
from rdkit.Chem import rdDepictor
from rdkit.Chem.Draw import rdMolDraw2D

# Suppress internal C++ parser warnings from RDKit
RDLogger.DisableLog('rdApp.*')

def _mol_to_svg(mol, width=180, height=130) -> str:
    if mol is None: return '<div style="color:#aaa;font-size:11px;text-align:center;">No 2D Structure</div>'
    try:
        m = Chem.Mol(mol)
        rdDepictor.Compute2DCoords(m)
        drawer = rdMolDraw2D.MolDraw2DSVG(width, height)
        opts = drawer.drawOptions()
        opts.clearBackground = True
        opts.bondLineWidth = 2
        drawer.DrawMolecule(m)
        drawer.FinishDrawing()
        svg = drawer.GetDrawingText()
        if "<?xml" in svg: svg = svg[svg.find("<svg"):]
        return svg
    except Exception:
        return '<div style="color:#e74c3c;font-size:11px;text-align:center;">Render Error</div>'

def _load_structures_for_hits(df_hits: pd.DataFrame, cwd: Path) -> dict:
    mols = {}
    needed = df_hits[["Title", "Frame", "Is_Native"]].to_dict("records")
    grouped = {}
    for item in needed:
        key = (item["Is_Native"], item["Frame"].replace("F", ""))
        grouped.setdefault(key, set()).add(item["Title"])
        
    for (is_native, frame_id), titles in grouped.items():
        if is_native:
            candidates = [
                cwd / "01_prep_receptor" / "export" / f"ref_native_f{frame_id}.sdf",
                cwd / "02_prep_ligands" / "export" / "native_prepared.sdf"
            ]
        else:
            dock_export = cwd / "04b_screening" / "export"
            candidates = [
                dock_export / f"glide_dock_f{frame_id}_poses.sdf",
                dock_export / f"gnina_dock_f{frame_id}_poses.sdf",
                dock_export / f"gnina_score_only_f{frame_id}_poses.sdf",
                dock_export / f"gnina_minimize_f{frame_id}_poses.sdf",
                cwd / "02_prep_ligands" / "export" / "library_prepared.sdf"
            ]
            
        found_titles = set()
        for sdf_path in candidates:
            if not sdf_path.exists(): continue
            suppl = Chem.SDMolSupplier(str(sdf_path), removeHs=True)
            for mol in suppl:
                if mol is None: continue
                name = mol.GetProp("_Name") if mol.HasProp("_Name") else mol.GetProp("s_m_title") if mol.HasProp("s_m_title") else ""
                if name in titles and name not in mols:
                    mols[name] = mol
                    found_titles.add(name)
            if titles.issubset(found_titles): break
    return mols

def generate_html_report(df: pd.DataFrame, cwd: Path, out_dir: Path):
    if df.empty: return
    
    mols_map = _load_structures_for_hits(df, cwd)
    fdir_env = os.environ.get("SBDD_FRAMEWORK_DIR")
    framework_dir = Path(fdir_env) if fdir_env else Path(__file__).resolve().parents[2]
    template_path = framework_dir / "templates" / "report_template.html"
    
    if not template_path.exists():
        logger.error(f"HTML Template missing at {template_path}. Report skipped.")
        return
        
    rows_html = []
    for rank, (_, row) in enumerate(df.iterrows(), start=1):
        raw_title = row["Title"]
        is_native = row.get("Is_Native", False)
        title = "Native Reference" if is_native else raw_title
        svg_code = _mol_to_svg(mols_map.get(raw_title))
        
        is_bench = row.get("Is_Benchmark", False)
        row_class = "native-row" if is_native else "control-row" if is_bench else ""
        badge = '<span class="badge-native">REF</span>' if is_native else '<span class="badge-control">CTRL</span>' if is_bench else f'<span class="badge-rank">#{rank}</span>'
        
        def fmt(val, dec=2): return f"{val:.{dec}f}" if pd.notna(val) else '<span class="na">—</span>'
        
        # Format Delta to Native
        delta_val = row.get("Δ_to_Native", np.nan)
        delta_class = "delta-pos" if pd.notna(delta_val) and delta_val > 0 else "delta-neutral"
        delta_str = f"+{delta_val:.2f}" if pd.notna(delta_val) and delta_val > 0 else fmt(delta_val)
        
        # Format Pose Consensus RMSD (Glide vs GNINA dock)
        p_rmsd = row.get("Pose_RMSD", np.nan)
        if pd.isna(p_rmsd):
            rmsd_badge = '<span class="na">—</span>'
        elif p_rmsd <= 2.0:
            rmsd_badge = f'<span class="col-num" style="color:#16a34a; font-weight:600;">{p_rmsd:.2f} Å</span>'
        else:
            rmsd_badge = f'<span class="col-num">{p_rmsd:.2f} Å</span>'

        # Color coding for internal strain (torsional deformation)
        strain_val = row.get("Lig_Strain", np.nan)
        if pd.isna(strain_val):
            strain_style = "color:#64748b;"
        elif strain_val <= 4.0:
            strain_style = "color:#16a34a; font-weight:600;"  # Low/Good strain (Green)
        elif strain_val <= 8.0:
            strain_style = "color:#d97706; font-weight:600;"  # Moderate strain (Amber)
        else:
            strain_style = "color:#dc2626; font-weight:700;"  # High strain (Red)
            
        # Extract best available CNN pose quality score (0.0 to 1.0)
        cnn_score = row.get("Resc_CNNscore") if pd.notna(row.get("Resc_CNNscore")) else row.get("Dock_CNNscore")

        rows_html.append(f"""
        <tr class="{row_class}">
            <td class="col-rank">{badge}</td>
            <td class="col-mol">{svg_code}</td>
            <td class="col-title"><strong>{title}</strong></td>
            <td class="col-frame"><span class="frame-tag">{row['Frame']}</span></td>
            <td class="col-score"><strong>{fmt(row.get('Ultimate_Score'))}</strong></td>
            <td class="col-delta {delta_class}">{delta_str}</td>
            <td>{rmsd_badge}</td>
            <td>{fmt(row.get('PR_Glide'))}</td>
            <td>{fmt(row.get('PR_MMGBSA'))}</td>
            <td>{fmt(row.get('PR_CNN_Resc'))}</td>
            <td>{fmt(row.get('PR_GNINA_Dock'))}</td>
            <td class="col-num">{fmt(row.get('GlideScore'))}</td>
            <td class="col-num">{fmt(row.get('MMGBSA'))}</td>
            <td class="col-num" style="{strain_style}">{fmt(strain_val)}</td>
            <td class="col-num">{fmt(row.get('Resc_CNNaffinity'))}</td>
            <td class="col-num">{fmt(row.get('Dock_CNNaffinity'))}</td>
            <td class="col-num">{fmt(cnn_score)}</td>
        </tr>
        """)

    with open(template_path, "r", encoding="utf-8") as f:
        html_template = f.read()
        
    html_content = html_template.replace("{table_body}", "\n".join(rows_html))
    
    out_file = out_dir / "00_summary_report.html"
    with open(out_file, "w", encoding="utf-8") as f:
        f.write(html_content)
        
    logger.success(f"Visual HTML report saved to: {out_file.name}")