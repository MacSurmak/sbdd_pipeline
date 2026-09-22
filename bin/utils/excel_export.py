"""
File: src/utils/excel_export.py
Description: Generates an Excel spreadsheet with embedded 2D chemical structures,
formatted table headers, cell borders, and customized scientific metrics.
"""
import io
from pathlib import Path
import pandas as pd
import numpy as np
from loguru import logger
from rdkit import Chem
from rdkit.Chem import rdDepictor
from rdkit.Chem.Draw import rdMolDraw2D

try:
    import openpyxl
    from openpyxl.styles import Font, Alignment, PatternFill, Border, Side
    from openpyxl.drawing.image import Image as OpenpyxlImage
    from openpyxl.utils import get_column_letter
    OPENPYXL_AVAILABLE = True
except ImportError:
    OPENPYXL_AVAILABLE = False


def _render_mol_png(mol, width: int = 180, height: int = 130) -> io.BytesIO:
    """Renders an RDKit Mol to a high-contrast transparent PNG in memory."""
    if mol is None:
        return None
    try:
        m = Chem.Mol(mol)
        rdDepictor.Compute2DCoords(m)
        drawer = rdMolDraw2D.MolDraw2DCairo(width, height)
        opts = drawer.drawOptions()
        opts.clearBackground = True
        opts.bondLineWidth = 2
        drawer.DrawMolecule(m)
        drawer.FinishDrawing()
        png_bytes = drawer.GetDrawingText()
        return io.BytesIO(png_bytes)
    except Exception:
        # Fallback if Cairo is not available in current RDKit build
        try:
            from rdkit.Chem import Draw
            img = Draw.MolToImage(m, size=(width, height))
            bio = io.BytesIO()
            img.save(bio, format="PNG")
            bio.seek(0)
            return bio
        except Exception:
            return None


def export_hits_to_excel(df: pd.DataFrame, cwd: Path, out_dir: Path):
    """
    Exports consensus hitlist to Excel with embedded 2D structures.
    """
    if not OPENPYXL_AVAILABLE:
        logger.warning("openpyxl is not installed. Skipping Excel report generation. Run: pip install openpyxl pillow")
        return

    if df.empty:
        return

    from utils.html_export import _load_structures_for_hits
    mols_map = _load_structures_for_hits(df, cwd)

    wb = openpyxl.Workbook()
    ws = wb.active
    ws.title = "Top Candidates"
    ws.views.sheetView[0].showGridLines = True

    # 1. Заголовки таблицы
    headers = [
        "Структура",
        "Соединение",
        "Фрейм",
        "GlideScore, ккал/моль",
        "Prime MM-GBSA, ккал/моль",
        "CNNaffinity, pKd",
        "Итоговый балл",
        "RMSD (Glide/GNINA), Å"
    ]

    # Стили
    header_fill = PatternFill(start_color="1E293B", end_color="1E293B", fill_type="solid")  # Slate 800
    header_font = Font(name="Calibri", size=11, bold=True, color="FFFFFF")
    cell_font = Font(name="Calibri", size=11)
    bold_cell_font = Font(name="Calibri", size=11, bold=True)
    
    center_align = Alignment(horizontal="center", vertical="center", wrap_text=True)
    left_align = Alignment(horizontal="left", vertical="center")
    right_align = Alignment(horizontal="right", vertical="center")

    thin_border = Border(
        left=Side(style="thin", color="CBD5E1"),
        right=Side(style="thin", color="CBD5E1"),
        top=Side(style="thin", color="CBD5E1"),
        bottom=Side(style="thin", color="CBD5E1")
    )

    # Запись строки заголовков
    ws.append(headers)
    ws.row_dimensions[1].height = 28
    for col_idx in range(1, len(headers) + 1):
        cell = ws.cell(row=1, column=col_idx)
        cell.fill = header_fill
        cell.font = header_font
        cell.alignment = center_align
        cell.border = thin_border

    # 2. Заполнение данными
    img_width, img_height = 160, 110
    row_height = 90

    for row_idx, (_, row) in enumerate(df.iterrows(), start=2):
        ws.row_dimensions[row_idx].height = row_height

        raw_title = str(row.get("Title", ""))
        is_native = bool(row.get("Is_Native", False))
        display_title = "Native Reference" if is_native else raw_title

        frame = str(row.get("Frame", ""))
        
        # Скоррелированные метрики
        gscore = row.get("GlideScore")
        mmgbsa = row.get("MMGBSA")
        
        # Выбираем лучший доступный скор CNN (Rescore или Dock)
        cnn_aff = row.get("Resc_CNNaffinity")
        if pd.isna(cnn_aff):
            cnn_aff = row.get("Dock_CNNaffinity")

        ultimate = row.get("Ultimate_Score")
        p_rmsd = row.get("Pose_RMSD")

        # Добавляем пустую ячейку под картинку в кол. 1
        ws.cell(row=row_idx, column=1, value="")
        
        c2 = ws.cell(row=row_idx, column=2, value=display_title)
        c2.font = bold_cell_font if is_native else cell_font
        c2.alignment = left_align

        c3 = ws.cell(row=row_idx, column=3, value=frame)
        c3.alignment = center_align

        c4 = ws.cell(row=row_idx, column=4, value=round(float(gscore), 2) if pd.notna(gscore) else "—")
        c4.alignment = right_align

        c5 = ws.cell(row=row_idx, column=5, value=round(float(mmgbsa), 2) if pd.notna(mmgbsa) else "—")
        c5.alignment = right_align

        c6 = ws.cell(row=row_idx, column=6, value=round(float(cnn_aff), 2) if pd.notna(cnn_aff) else "—")
        c6.alignment = right_align

        c7 = ws.cell(row=row_idx, column=7, value=round(float(ultimate), 2) if pd.notna(ultimate) else "—")
        c7.font = bold_cell_font
        c7.alignment = right_align

        c8 = ws.cell(row=row_idx, column=8, value=round(float(p_rmsd), 2) if pd.notna(p_rmsd) else "—")
        c8.alignment = right_align

        # Стили границ
        for col_idx in range(1, len(headers) + 1):
            cell = ws.cell(row=row_idx, column=col_idx)
            cell.border = thin_border
            if is_native:
                cell.fill = PatternFill(start_color="F0FDF4", end_color="F0FDF4", fill_type="solid")  # Light green for native

        # Вставка картинки молекулы
        mol = mols_map.get(raw_title)
        if mol is not None:
            bio = _render_mol_png(mol, width=img_width, height=img_height)
            if bio is not None:
                img = OpenpyxlImage(bio)
                img.width = img_width
                img.height = img_height
                # Привязываем изображение к ячейке A{row_idx}
                cell_coord = f"A{row_idx}"
                ws.add_image(img, cell_coord)

    # 3. Автоширина колонок
    col_widths = {
        "A": 24,  # Структура
        "B": 28,  # Соединение
        "C": 12,  # Фрейм
        "D": 22,  # GlideScore
        "E": 24,  # MM-GBSA
        "F": 18,  # CNNaffinity
        "G": 16,  # Итоговый балл
        "H": 22   # RMSD
    }
    for col_letter, width in col_widths.items():
        ws.column_dimensions[col_letter].width = width

    excel_path = out_dir / "00_summary_table.xlsx"
    wb.save(excel_path)
    logger.success(f"Excel report with 2D structures generated: {excel_path.name}")