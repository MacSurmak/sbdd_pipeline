"""
File: src/utils/ui.py
Description: Centralized UI drawing for rich tables.
"""
from rich.table import Table
from core.logger import console
from typing import List, Dict

def print_summary_table(title: str, columns: List[str], data: List[Dict[str, str]]):
    """
    Draws a stylized Rich table from a list of dictionaries.
    """
    if not data:
        return
        
    table = Table(title=title, header_style="bold cyan")
    
    for i, col in enumerate(columns):
        justify = "center" if i == 0 else "right"
        table.add_column(col, justify=justify, style="bold" if i == 0 else "")
        
    for row in data:
        table.add_row(*[str(row.get(col, "")) for col in columns])
        
    console.print("")
    console.print(table)
    console.print("")