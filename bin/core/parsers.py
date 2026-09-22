"""
File: src/core/parsers.py
Description: Modular strategy-based parsers for external CLI tool logs.
Standardizes progress updates for Rich UI and isolates tool-specific regex.
"""
import re
from dataclasses import dataclass
from typing import Optional, Any

@dataclass
class ProgressUpdate:
    completed: Optional[int] = None
    total: Optional[int] = None
    description: Optional[str] = None

class BaseLogParser:
    def parse(self, text: str) -> Optional[ProgressUpdate]:
        raise NotImplementedError

class PackmolParser(BaseLogParser):
    def __init__(self):
        self.progress_re = re.compile(r"(\d+)/(\d+)\s*\[")
        self.pct_re = re.compile(r"(\d+)%")

    def parse(self, text: str) -> Optional[ProgressUpdate]:
        lines = [line.strip() for line in text.replace("\r", "\n").split("\n") if line.strip()]
        for line in reversed(lines[-30:]):
            m = self.progress_re.search(line)
            if m:
                comp, tot = int(m.group(1)), int(m.group(2))
                pct_m = self.pct_re.search(line)
                pct = pct_m.group(1) if pct_m else str(int(comp / tot * 100))
                stage = line.split(":")[0].strip() if ":" in line else "Packing"
                return ProgressUpdate(
                    completed=comp,
                    total=tot,
                    description=f"[bold magenta]Packmol: {stage} ({pct}%)[/bold magenta]"
                )
        return None

class GninaParser(BaseLogParser):
    def __init__(self):
        self.header_re = re.compile(r"mode\s+\|\s+affinity", re.IGNORECASE)

    def parse(self, text: str) -> Optional[ProgressUpdate]:
        matches = len(self.header_re.findall(text))
        if matches > 0:
            return ProgressUpdate(completed=matches)
        return None

class DeclashParser(BaseLogParser):
    def __init__(self):
        # Matches iteration and remaining clashes: "iter  12: clashes=   45 applied= ..."
        self.iter_re = re.compile(r"iter\s+(\d+):\s+clashes=\s*(\d+)", re.IGNORECASE)

    def parse(self, text: str) -> Optional[ProgressUpdate]:
        lines = [l.strip() for l in text.splitlines() if l.strip()]
        for line in reversed(lines[-20:]):
            m = self.iter_re.search(line)
            if m:
                it = int(m.group(1)) + 1
                clashes = int(m.group(2))
                return ProgressUpdate(
                    completed=it,
                    total=300,
                    description=f"[bold cyan]Declash: Step {it}/300 | Clashes left: {clashes}[/bold cyan]"
                )
        return None

class GlideParser(BaseLogParser):
    def __init__(self):
        self.lig_re = re.compile(r"LIG\s+(\d+)")

    def parse(self, text: str) -> Optional[ProgressUpdate]:
        lines = [l.strip() for l in text.splitlines() if l.strip()]
        for line in reversed(lines[-20:]):
            m = self.lig_re.search(line)
            if m:
                return ProgressUpdate(completed=int(m.group(1)))
        return None

class PrimeParser(BaseLogParser):
    def __init__(self):
        # Matches JobDJ table rows: "  8  8  0 | finished ..."
        self.jobdj_re = re.compile(r"^\s*(\d+)\s+(\d+)\s+(\d+)\s+\|\s+(finished|launched)", re.MULTILINE)

    def parse(self, text: str) -> Optional[ProgressUpdate]:
        lines = [l for l in text.splitlines() if "|" in l]
        for line in reversed(lines[-20:]):
            m = self.jobdj_re.search(line)
            if m:
                c, a, w = int(m.group(1)), int(m.group(2)), int(m.group(3))
                tot = c + a + w
                if tot > 0:
                    pct = int((c / tot) * 100)
                    return ProgressUpdate(
                        completed=c,
                        total=tot,
                        description=f"[yellow]Prime MM-GBSA Subjobs: {c}/{tot} ({pct}%)[/yellow]"
                    )
        return None

class GromacsParser(BaseLogParser):
    def __init__(self):
        # Matches inline steps: "step 220000", "Step= 500", "Writing checkpoint, step 220000"
        self.inline_re = re.compile(r"(?:step\s*[:=]?\s*|writing checkpoint,\s+step\s+)(\d+)", re.IGNORECASE)
        # Matches GROMACS energy table header: "Step           Time"
        self.header_re = re.compile(r"^\s*Step\s+Time\s*$", re.IGNORECASE)
        # Matches step and time values: "220000      440.00000"
        self.step_val_re = re.compile(r"^\s*(\d+)\s+[\d\.]+\s*$")

    def parse(self, text: str) -> Optional[ProgressUpdate]:
        lines = [l.strip() for l in text.splitlines() if l.strip()]
        for idx in range(len(lines) - 1, -1, -1):
            line = lines[idx]

            # 1. Parse two-line step block from .log files
            m_val = self.step_val_re.match(line)
            if m_val and idx > 0 and self.header_re.match(lines[idx - 1]):
                return ProgressUpdate(completed=int(m_val.group(1)))

            # 2. Parse minimization live metrics: "Step= 420, Dmax= ..., Epot= -9.8e+05 Fmax= 8.4e+02"
            if "fmax=" in line.lower() and "step=" in line.lower():
                try:
                    parts = line.split(",")
                    step_val = int(parts[0].split("=")[1].strip())
                    epot_val = parts[2].split("=")[1].split()[0].strip()
                    fmax_val = float(parts[2].split("Fmax=")[1].split()[0].strip())
                    return ProgressUpdate(
                        completed=step_val,
                        description=f"[bold green]EM: Step {step_val} | Fmax: {fmax_val:.1f} -> 500 | Epot: {epot_val}[/bold green]"
                    )
                except Exception:
                    pass

            # 3. Parse inline MD progress
            m_inline = self.inline_re.search(line)
            if m_inline:
                return ProgressUpdate(completed=int(m_inline.group(1)))

        return None

class MmpbsaParser(BaseLogParser):
    def __init__(self):
        self.tqdm_re = re.compile(r"(\d+)%\|.*?\|\s*(\d+)/(\d+)")
        self.stage_re = re.compile(r"calculating\s+(complex|receptor|ligand)\s+contribution", re.IGNORECASE)
        self.model_re = re.compile(r"Beginning\s+(GB|PB)\s+calculations", re.IGNORECASE)
        self.current_model = "PB"
        self.substage_name = "Complex"
        self.substage_idx = 1

    def parse(self, text: str) -> Optional[ProgressUpdate]:
        # Scan entire buffer forward so we never lose the state
        for line in text.splitlines():
            m_mod = self.model_re.search(line)
            if m_mod:
                self.current_model = m_mod.group(1).upper()
                
            m_stg = self.stage_re.search(line)
            if m_stg:
                name = m_stg.group(1).capitalize()
                self.substage_name = name
                self.substage_idx = 1 if name == "Complex" else (2 if name == "Receptor" else 3)

        # Parse latest frame progress from the tail
        lines = [l.strip() for l in text.splitlines() if l.strip()]
        for line in reversed(lines[-20:]):
            m = self.tqdm_re.search(line)
            if m:
                pct = int(m.group(1))
                c = int(m.group(2))
                tot = int(m.group(3))
                return ProgressUpdate(
                    completed=c,
                    total=tot,
                    description=f"[bold green]MMPBSA-{self.current_model} [{self.substage_name} {self.substage_idx}/3]: {c}/{tot} ({pct}%)[/bold green]"
                )
        return None

def get_parser(cmd: list) -> Optional[BaseLogParser]:
    """Factory function to auto-detect the parser based on command context."""
    cmd_str = " ".join(str(x) for x in cmd)
    cmd_name = cmd[0].split("/")[-1]

    if "gmx_MMPBSA" in cmd_name or "gmx_MMPBSA" in cmd_str:
        return MmpbsaParser()
    elif cmd_name == "packmol-memgen":
        return PackmolParser()
    elif cmd_name == "gnina":
        return GninaParser()
    elif "declash" in cmd_str:
        return DeclashParser()
    elif cmd_name.startswith("gmx"):
        return GromacsParser()
    elif cmd_name in ["glide", "ligprep"] or any(str(x).endswith(".in") for x in cmd):
        return GlideParser()
    elif "prime" in cmd_name or "prime_mmgbsa" in cmd_str:
        return PrimeParser()
    return None