"""
bin/utils/atom_selection.py
Deterministic Atom Selection Parser for MDAnalysis trajectories.
Enforces the strict physical definition of Pocket (sidechains only),
standard biophysical macros, NDX group lookups, and combinatorial set algebra.
"""

import re
from typing import Dict, List, Optional
import MDAnalysis as mda

STANDARD_BACKBONE_NAMES = {"CA", "C", "N", "O", "OXT"}


class AtomSelectionParser:
    """
    Evaluates hierarchical atom selections:
    1. Exact match in index.ndx (case-insensitive)
    2. Canonical biophysical macros (Pocket=sidechains, Backbone, Sidechains, Ligand)
    3. Combinatorial set algebra (+ for union, - for difference, & for intersection)
    4. Direct MDAnalysis selection syntax fallback
    Strict Fail-Fast on empty resolutions or syntax errors.
    """

    def __init__(
        self,
        universe: mda.Universe,
        ligand_resname: str = "",
        pocket_radius: float = 5.0,
        ndx_groups: Optional[Dict[str, List[int]]] = None,
    ):
        self.u = universe
        self.ligand_resname = ligand_resname.strip().upper() if ligand_resname else ""
        self.pocket_radius = float(pocket_radius)
        self.ndx_groups = {k.strip().lower(): v for k, v in (ndx_groups or {}).items()}

    def _resolve_primitive(self, token: str, desc: str) -> mda.AtomGroup:
        t_clean = token.strip()
        t_lower = t_clean.lower()

        # 1. Exact match in index.ndx (case-insensitive)
        if t_lower in self.ndx_groups:
            atom_indices = self.ndx_groups[t_lower]
            ag = self.u.atoms[atom_indices]
            if len(ag) == 0:
                raise ValueError(
                    f"FATAL: Group '[ {t_clean} ]' from index.ndx resolved to 0 atoms for {desc}!"
                )
            return ag

        # 2. Canonical Biophysical Macros
        if t_lower in ["backbone", "bb"]:
            ag = self.u.select_atoms(
                "protein and (name CA or name C or name N or name O or name OXT)"
            )
            if len(ag) == 0:
                raise ValueError(
                    f"FATAL: Macro 'Backbone' resolved to 0 atoms for {desc}!"
                )
            return ag

        elif t_lower in ["c-alpha", "ca"]:
            ag = self.u.select_atoms("protein and name CA")
            if len(ag) == 0:
                raise ValueError(
                    f"FATAL: Macro 'C-alpha' resolved to 0 atoms for {desc}!"
                )
            return ag

        elif t_lower in ["sidechains", "sc"]:
            ag = self.u.select_atoms(
                "protein and not (name CA or name C or name N or name O or name OXT) and not name H* and not type H"
            )
            if len(ag) == 0:
                raise ValueError(
                    f"FATAL: Macro 'Sidechains' resolved to 0 atoms for {desc}!"
                )
            return ag

        elif t_lower in ["ligand", "lig"]:
            if not self.ligand_resname:
                raise ValueError(
                    f"FATAL: Macro 'Ligand' requires --ligand_resname for {desc}, but none was provided!"
                )
            ag = self.u.select_atoms(
                f"resname {self.ligand_resname} and not name H* and not type H"
            )
            if len(ag) == 0:
                raise ValueError(
                    f"FATAL: Macro 'Ligand' (resname {self.ligand_resname}) resolved to 0 atoms in {desc}!"
                )
            return ag

        elif t_lower in ["pocket"]:
            # STRICT PHYSICAL DEFINITION: Pocket consists EXCLUSIVELY of SIDECHAINS of residues within pocket_radius of ligand!
            if not self.ligand_resname:
                raise ValueError(
                    f"FATAL: Macro 'Pocket' requires --ligand_resname for {desc}, but none was provided!"
                )
            lig_atoms = self.u.select_atoms(
                f"resname {self.ligand_resname} and not name H* and not type H"
            )
            if len(lig_atoms) == 0:
                raise ValueError(
                    f"FATAL: Reference ligand '{self.ligand_resname}' not found for Pocket selection in {desc}!"
                )
            ag = self.u.select_atoms(
                f"(same residue as (protein and around {self.pocket_radius} resname {self.ligand_resname})) "
                f"and not (name CA or name C or name N or name O or name OXT) and not name H* and not type H"
            )
            if len(ag) == 0:
                raise ValueError(
                    f"FATAL: Macro 'Pocket' (sidechains around {self.pocket_radius} Å of {self.ligand_resname}) "
                    f"resolved to 0 atoms in {desc}!"
                )
            return ag

        elif t_lower in ["pocket_all", "pocket_full"]:
            # Both backbone and sidechains within pocket
            if not self.ligand_resname:
                raise ValueError(
                    f"FATAL: Macro 'Pocket_All' requires --ligand_resname for {desc}!"
                )
            ag = self.u.select_atoms(
                f"(same residue as (protein and around {self.pocket_radius} resname {self.ligand_resname})) "
                f"and not name H* and not type H"
            )
            if len(ag) == 0:
                raise ValueError(
                    f"FATAL: Macro 'Pocket_All' resolved to 0 atoms in {desc}!"
                )
            return ag

        # 3. Direct MDAnalysis Selection Syntax Fallback
        try:
            ag = self.u.select_atoms(t_clean)
            if len(ag) > 0:
                return ag
        except Exception as e:
            raise ValueError(
                f"FATAL: Selection expression '{t_clean}' failed to parse for {desc}: {e}"
            )

        avail_ndx = list(self.ndx_groups.keys())
        raise ValueError(
            f"FATAL: Selection token '{t_clean}' for {desc} resolved to 0 atoms!\n"
            f"Available NDX groups: {avail_ndx}\n"
            f"Recognized macros: Backbone, C-alpha, Sidechains, Ligand, Pocket, Pocket_All"
        )

    def parse(self, query: str, desc: str = "Atom Selection") -> mda.AtomGroup:
        """
        Parses composite query with set algebra (+ for union, - for difference, & for intersection).
        """
        if not query or not query.strip():
            raise ValueError(
                f"FATAL: Empty atom selection query provided for {desc}!"
            )

        raw_query = query.strip()
        tokens = re.split(r"(\+|\-|\&)", raw_query)

        if len(tokens) == 1:
            ag = self._resolve_primitive(tokens[0], desc)
        else:
            ag = self._resolve_primitive(tokens[0], desc)
            i = 1
            while i < len(tokens):
                op = tokens[i].strip()
                rhs_token = tokens[i + 1].strip()
                rhs_ag = self._resolve_primitive(rhs_token, desc)

                if op == "+":
                    ag = ag | rhs_ag
                elif op == "-":
                    ag = ag - rhs_ag
                elif op == "&":
                    ag = ag & rhs_ag
                else:
                    raise ValueError(
                        f"FATAL: Unknown set operator '{op}' in query '{raw_query}'!"
                    )
                i += 2

        if len(ag) == 0:
            raise ValueError(
                f"FATAL: Evaluated selection query '{raw_query}' for {desc} resulted in 0 atoms!"
            )

        return ag
