"""
File: cli/commands/add.py
Description: Clean 3NF CLI handlers for importing targets, structural frames, ligands, and trajectories.
Enforces InChIKey deduplication, atomic poses, and multi-model frame transparent splitting.
"""

from pathlib import Path
import sys
import shutil
from rich.console import Console
from cli.core.db import (
    init_db,
    Target,
    ReceptorFrame,
    Substance,
    SubstanceState,
    Pose,
    System,
    Trajectory,
)
from cli.core.chem import (
    _format_stereo_tag,
    get_microstate_identity,
    get_mol_identity,
)

console = Console()


import math


def _detect_gaps(struct_path: Path) -> bool:
    """Rigorous gap detection supporting both PDB and mmCIF via Bio.PDB. Strictly forbids silent fallbacks."""
    from Bio.PDB import PDBParser, MMCIFParser, is_aa

    suffix = struct_path.suffix.lower()
    parser = (
        MMCIFParser(QUIET=True)
        if suffix in [".cif", ".mmcif"]
        else PDBParser(QUIET=True)
    )
    st = parser.get_structure("struct", str(struct_path))

    for model in st:
        for chain in model:
            prev_ca = None
            prev_res = None
            for res in chain:
                if is_aa(res) and "CA" in res:
                    res_id = res.get_id()
                    res_seq = res_id[1]
                    coord = res["CA"].get_coord()
                    if prev_ca is not None and prev_res is not None:
                        d2 = sum((coord[k] - prev_ca[k]) ** 2 for k in range(3))
                        if (res_seq - prev_res > 1) or d2 > 20.25:
                            return True
                    prev_ca = coord
                    prev_res = res_seq
    return False


STANDARD_PROTEIN_RESIDUES = {
    "ALA",
    "ARG",
    "ASN",
    "ASP",
    "CYS",
    "CYX",
    "GLN",
    "GLU",
    "GLY",
    "HIS",
    "HID",
    "HIE",
    "HIP",
    "ILE",
    "LEU",
    "LYS",
    "LYN",
    "MET",
    "PHE",
    "PRO",
    "SER",
    "THR",
    "TRP",
    "TYR",
    "VAL",
    "ASH",
    "GLH",
    "ACE",
    "NME",
    "NMA",
    "MSE",
    "CSE",
}


class CleanStructureSelect:
    """Strict filter: retains standard protein residues AND explicitly declared heterogens only."""

    def __init__(
        self,
        allowed_ligand: str | None,
        allowed_cofactors: list[str],
        allowed_ions: list[str],
    ):
        self.allowed_resnames = set()
        if allowed_ligand:
            self.allowed_resnames.add(allowed_ligand.strip().upper())
        for c in allowed_cofactors:
            self.allowed_resnames.add(c.strip().upper())
        for i in allowed_ions:
            self.allowed_resnames.add(i.strip().upper())

    def accept_model(self, model):
        return 1

    def accept_chain(self, chain):
        return 1

    def accept_atom(self, atom):
        return 1

    def accept_residue(self, res):
        from Bio.PDB import is_aa

        r_name = res.get_resname().strip().upper()
        if r_name in STANDARD_PROTEIN_RESIDUES or (
            is_aa(res) and res.get_id()[0] == " "
        ):
            return 1
        if r_name in self.allowed_resnames:
            return 1
        return 0


def _sanitize_and_split_models(
    struct_path: Path,
    dest_dir: Path,
    base_name: str,
    allowed_ligand: str | None,
    allowed_cofactors: list[str],
    allowed_ions: list[str],
) -> list[tuple[str, Path]]:
    """
    Parses macromolecular structure, executes strict Fail-Fast verification of explicitly
    requested ligand, cofactors, and ions, and produces sanitized complex frames in dedicated
    frame subdirectories: targets/<target>/frames/<frame>/<frame>_raw.pdb.
    """
    from Bio.PDB import PDBParser, MMCIFParser, PDBIO

    suffix = struct_path.suffix.lower()
    if suffix in [".cif", ".mmcif"]:
        parser = MMCIFParser(QUIET=True)
    elif suffix == ".pdb":
        parser = PDBParser(QUIET=True)
    else:
        raise ValueError(
            f"FATAL: Unsupported structure format '{suffix}'. Expected .pdb or .cif/.mmcif"
        )

    st = parser.get_structure("struct", str(struct_path))

    # Strict Fail-Fast verification of explicitly requested molecules
    all_resnames = set()
    for model in st:
        for chain in model:
            for res in chain:
                all_resnames.add(res.get_resname().strip().upper())

    if allowed_ligand and allowed_ligand.upper() not in all_resnames:
        raise ValueError(
            f"FATAL: Explicitly specified ligand '{allowed_ligand}' was not found in {struct_path.name}! "
            f"Available residues: {sorted(list(all_resnames))}"
        )

    for cof in allowed_cofactors:
        if cof.upper() not in all_resnames:
            raise ValueError(
                f"FATAL: Explicitly specified cofactor '{cof}' was not found in {struct_path.name}! "
                f"Available residues: {sorted(list(all_resnames))}"
            )

    for ion in allowed_ions:
        if ion.upper() not in all_resnames:
            raise ValueError(
                f"FATAL: Explicitly specified ion '{ion}' was not found in {struct_path.name}! "
                f"Available residues: {sorted(list(all_resnames))}"
            )

    io = PDBIO()
    clean_selector = CleanStructureSelect(
        allowed_ligand, allowed_cofactors, allowed_ions
    )

    models = list(st.get_models())
    results = []

    if len(models) > 1:
        console.print(
            f"[cyan]Detected {len(models)} models in ensemble. Sanitizing complete complex frames...[/cyan]"
        )
        for idx, model in enumerate(models, start=1):
            f_name = f"{base_name}_m{idx}"
            frame_dir = dest_dir / f_name
            frame_dir.mkdir(parents=True, exist_ok=True)
            out_complex = frame_dir / f"{f_name}_raw.pdb"

            io.set_structure(model)
            io.save(str(out_complex), clean_selector)
            results.append((f_name, out_complex))
    else:
        f_name = base_name
        frame_dir = dest_dir / f_name
        frame_dir.mkdir(parents=True, exist_ok=True)
        out_complex = frame_dir / f"{f_name}_raw.pdb"

        io.set_structure(st[0])
        io.save(str(out_complex), clean_selector)
        results.append((f_name, out_complex))

    return results


def _extract_ligand_atoms(rec_path: Path, lig_resname: str, out_sdf: Path) -> None:
    """
    Extracts ligand residue from PDB, CIF, MMCIF or MAE into a standalone SDF.
    Perceives bond orders and aromaticity chemically via Schrodinger assignbondorders
    so that identical molecules produce the exact same topological InChIKey.
    """
    import os
    import subprocess
    import shutil

    lig_clean = lig_resname.strip().upper()
    out_sdf.parent.mkdir(parents=True, exist_ok=True)

    schrod_bin = shutil.which("schrodinger")
    if not schrod_bin and "SCHRODINGER" in os.environ:
        cand = Path(os.environ["SCHRODINGER"]) / "run"
        if cand.exists():
            schrod_bin = str(cand)

    if not schrod_bin:
        raise RuntimeError(
            "FATAL: Schrödinger runner ('schrodinger run') is required for robust bond order perception."
        )

    schrod_cmd = (
        [schrod_bin, "run", "python3"]
        if not schrod_bin.endswith("/run")
        else [schrod_bin, "python3"]
    )

    script_path = (
        Path(__file__).parent.parent.parent
        / "bin"
        / "utils"
        / "extract_ligand_bonds.py"
    )
    if not script_path.exists():
        raise FileNotFoundError(f"FATAL: Extraction script not found at {script_path}")

    cmd = schrod_cmd + [
        str(script_path.resolve()),
        str(rec_path.resolve()),
        lig_clean,
        str(out_sdf.resolve()),
    ]
    res = subprocess.run(cmd, capture_output=True, text=True)

    if res.returncode == 2:
        raise ValueError(
            f"FATAL: Ligand residue '{lig_clean}' not found in {rec_path.name}! {res.stderr.strip()}"
        )
    elif res.returncode != 0:
        raise RuntimeError(
            f"FATAL: Extraction failed with exit code {res.returncode}:\n{res.stderr}\n{res.stdout}"
        )

    if not out_sdf.exists() or out_sdf.stat().st_size == 0:
        raise ValueError(
            f"FATAL: Ligand extraction failed to produce valid SDF at {out_sdf}"
        )


from Bio.PDB.Polypeptide import protein_letters_3to1


def _extract_resolved_sequences_by_chain(st) -> dict[str, str]:
    """Extracts 1-letter amino acid sequences from resolved protein chains in structure."""
    from Bio.PDB import is_aa

    chain_seqs = {}
    for model in st:
        for chain in model:
            letters = []
            for res in chain:
                if is_aa(res) and "CA" in res:
                    rname = res.get_resname().strip().upper()
                    letters.append(protein_letters_3to1.get(rname, "X"))
            if letters:
                chain_seqs[chain.id] = "".join(letters)
        break  # Use first model for primary sequence definition
    return chain_seqs


def _build_multi_fasta(chain_dict: dict[str, str]) -> str:
    """Encodes dictionary of chain sequences into canonical multi-FASTA string."""
    return "\n".join(f">{c}\n{s}" for c, s in sorted(chain_dict.items()))


def _align_and_clamp_sequence(struct_seq: str, full_fasta_seq: str) -> str:
    """
    Performs pairwise alignment to map an isolated construct (e.g. residues 600-800)
    against full-length gene FASTA and clamps target sequence to construct boundaries.
    """
    from Bio import Align

    aligner = Align.PairwiseAligner()
    aligner.mode = "global"
    aligner.match_score = 2
    aligner.mismatch_score = -1
    aligner.open_gap_score = -4
    aligner.extend_gap_score = -0.5

    alignments = aligner.align(full_fasta_seq, struct_seq)
    if not alignments:
        raise ValueError(
            "FATAL: Provided FASTA sequence cannot be aligned to experimental structure!"
        )

    best = alignments[0]
    blocks_fasta = best.aligned[0]
    f_start = blocks_fasta[0][0]
    f_end = blocks_fasta[-1][1]

    # Return clamped sequence corresponding strictly to resolved construct domain
    return full_fasta_seq[f_start:f_end]


def _resolve_target_sequence(
    target_name: str, struct_path: Path, fasta_arg: str | None, has_gaps: bool
) -> str | None:
    """
    Deterministically resolves reference sequence without blocking interactive prompts.
    - If continuous (no gaps): automatically extracts sequence directly from 3D coordinates.
    - If gaps exist and FASTA provided: aligns and clamps reference construct sequence.
    - If gaps exist and NO FASTA: logs warning and sets sequence to None (awaiting sbdd heal).
    """
    from Bio.PDB import PDBParser, MMCIFParser

    suffix = struct_path.suffix.lower()
    parser = (
        MMCIFParser(QUIET=True)
        if suffix in [".cif", ".mmcif"]
        else PDBParser(QUIET=True)
    )
    st = parser.get_structure("struct", str(struct_path))
    resolved_chains = _extract_resolved_sequences_by_chain(st)

    if not resolved_chains:
        raise ValueError(f"FATAL: No amino acid chains found in {struct_path.name}!")

    if not fasta_arg:
        if not has_gaps:
            # Continuous chain(s): 100% deterministic sequence extraction directly from atoms
            return _build_multi_fasta(resolved_chains)
        else:
            # Broken chain with gaps: coordinates cannot guess missing amino acids
            return None

    # FASTA file explicitly provided by user
    f_path = Path(fasta_arg).resolve()
    if not f_path.exists():
        raise FileNotFoundError(f"FATAL: Specified FASTA file not found: {f_path}")

    raw_fasta = f_path.read_text(encoding="utf-8", errors="ignore")
    fasta_records = {}
    curr_id = None
    curr_seq = []

    for line in raw_fasta.splitlines():
        line = line.strip()
        if not line:
            continue
        if line.startswith(">"):
            if curr_id and curr_seq:
                fasta_records[curr_id] = "".join(curr_seq).upper()
            curr_id = line[1:].split()[0].strip()
            curr_seq = []
        else:
            curr_seq.append(line)
    if curr_id and curr_seq:
        fasta_records[curr_id] = "".join(curr_seq).upper()

    if not fasta_records:
        raise ValueError(f"FATAL: FASTA file {f_path} is empty or invalid.")

    # Match and clamp chains
    clamped_chains = {}
    if len(resolved_chains) == 1 and len(fasta_records) == 1:
        c_id = list(resolved_chains.keys())[0]
        full_seq = list(fasta_records.values())[0]
        clamped_chains[c_id] = _align_and_clamp_sequence(
            resolved_chains[c_id], full_seq
        )
    else:
        for c_id, r_seq in resolved_chains.items():
            if c_id in fasta_records:
                clamped_chains[c_id] = _align_and_clamp_sequence(
                    r_seq, fasta_records[c_id]
                )
            elif len(fasta_records) == 1:
                full_seq = list(fasta_records.values())[0]
                clamped_chains[c_id] = _align_and_clamp_sequence(r_seq, full_seq)
            else:
                raise KeyError(
                    f"FATAL: Chain '{c_id}' in structure has no matching header in multi-FASTA {f_path.name}! "
                    f"Available headers: {list(fasta_records.keys())}"
                )

    return _build_multi_fasta(clamped_chains)


def add_receptor(
    paths: list[str],
    target: str | None = None,
    name: str | None = None,
    cofactors: list[str] | None = None,
    ions: list[str] | None = None,
    ligand: str | None = None,
    fasta: str | None = None,
) -> None:
    """
    Imports macromolecular structure(s) (PDB/CIF) into a Target with mandatory reference sequence.
    Enforces 100% APO protein frames and extracts explicit heterogens into atomic 3D poses.
    """
    session = init_db()
    cwd = Path.cwd()
    target_cli_name = target

    raw_inputs = paths if isinstance(paths, list) else [paths]
    struct_paths = [Path(p).resolve() for p in raw_inputs]

    for struct_path in struct_paths:
        if not struct_path.exists():
            console.print(
                f"[bold red]FATAL: Structure file not found: {struct_path}[/bold red]"
            )
            sys.exit(1)

        target_name = target_cli_name if target_cli_name else struct_path.stem
        target_obj = session.query(Target).filter_by(name=target_name).first()

        base_name = name if (len(struct_paths) == 1 and name) else struct_path.stem
        safe_cofactors = [c.strip().upper() for c in (cofactors or [])]
        safe_ions = [i.strip().upper() for i in (ions or [])]
        lig_code = ligand.strip().upper() if ligand else None
        source_type = "crystal"

        frames_dir = cwd / "targets" / target_name / "frames"
        frames = _sanitize_and_split_models(
            struct_path, frames_dir, base_name, lig_code, safe_cofactors, safe_ions
        )

        # Check overall gap status across models
        ensemble_has_gaps = any(
            _detect_gaps(complex_file) for _, complex_file in frames
        )

        target_dir = cwd / "targets" / target_name
        target_dir.mkdir(parents=True, exist_ok=True)

        if not target_obj:
            ref_seq = _resolve_target_sequence(
                target_name, struct_path, fasta, ensemble_has_gaps
            )
            target_obj = Target(
                name=target_name,
                reference_sequence=ref_seq,
                description=f"Target {target_name}",
            )
            session.add(target_obj)
            session.flush()

            if ref_seq:
                fasta_file = target_dir / f"{target_name}_sequence.fasta"
                fasta_file.write_text(ref_seq + "\n", encoding="utf-8")

            seq_info = (
                f"{len(ref_seq.splitlines()[-1])} aa"
                if ref_seq
                else "[yellow]pending FASTA[/yellow]"
            )
            console.print(
                f"[bold green]✓ Created Target '{target_name}'[/bold green] [dim](Sequence: {seq_info})[/dim]"
            )
        elif not target_obj.reference_sequence:
            ref_seq = _resolve_target_sequence(
                target_name, struct_path, fasta, ensemble_has_gaps
            )
            if ref_seq:
                target_obj.reference_sequence = ref_seq
                fasta_file = target_dir / f"{target_name}_sequence.fasta"
                fasta_file.write_text(ref_seq + "\n", encoding="utf-8")
                session.flush()
                console.print(
                    f"[bold green]✓ Associated Sequence with existing Target '{target_name}'[/bold green]"
                )

        for frame_name, complex_file in frames:
            has_gaps = _detect_gaps(complex_file)
            initial_state = "broken" if has_gaps else "raw"

            ex_frame = session.query(ReceptorFrame).filter_by(name=frame_name).first()
            if ex_frame:
                ex_frame.path = str(complex_file.resolve())
                ex_frame.state = initial_state
                frame_obj = ex_frame
                console.print(
                    f"[yellow]Updated Frame '{frame_name}' in Target '{target_name}' (state: {initial_state})[/yellow]"
                )
            else:
                frame_obj = ReceptorFrame(
                    target_id=target_obj.id,
                    name=frame_name,
                    source_type=source_type,
                    state=initial_state,
                    path=str(complex_file.resolve()),
                )
                session.add(frame_obj)
                session.flush()

                status_color = "red" if has_gaps else "green"
                console.print(
                    f"[bold green]✓ Added Frame '{frame_name}'[/bold green] [dim]([{status_color}]{initial_state}[/{status_color}], {source_type})[/dim]"
                )

            if has_gaps:
                console.print(
                    f"[!] [bold yellow]Warning: Internal missing loops detected in '{frame_name}'.[/bold yellow]\n"
                    f"    It cannot be processed by 'sbdd prep'. Provide FASTA and run 'sbdd heal' to physically repair it."
                )

            # Enumerate all physical residue instances matching resname in complex_file
            from Bio.PDB import PDBParser as _BioPDBParser

            _st_parsed = _BioPDBParser(QUIET=True).get_structure("frame", str(complex_file))

            def _process_heterogen_instances(resname: str, sub_type: str) -> list[Pose]:
                import tempfile

                clean_res_name = resname.strip().upper()
                matching_instances = []
                for _chain in _st_parsed[0]:
                    for _res in _chain:
                        if _res.get_resname().strip().upper() == clean_res_name:
                            matching_instances.append((_chain.id.strip(), int(_res.get_id()[1])))

                if not matching_instances:
                    raise ValueError(
                        f"FATAL: Residue '{clean_res_name}' not found in sanitized frame {complex_file}!"
                    )

                created_poses = []
                multi_instance = len(matching_instances) > 1

                for c_id, r_num in matching_instances:
                    with tempfile.NamedTemporaryFile(suffix=".sdf", delete=False) as tmp_f:
                        tmp_extract_sdf = Path(tmp_f.name)

                    spec_token = f"{c_id}:{r_num}:{clean_res_name}"
                    _extract_ligand_atoms(complex_file, spec_token, tmp_extract_sdf)
                    (
                        parent_inchikey,
                        parent_smiles,
                        mol_name,
                        state_hash,
                        isomeric_smiles,
                        formal_charge,
                    ) = get_microstate_identity(tmp_extract_sdf)

                    sub_dir = cwd / "substances" / parent_inchikey
                    sub_dir.mkdir(parents=True, exist_ok=True)

                    base_sdf = sub_dir / f"{clean_res_name}_base.sdf"
                    if not base_sdf.exists():
                        shutil.copy2(tmp_extract_sdf, base_sdf)

                    sub_entity = (
                        session.query(Substance).filter_by(inchikey=parent_inchikey).first()
                    )
                    if not sub_entity:
                        sub_entity = Substance(
                            inchikey=parent_inchikey,
                            smiles=parent_smiles,
                            name=mol_name or clean_res_name,
                            res_name=clean_res_name,
                            substance_type=sub_type,
                            state="raw",
                            batch_id="crystal_extract",
                            path_sdf=str(base_sdf.resolve()),
                        )
                        session.add(sub_entity)
                        session.flush()
                        console.print(
                            f"[bold green]✓ Registered {sub_type.capitalize()} '{sub_entity.name}'[/bold green] [dim](raw, {parent_inchikey})[/dim]"
                        )

                    state_entity = (
                        session.query(SubstanceState)
                        .filter_by(substance_id=sub_entity.id, state_hash=state_hash)
                        .first()
                    )
                    if not state_entity:
                        state_entity = SubstanceState(
                            substance_id=sub_entity.id,
                            state_hash=state_hash,
                            isomeric_smiles=isomeric_smiles,
                            formal_charge=formal_charge,
                            state="raw",
                            path_sdf=str(base_sdf.resolve()),
                        )
                        session.add(state_entity)
                        session.flush()

                    if multi_instance:
                        pose_name = f"{frame_name}_{clean_res_name}_{c_id}{r_num}_crystal"
                    else:
                        pose_name = f"{frame_name}_{clean_res_name}_crystal"

                    pose_dir = sub_dir / "poses" / pose_name
                    pose_dir.mkdir(parents=True, exist_ok=True)
                    pose_sdf = pose_dir / f"{clean_res_name}_pose.sdf"

                    shutil.move(str(tmp_extract_sdf), str(pose_sdf))
                    if tmp_extract_sdf.exists():
                        tmp_extract_sdf.unlink()

                    meta_scores = {"chain_id": c_id, "resnum": r_num}
                    p_obj = session.query(Pose).filter_by(name=pose_name).first()
                    if p_obj:
                        p_obj.substance_state_id = state_entity.id
                        p_obj.path = str(pose_sdf.resolve())
                        p_obj.path_sdf = str(pose_sdf.resolve())
                        p_obj.res_name = clean_res_name
                        p_obj.scores_dict = meta_scores
                    else:
                        p_obj = Pose(
                            substance_state_id=state_entity.id,
                            receptor_frame_id=frame_obj.id,
                            name=pose_name,
                            res_name=clean_res_name,
                            origin="crystal",
                            state="ready",
                            path=str(pose_sdf.resolve()),
                            path_sdf=str(pose_sdf.resolve()),
                        )
                        p_obj.scores_dict = meta_scores
                        session.add(p_obj)
                        session.flush()
                        ch_label = (
                            f"q={formal_charge:+d}" if formal_charge != 0 else "neutral"
                        )
                        console.print(
                            f"[bold green]✓ Registered Pose '{pose_name}'[/bold green] [dim]({clean_res_name} {c_id}:{r_num}, {ch_label}, state={state_hash})[/dim]"
                        )
                    created_poses.append(p_obj)
                return created_poses

            main_ligand_pose = None
            if lig_code:
                lig_poses = _process_heterogen_instances(lig_code, "ligand")
                if len(lig_poses) > 1:
                    raise ValueError(
                        f"FATAL: Active ligand '{lig_code}' has {len(lig_poses)} copies in frame '{frame_name}'! "
                        "Ensure the structure contains a single primary chain or unique ligand residue ID."
                    )
                main_ligand_pose = lig_poses[0]

            hetero_poses = []
            for cof in safe_cofactors:
                hetero_poses.extend(_process_heterogen_instances(cof, "cofactor"))

            for ion in safe_ions:
                hetero_poses.extend(_process_heterogen_instances(ion, "ion"))

            # System assembly (Pure relational binding of Frame + Poses)
            sys_name = f"{frame_name}_{lig_code}" if lig_code else frame_name
            ex_sys = session.query(System).filter_by(name=sys_name).first()
            if ex_sys:
                ex_sys.pose_id = main_ligand_pose.id if main_ligand_pose else None
                ex_sys.state = "raw"
                ex_sys.cofactors = hetero_poses
            else:
                new_sys = System(
                    target_id=target_obj.id,
                    receptor_frame_id=frame_obj.id,
                    pose_id=main_ligand_pose.id if main_ligand_pose else None,
                    name=sys_name,
                    state="raw",
                    origin="crystal",
                )
                new_sys.cofactors = hetero_poses
                session.add(new_sys)
                sys_type = "Holo" if main_ligand_pose else "Apo"
                console.print(
                    f"[bold green]✓ Created {sys_type} System '{sys_name}'[/bold green] [dim]with {len(hetero_poses)} cofactor(s)/ion(s)[/dim]"
                )

    session.commit()
    session.close()


def add_fasta(path: str) -> None:
    """Imports FASTA reference sequences for Targets."""
    session = init_db()
    fasta_path = Path(path).resolve()

    if not fasta_path.exists():
        console.print(f"[bold red]FATAL: FASTA file not found: {fasta_path}[/bold red]")
        sys.exit(1)

    updated = 0
    with open(fasta_path, "r", encoding="utf-8", errors="ignore") as f:
        current_target = None
        seq_lines = []

        def save_current():
            nonlocal updated
            if current_target and seq_lines:
                t = session.query(Target).filter_by(name=current_target).first()
                if t:
                    t.reference_sequence = "".join(seq_lines).replace(" ", "").upper()
                    updated += 1
                else:
                    console.print(
                        f"[yellow]Warning: Target '{current_target}' not found in DB, skipping sequence.[/yellow]"
                    )

        for line in f:
            line = line.strip()
            if not line:
                continue
            if line.startswith(">"):
                save_current()
                # Parse target name from fasta header (e.g. >EGFR_human)
                current_target = line[1:].split()[0].strip()
                seq_lines = []
            else:
                seq_lines.append(line)
        save_current()

    session.commit()
    session.close()
    console.print(
        f"[bold green]✓ Updated reference sequences for {updated} target(s).[/bold green]"
    )


def add_ligand(
    path: str,
    name: str | None = None,
    batch: str | None = None,
    force: bool = False,
) -> None:
    """
    Imports chemical substances (SDF/SMILES) into the Substance table.
    Performs topological deduplication via Standard InChIKey and records provenance batch_id.
    """
    session = init_db()
    target_path = Path(path).resolve()
    if not target_path.exists():
        console.print(f"[bold red]FATAL: File does not exist: {target_path}[/bold red]")
        sys.exit(1)

    effective_batch = batch.strip() if batch and batch.strip() else target_path.stem

    try:
        from rdkit import Chem
        from rdkit.Chem.MolStandardize import rdMolStandardize
    except ImportError:
        console.print(
            "[bold red]FATAL: RDKit is required for ligand import.[/bold red]"
        )
        sys.exit(1)

    suppl = Chem.SDMolSupplier(str(target_path), removeHs=False, sanitize=False)
    raw_mols = []
    for idx, m in enumerate(suppl, start=1):
        if m is None:
            console.print(
                f"[bold red]FATAL: Corrupted molecular record at index {idx} in {target_path}[/bold red]"
            )
            sys.exit(1)
        raw_mols.append(m)

    if not raw_mols:
        console.print(
            f"[bold red]FATAL: No valid molecules could be parsed from {target_path}[/bold red]"
        )
        sys.exit(1)

    imported = 0
    skipped = 0
    cwd = Path.cwd()

    # STRICT STEREOCHEMICAL INVARIANT:
    # Input stereochemistry (tetrahedral R/S centers and E/Z double-bond geometries such as syn/anti oximes)
    # is taken strictly from the input structure and MUST NEVER be stripped, re-embedded, or enumerated
    # during import. Any conformer removal or vacuum UFF re-optimization would corrupt defined E/Z and R/S parities.
    for idx, raw_mol in enumerate(raw_mols, start=1):
        try:
            if len(raw_mols) == 1 and name:
                base_title = name.strip()
            elif raw_mol.HasProp("_Name") and raw_mol.GetProp("_Name").strip():
                base_title = raw_mol.GetProp("_Name").strip()
            elif name:
                base_title = f"{name}_{idx:05d}"
            else:
                base_title = f"{target_path.stem}_{idx:05d}"

            mol = Chem.Mol(raw_mol)
            (
                parent_inchikey,
                parent_smiles,
                _,
                state_hash,
                isomeric_smiles,
                formal_charge,
            ) = get_microstate_identity(mol)

            mol_title = base_title
            # Disambiguate if another distinct substance (different stereoisomer/InChIKey) already uses mol_title
            name_collision = (
                session.query(Substance)
                .filter(
                    Substance.name == mol_title,
                    Substance.inchikey != parent_inchikey,
                )
                .first()
            )
            if name_collision:
                stereo_suffix = _format_stereo_tag(mol, parent_inchikey)
                mol_title = f"{mol_title}_{stereo_suffix}"
                if (
                    session.query(Substance)
                    .filter(
                        Substance.name == mol_title,
                        Substance.inchikey != parent_inchikey,
                    )
                    .first()
                ):
                    mol_title = f"{base_title}_{parent_inchikey.split('-')[1]}"

            mol.SetProp("_Name", mol_title)

            # Canonical physical substance storage: substances/<inchikey>/<res>_base.sdf
            sub_dir = cwd / "substances" / parent_inchikey
            sub_dir.mkdir(parents=True, exist_ok=True)
            clean_res = (
                "".join(c for c in mol_title if c.isalnum())[:4].upper() or "LIG"
            )
            base_sdf = sub_dir / f"{clean_res}_base.sdf"
            if not base_sdf.exists() or force:
                with Chem.SDWriter(str(base_sdf)) as w:
                    w.write(mol)

            ex_sub = (
                session.query(Substance)
                .filter_by(inchikey=parent_inchikey)
                .first()
            )
            if ex_sub:
                sub_entity = ex_sub
                if force:
                    sub_entity.name = mol_title
                    sub_entity.smiles = parent_smiles
                    sub_entity.path_sdf = str(base_sdf.resolve())
                    sub_entity.state = "raw"
                    sub_entity.batch_id = effective_batch
            else:
                sub_entity = Substance(
                    inchikey=parent_inchikey,
                    smiles=parent_smiles,
                    name=mol_title,
                    res_name=clean_res,
                    substance_type="ligand",
                    state="raw",
                    batch_id=effective_batch,
                    path_sdf=str(base_sdf.resolve()),
                )
                session.add(sub_entity)
                session.flush()

            if ex_sub and not force:
                skipped += 1
                continue

            imported += 1
            if imported % 500 == 0:
                session.commit()
        except Exception as e:
            console.print(
                f"[bold red]FATAL: Substance parsing failed at index {idx}: {e}[/bold red]"
            )
            sys.exit(1)

    session.commit()
    session.close()
    console.print(
        f"[bold green]✓ Imported {imported} unique stereospecific substance(s) into batch '{effective_batch}'[/bold green] "
        f"[dim]({skipped} duplicates skipped)[/dim]"
    )


def add_traj(
    system: str,
    tpr: str,
    xtc: str | None = None,
    raw_dir: str | None = None,
    replicate: int = 1,
    name: str | None = None,
    length_ns: float = 0.0,
    clean_xtc: str | None = None,
    ref_pdb: str | None = None,
) -> None:
    """
    Registers an MD simulation replicate linked to a physical System.
    """
    session = init_db()

    sys_obj = session.query(System).filter_by(name=system).first()
    if not sys_obj:
        console.print(
            f"[bold red]FATAL: Target System '{system}' not found in DB![/bold red]"
        )
        sys.exit(1)

    if raw_dir:
        raw_dir_path = Path(raw_dir).resolve()
    elif xtc:
        raw_dir_path = Path(xtc).resolve().parent
    else:
        console.print(
            "[bold red]FATAL: Either raw_dir or xtc must be provided.[/bold red]"
        )
        sys.exit(1)

    tpr_path = Path(tpr).resolve()

    if not tpr_path.exists():
        console.print(f"[bold red]FATAL: TPR file not found: {tpr_path}[/bold red]")
        sys.exit(1)

    rep_idx = replicate or 1
    traj_name = name or f"{sys_obj.name}_rep{rep_idx}"

    clean_xtc_path_resolved = str(Path(clean_xtc).resolve()) if clean_xtc else None
    ref_pdb_path_resolved = str(Path(ref_pdb).resolve()) if ref_pdb else None
    status = "cleaned" if clean_xtc_path_resolved and ref_pdb_path_resolved else "raw"

    ex = session.query(Trajectory).filter_by(name=traj_name).first()
    if ex:
        ex.system_id = sys_obj.id
        ex.replicate = rep_idx
        ex.raw_dir = str(raw_dir_path)
        ex.path_tpr = str(tpr_path)
        ex.path_clean_xtc = clean_xtc_path_resolved
        ex.path_ref_pdb = ref_pdb_path_resolved
        ex.status = status
        console.print(f"[yellow]Updated Trajectory '{traj_name}' ({status})[/yellow]")
    else:
        new_traj = Trajectory(
            system_id=sys_obj.id,
            name=traj_name,
            replicate=rep_idx,
            length_ns=float(length_ns or 0.0),
            raw_dir=str(raw_dir_path),
            path_tpr=str(tpr_path),
            path_clean_xtc=clean_xtc_path_resolved,
            path_ref_pdb=ref_pdb_path_resolved,
            status=status,
        )
        session.add(new_traj)
        console.print(
            f"[bold green]✓ Registered Trajectory '{traj_name}' ({status}) linked to System '{sys_obj.name}'[/bold green]"
        )

    session.commit()
    session.close()
