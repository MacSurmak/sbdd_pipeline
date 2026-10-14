#!/usr/bin/env python3
import argparse
import sys
import subprocess
import os
import re
from pathlib import Path


def get_tool(tool_name: str) -> list:
    if "SCHRODINGER" not in os.environ:
        print(
            "FATAL: Environment variable $SCHRODINGER is not defined! "
            "Process must be executed within the isolated Schrödinger Apptainer/host runtime.",
            file=sys.stderr,
        )
        sys.exit(1)

    schrod_dir = Path(os.environ["SCHRODINGER"]).resolve()
    if not schrod_dir.exists():
        print(
            f"FATAL: Directory $SCHRODINGER={schrod_dir} does not exist on this node.",
            file=sys.stderr,
        )
        sys.exit(1)

    bin_path = schrod_dir / tool_name
    if bin_path.exists() and os.access(bin_path, os.X_OK):
        return [str(bin_path)]

    util_path = schrod_dir / "utilities" / tool_name
    if util_path.exists() and os.access(util_path, os.X_OK):
        return [str(util_path)]

    print(
        f"FATAL: Schrödinger binary '{tool_name}' was not found at canonical paths:\n"
        f"  - {bin_path}\n"
        f"  - {util_path}\n"
        f"Execution aborted by Fail-Fast law.",
        file=sys.stderr,
    )
    sys.exit(1)


def distill_schrodinger_errors(job_name: str, raw_stderr: str = "") -> str:
    critical_markers = [
        r"(lib.*\.so\..*?: cannot open shared object file: No such file or directory)",
        r"(Segmentation fault.*)",
        r"(License error.*)",
        r"(FlexLM Error.*)",
        r"(FATAL:.*)",
        r"(Error:.*)",
    ]
    logs_to_check = [raw_stderr]
    if job_name:
        prime_log = Path(f"{job_name}.Prime.log")
        gen_log = Path(f"{job_name}.log")
        if prime_log.exists():
            logs_to_check.append(prime_log.read_text(errors="ignore"))
        if gen_log.exists():
            logs_to_check.append(gen_log.read_text(errors="ignore"))

    for log_text in logs_to_check:
        if not log_text:
            continue
        for marker in critical_markers:
            matches = re.findall(marker, log_text, re.IGNORECASE)
            if matches:
                return matches[0].strip()
    return "Unknown internal Schrödinger error. See raw diagnostic logs."


def run_cmd(cmd, env=None, job_name=None):
    print(f"Exec: {' '.join(cmd)}")
    # Жесткий перехват потоков. Мусор больше не льется напрямую в Nextflow.
    res = subprocess.run(cmd, capture_output=True, text=True, env=env)
    if res.returncode != 0:
        distilled = distill_schrodinger_errors(job_name, res.stderr)
        print(
            f"FATAL: Schrödinger process failed! Exit Code {res.returncode}\nRoot Cause: {distilled}",
            file=sys.stderr,
        )
        sys.exit(res.returncode)
    return res


def cmd_gridgen(args):
    from schrodinger import structure

    st = structure.StructureReader.read(args.ref_lig_mae)
    if st.atom_total == 0:
        print("FATAL: Empty reference ligand", file=sys.stderr)
        sys.exit(1)
    cx = sum(a.x for a in st.atom) / st.atom_total
    cy = sum(a.y for a in st.atom) / st.atom_total
    cz = sum(a.z for a in st.atom) / st.atom_total

    with open("gridgen.in", "w") as f:
        f.write(f"GRIDFILE {args.out_zip}\n")
        f.write(f"RECEP_FILE {args.receptor_mae}\n")
        f.write(f"GRID_CENTER {cx:.4f}, {cy:.4f}, {cz:.4f}\n")
        f.write(f"REF_LIGAND_FILE {args.ref_lig_mae}\n")
        f.write(f"INNERBOX {args.inner_box}\n")
        f.write(f"OUTERBOX {args.outer_box}\n")
        f.write(f"FORCEFIELD {args.forcefield}\n")
        f.write(f"RECEP_VSCALE {args.recep_vscale}\n")
        f.write(f"RECEP_CCUT {args.recep_ccut}\n")

    run_cmd(
        get_tool("glide") + ["gridgen.in", "-OVERWRITE", "-WAIT", "-HOST", "localhost"]
    )


def cmd_glide_dock(args):
    tmp_dir = os.environ.get("SLURM_TMPDIR", "/tmp")
    os.environ["SCHRODINGER_JOBDB"] = f"{tmp_dir}/schrodinger_jobdb_{os.getpid()}"
    Path(os.environ["SCHRODINGER_JOBDB"]).mkdir(parents=True, exist_ok=True)
    with open("dock.in", "w") as f:
        f.write(f"GRIDFILE {args.grid_zip}\n")
        f.write(f"LIGANDFILE {args.ligands_mae}\n")
        f.write(f"PRECISION {args.precision}\n")
        f.write(f"POSES_PER_LIG {args.poses_per_lig}\n")
        f.write("POSTDOCK yes\nEPIK_PENALTIES yes\nSAMPLE_RINGS yes\nWRITE_CSV yes\n")
        f.write(f"FORCEFIELD {args.forcefield}\n")
        f.write(f"LIG_VSCALE {args.lig_vscale}\nLIG_CCUT {args.lig_ccut}\n")
    run_cmd(
        get_tool("glide") + ["dock.in", "-OVERWRITE", "-WAIT", "-HOST", "localhost"]
    )

    if Path("dock_pv.maegz").exists():
        Path("dock_pv.maegz").rename(args.out_pv)
    if Path("dock.csv").exists():
        Path("dock.csv").rename(args.out_csv)

    run_cmd(
        get_tool("run")
        + [
            "python3",
            f"{Path(__file__).parent}/pv_to_sdf.py",
            "--pv",
            args.out_pv,
            "--sdf",
            args.out_sdf,
            "--max_poses",
            args.poses_per_lig,
        ]
    )


def cmd_prime_mmgbsa(args):
    env = os.environ.copy()
    tmp_dir = os.environ.get("SLURM_TMPDIR", "/tmp")
    env["SCHRODINGER_JOBDB"] = f"{tmp_dir}/schrodinger_jobdb_{os.getpid()}"
    Path(env["SCHRODINGER_JOBDB"]).mkdir(parents=True, exist_ok=True)

    cmd = get_tool("prime_mmgbsa") + [
        "-job_type",
        args.job_type,
        "-csv_output",
        "yes",
        "-out_type",
        "PV",
        "-jobname",
        f"{args.job_id}_mmgbsa",
        "-HOST",
        "localhost",
        "-NJOBS",
        "1",
        "-OVERWRITE",
        "-WAIT",
    ]
    if args.membrane:
        cmd.append("-membrane")
    cmd.append(args.pv_maegz)

    job_name = f"{args.job_id}_mmgbsa"
    run_cmd(cmd, env=env, job_name=job_name)

    canonical_csv = Path(f"{args.job_id}_mmgbsa-out.csv")
    canonical_mae = Path(f"{args.job_id}_mmgbsa-out.maegz")

    # Deterministic mapping for Prime output naming conventions
    gen_csv_1 = Path(f"{job_name}-out.csv")
    gen_csv_2 = Path(f"{job_name}.csv")
    if not canonical_csv.exists():
        if gen_csv_1.exists():
            gen_csv_1.rename(canonical_csv)
        elif gen_csv_2.exists():
            gen_csv_2.rename(canonical_csv)

    gen_mae_1 = Path(f"{job_name}-out.maegz")
    gen_mae_2 = Path(f"{job_name}.maegz")
    if not canonical_mae.exists():
        if gen_mae_1.exists():
            gen_mae_1.rename(canonical_mae)
        elif gen_mae_2.exists():
            gen_mae_2.rename(canonical_mae)

    # Fail-Fast inspection: dump Prime diagnostic log on missing artifacts
    if not canonical_csv.exists() or not canonical_mae.exists():
        distilled = distill_schrodinger_errors(job_name, "")
        print(
            f"FATAL: Prime MM-GBSA failed to generate deliverables for task '{args.job_id}' despite exit code 0!\n"
            f"Expected: {canonical_csv.name} and {canonical_mae.name}\n"
            f"Root Cause: {distilled}",
            file=sys.stderr,
        )
        sys.exit(1)


def cmd_light_prep(args):
    in_path = Path(args.input_pdb)
    raw_mae = "raw_cluster.mae"

    if in_path.suffix.lower() in [".mae", ".maegz"]:
        import shutil

        shutil.copy2(in_path, raw_mae)
    else:
        run_cmd(get_tool("structconvert") + [str(in_path), raw_mae])

    from schrodinger import structure

    st = structure.StructureReader.read(raw_mae)

    if args.native_ligand:
        lig_res = args.native_ligand.split("_")[-1].strip().upper()
        lig_atoms = [a.index for a in st.atom if a.pdbres.strip().upper() == lig_res]
        if lig_atoms:
            st.deleteAtoms(lig_atoms)

    # Pure text normalization of AMBER residue names: preserves 100% of MD coordinates, hydrogens, and charges
    res_translation = {
        "HIE": "HIS",
        "HID": "HIS",
        "HIP": "HIS",
        "CYX": "CYS",
        "CYM": "CYS",
        "ASH": "ASP",
        "GLH": "GLU",
        "LYN": "LYS",
    }
    for res in st.residue:
        r_name = res.pdbres.strip().upper()
        if r_name in res_translation:
            res.pdbres = res_translation[r_name]

    st.write(args.out_mae)
    run_cmd(get_tool("structconvert") + [args.out_mae, args.out_pdb])


def cmd_ligprep(args):
    in_flag = (
        "-imae"
        if args.ligand_file.endswith(".mae") or args.ligand_file.endswith(".maegz")
        else "-isd"
    )
    run_cmd(
        get_tool("ligprep")
        + [
            in_flag,
            args.ligand_file,
            "-omae",
            args.out_maegz,
            "-epik",
            "-ph",
            args.target_ph,
            "-pht",
            args.ph_tolerance,
            "-s",
            args.max_states,
            "-bff",
            args.ff_version,
            "-WAIT",
        ]
    )
    if not Path(args.out_maegz).exists():
        print(f"FATAL: Ligprep output missing: {args.out_maegz}", file=sys.stderr)
        sys.exit(1)
    run_cmd(get_tool("structconvert") + [args.out_maegz, args.out_sdf])


def cmd_convert_ligand(args):
    from schrodinger import structure

    st = structure.StructureReader.read(args.in_sdf)
    st.title = args.title.strip()
    st.property["s_m_title"] = args.title.strip()
    st.write(args.out_maegz)
    st.write(args.out_sdf)


def cmd_prepwizard(args):
    complex_path = Path(args.complex_mae).resolve()
    out_path = Path(args.out_mae).resolve()
    target_ph = float(args.target_ph)
    ph_tolerance = float(args.ph_tolerance)
    min_rmsd = float(args.min_rmsd)
    cap_termini = args.cap_termini

    from schrodinger import structure
    from schrodinger.structutils import analyze

    st = structure.StructureReader.read(str(complex_path))

    # Strict chemical fractionation: partition complex into scaffold and organic hetero molecules
    # ASL: "not (protein or water or ions)" identifies organic ligands and organic cofactors
    hetero_indices = analyze.evaluate_asl(st, "not (protein or water or ions)")

    working_complex = complex_path

    if hetero_indices:
        print(
            f"Chemical Fractionation: Detected {len(hetero_indices)} hetero atom(s) (ligand / organic cofactors)."
        )
        hetero_set = set(hetero_indices)
        st_hetero = st.extract(hetero_indices, copy_props=True)
        st_scaffold = st.extract(
            [a.index for a in st.atom if a.index not in hetero_set], copy_props=True
        )

        raw_hetero_file = "raw_hetero.mae"
        # Write each hetero molecule as an independent structure so Epik processes them cleanly
        with structure.StructureWriter(raw_hetero_file) as writer:
            for mol in st_hetero.molecule:
                writer.append(mol.extractStructure(copy_props=True))

        # Run isolated Epik in aqueous bulk solvent mode (guarantees aqueous pKa without rigid pocket penalty)
        # Coordinates of heavy atoms are strictly preserved by Epik
        epik_out_file = "epik_hetero.mae"
        epik_cmd = get_tool("epik") + [
            "-imae",
            raw_hetero_file,
            "-omae",
            epik_out_file,
            "-ph",
            str(target_ph),
            "-pht",
            str(ph_tolerance),
            "-ms",
            "1",
            "-WAIT",
        ]
        run_cmd(epik_cmd, job_name="epik_hetero")

        if not Path(epik_out_file).exists() or Path(epik_out_file).stat().st_size == 0:
            print(
                f"FATAL: Epik failed to generate aqueous protonated states: {epik_out_file}",
                file=sys.stderr,
            )
            sys.exit(1)

        # Re-integrate aqueous-protonated ligand and cofactors into the scaffold (protein + ions)
        for h_st in structure.StructureReader(epik_out_file):
            st_scaffold.extend(h_st)

        fractionated_complex = "complex_epik_aqueous.mae"
        st_scaffold.write(fractionated_complex)
        working_complex = Path(fractionated_complex).resolve()
        print(
            "✓ Successfully integrated aqueous-protonated molecules into complex with native coordinates preserved."
        )

    # Execute PrepWizard:
    # If hetero atoms were pre-protonated via aqueous Epik, we strictly forbid prepwizard
    # from altering ligand/cofactor charges (-noepik).
    # PROPKA will optimize the protein H-bond network around these correctly charged groups.
    prep_cmd = get_tool("prepwizard") + [
        "-fix",
        "-fillsidechains",
        "-propka_pH",
        str(target_ph),
        "-rmsd",
        str(min_rmsd),
        "-preserve_st_titles",
        "-WAIT",
    ]
    if hetero_indices:
        prep_cmd.append("-noepik")
    else:
        prep_cmd.extend(["-epik_pH", str(target_ph), "-epik_pHt", str(ph_tolerance)])

    if cap_termini:
        prep_cmd.append("-captermini")

    prep_cmd.extend([str(working_complex), str(out_path)])
    run_cmd(prep_cmd, job_name="prepwizard")

    if not out_path.exists() or out_path.stat().st_size == 0:
        print(
            f"FATAL: PrepWizard output deliverable missing: {out_path}", file=sys.stderr
        )
        sys.exit(1)


def main():
    parser = argparse.ArgumentParser()
    subparsers = parser.add_subparsers(dest="subcommand", required=True)

    p_grid = subparsers.add_parser("gridgen")
    p_grid.add_argument("--receptor_mae", required=True)
    p_grid.add_argument("--ref_lig_mae", required=True)
    p_grid.add_argument("--inner_box", required=True)
    p_grid.add_argument("--outer_box", required=True)
    p_grid.add_argument("--forcefield", required=True)
    p_grid.add_argument("--recep_vscale", required=True)
    p_grid.add_argument("--recep_ccut", required=True)
    p_grid.add_argument("--out_zip", required=True)

    p_dock = subparsers.add_parser("glide_dock")
    p_dock.add_argument("--grid_zip", required=True)
    p_dock.add_argument("--ligands_mae", required=True)
    p_dock.add_argument("--precision", required=True)
    p_dock.add_argument("--poses_per_lig", required=True)
    p_dock.add_argument("--forcefield", required=True)
    p_dock.add_argument("--lig_vscale", required=True)
    p_dock.add_argument("--lig_ccut", required=True)
    p_dock.add_argument("--out_pv", required=True)
    p_dock.add_argument("--out_csv", required=True)
    p_dock.add_argument("--out_sdf", required=True)

    p_prime = subparsers.add_parser("prime_mmgbsa")
    p_prime.add_argument("--pv_maegz", required=True)
    p_prime.add_argument("--job_type", required=True)
    p_prime.add_argument("--job_id", required=True)
    p_prime.add_argument("--membrane", action="store_true")

    p_light = subparsers.add_parser("light_prep")
    p_light.add_argument("--input_pdb", required=True)
    p_light.add_argument("--out_mae", required=True)
    p_light.add_argument("--out_pdb", required=True)
    p_light.add_argument("--native_ligand")

    p_ligprep = subparsers.add_parser("ligprep")
    p_ligprep.add_argument("--ligand_file", required=True)
    p_ligprep.add_argument("--target_ph", required=True)
    p_ligprep.add_argument("--ph_tolerance", required=True)
    p_ligprep.add_argument("--ff_version", required=True)
    p_ligprep.add_argument("--max_states", required=True)
    p_ligprep.add_argument("--out_maegz", required=True)
    p_ligprep.add_argument("--out_sdf", required=True)

    p_conv = subparsers.add_parser("convert_ligand")
    p_conv.add_argument("--in_sdf", required=True)
    p_conv.add_argument("--out_maegz", required=True)
    p_conv.add_argument("--out_sdf", required=True)
    p_conv.add_argument("--title", required=True)

    p_prep = subparsers.add_parser("prepwizard")
    p_prep.add_argument("--complex_mae", required=True)
    p_prep.add_argument("--out_mae", required=True)
    p_prep.add_argument("--target_ph", required=True)
    p_prep.add_argument("--ph_tolerance", required=True)
    p_prep.add_argument("--min_rmsd", required=True)
    p_prep.add_argument("--cap_termini", action="store_true")

    args = parser.parse_args()
    if args.subcommand == "gridgen":
        cmd_gridgen(args)
    elif args.subcommand == "glide_dock":
        cmd_glide_dock(args)
    elif args.subcommand == "prime_mmgbsa":
        cmd_prime_mmgbsa(args)
    elif args.subcommand == "light_prep":
        cmd_light_prep(args)
    elif args.subcommand == "ligprep":
        cmd_ligprep(args)
    elif args.subcommand == "convert_ligand":
        cmd_convert_ligand(args)
    elif args.subcommand == "prepwizard":
        cmd_prepwizard(args)


if __name__ == "__main__":
    main()
