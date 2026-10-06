#!/usr/bin/env python3
import argparse
import sys
import subprocess
import os
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

def run_cmd(cmd, env=None):
    print(f"Exec: {' '.join(cmd)}")
    res = subprocess.run(cmd, check=False, env=env)
    if res.returncode != 0:
        print(f"FATAL: Command failed with exit code {res.returncode}", file=sys.stderr)
        sys.exit(res.returncode)

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
        f.write("RECEP_CCUT 0.25\n")
        
    run_cmd(get_tool("glide") + ["gridgen.in", "-OVERWRITE", "-WAIT", "-HOST", "localhost"])

def cmd_glide_dock(args):
    os.environ["SCHRODINGER_JOBDB"] = f"/tmp/schrodinger_jobdb_{os.getpid()}"
    Path(os.environ["SCHRODINGER_JOBDB"]).mkdir(parents=True, exist_ok=True)
    with open("dock.in", "w") as f:
        f.write(f"GRIDFILE {args.grid_zip}\n")
        f.write(f"LIGANDFILE {args.ligands_mae}\n")
        f.write(f"PRECISION {args.precision}\n")
        f.write(f"POSES_PER_LIG {args.poses_per_lig}\n")
        f.write("POSTDOCK yes\nEPIK_PENALTIES yes\nSAMPLE_RINGS yes\nWRITE_CSV yes\n")
        f.write(f"FORCEFIELD {args.forcefield}\n")
        f.write(f"LIG_VSCALE {args.lig_vscale}\nLIG_CCUT 0.15\n")
    run_cmd(get_tool("glide") + ["dock.in", "-OVERWRITE", "-WAIT", "-HOST", "localhost"])
    
    if Path("dock_pv.maegz").exists():
        Path("dock_pv.maegz").rename(args.out_pv)
    if Path("dock.csv").exists():
        Path("dock.csv").rename(args.out_csv)
        
    run_cmd(get_tool("run") + ["python3", f"{Path(__file__).parent}/pv_to_sdf.py", "--pv", args.out_pv, "--sdf", args.out_sdf, "--max_poses", args.poses_per_lig])

def cmd_prime_mmgbsa(args):
    libs = ["libtirpc.so.3", "libkeyutils.so.1", "libgssapi_krb5.so.2", "libkrb5.so.3", "libk5crypto.so.3", "libcom_err.so.2", "libkrb5support.so.0"]
    binds = []
    for lib in libs:
        for loc in ["/lib64", "/usr/lib64", "/lib/x86_64-linux-gnu", "/usr/lib/x86_64-linux-gnu", "/usr/lib"]:
            p = Path(loc) / lib
            if p.exists():
                binds.extend([f"{p}:/usr/lib/x86_64-linux-gnu/{lib}", f"{p}:/lib/x86_64-linux-gnu/{lib}", f"{p}:/usr/lib64/{lib}"])
                break
    env = os.environ.copy()
    if binds:
        new_binds = ",".join(binds)
        for var in ["APPTAINER_BINDPATH", "SINGULARITY_BINDPATH"]:
            ext = env.get(var, "")
            env[var] = f"{new_binds},{ext}" if ext else new_binds
            
    env["SCHRODINGER_JOBDB"] = f"/tmp/schrodinger_jobdb_{os.getpid()}"
    Path(env["SCHRODINGER_JOBDB"]).mkdir(parents=True, exist_ok=True)
    
    cmd = get_tool("prime_mmgbsa") + ["-job_type", args.job_type, "-csv_output", "yes", "-out_type", "PV", "-jobname", f"{args.job_id}_mmgbsa", "-HOST", "localhost", "-NJOBS", "1", "-OVERWRITE", "-WAIT"]
    if args.membrane:
        cmd.append("-membrane")
    cmd.append(args.pv_maegz)
    
    run_cmd(cmd, env=env)
        
    csv_in, csv_out = Path(f"{args.job_id}_mmgbsa.csv"), Path(f"{args.job_id}_mmgbsa-out.csv")
    if csv_in.exists(): csv_in.rename(csv_out)
    
    mae_in, mae_out = Path(f"{args.job_id}_mmgbsa.maegz"), Path(f"{args.job_id}_mmgbsa-out.maegz")
    if mae_in.exists() and not mae_out.exists(): mae_in.rename(mae_out)
        
    if not csv_out.exists() or not mae_out.exists():
        print("FATAL: Prime outputs not found", file=sys.stderr)
        sys.exit(1)

def cmd_light_prep(args):
    run_cmd(get_tool("structconvert") + [args.input_pdb, args.out_mae])
    if args.native_ligand:
        from schrodinger import structure
        st = structure.StructureReader.read(args.out_mae)
        lig_res = args.native_ligand.split('_')[-1].strip().upper()
        lig_atoms = [a.index for a in st.atom if a.pdbres.strip().upper() == lig_res]
        if lig_atoms:
            st.deleteAtoms(lig_atoms)
            st.write(args.out_mae)
    run_cmd(get_tool("structconvert") + [args.out_mae, args.out_pdb])

def cmd_ligprep(args):
    in_flag = "-imae" if args.ligand_file.endswith(".mae") or args.ligand_file.endswith(".maegz") else "-isd"
    run_cmd(get_tool("ligprep") + [in_flag, args.ligand_file, "-omae", args.out_maegz, "-epik", "-ph", args.target_ph, "-pht", args.ph_tolerance, "-s", args.max_states, "-bff", args.ff_version, "-WAIT"])
    if not Path(args.out_maegz).exists():
        print(f"FATAL: Ligprep output missing: {args.out_maegz}", file=sys.stderr)
        sys.exit(1)
    run_cmd(get_tool("structconvert") + [args.out_maegz, args.out_sdf])

def cmd_convert_ligand(args):
    from schrodinger import structure
    st = structure.StructureReader.read(args.in_sdf)
    st.title = args.title.strip()
    st.property['s_m_title'] = args.title.strip()
    st.write(args.out_maegz)
    st.write(args.out_sdf)

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
    p_grid.add_argument("--out_zip", required=True)

    p_dock = subparsers.add_parser("glide_dock")
    p_dock.add_argument("--grid_zip", required=True)
    p_dock.add_argument("--ligands_mae", required=True)
    p_dock.add_argument("--precision", required=True)
    p_dock.add_argument("--poses_per_lig", required=True)
    p_dock.add_argument("--forcefield", required=True)
    p_dock.add_argument("--lig_vscale", required=True)
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

    args = parser.parse_args()
    if args.subcommand == "gridgen": cmd_gridgen(args)
    elif args.subcommand == "glide_dock": cmd_glide_dock(args)
    elif args.subcommand == "prime_mmgbsa": cmd_prime_mmgbsa(args)
    elif args.subcommand == "light_prep": cmd_light_prep(args)
    elif args.subcommand == "ligprep": cmd_ligprep(args)
    elif args.subcommand == "convert_ligand": cmd_convert_ligand(args)

if __name__ == "__main__":
    main()
