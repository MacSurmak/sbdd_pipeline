#!/usr/bin/env python3
import argparse
import sys
import subprocess
import shutil
import re
from pathlib import Path

def run_gmx(cmd, input_str=None):
    print(f"Exec: {' '.join(cmd)}")
    res = subprocess.run(cmd, input=input_str, text=True, check=False)
    if res.returncode != 0:
        print(f"FATAL: GROMACS command failed: {' '.join(cmd)}", file=sys.stderr)
        sys.exit(res.returncode)

def cmd_cluster(args):
    tgt = args.xtcs[0]
    if len(args.xtcs) > 1:
        run_gmx(["gmx", "trjcat", "-f"] + args.xtcs + ["-o", "combined_clean.xtc", "-cat"])
        tgt = "combined_clean.xtc"
        
    run_gmx(["gmx", "make_ndx", "-f", args.ref_pdb, "-o", "index.ndx"], input_str="q\n")
    subprocess.run([sys.executable, f"{Path(__file__).parent}/gromacs_ndx_fixer.py", "--ndx", "index.ndx", "--action", "ensure_complex"], check=True)
    
    grp = args.group
    if grp.lower() == "pocket":
        subprocess.run([sys.executable, f"{Path(__file__).parent}/gromacs_ndx_fixer.py", "--ndx", "index.ndx", "--pdb", args.ref_pdb, "--ligand", args.ligand, "--action", "add_pocket"], check=True)
        grp = "Pocket"
    elif grp.lower() == "ligand":
        grp = args.ligand
        
    cmd = ["gmx", "cluster", "-f", tgt, "-s", args.ref_pdb, "-n", "index.ndx", "-method", args.method, "-cutoff", str(args.cutoff)]
    if args.dt and args.dt != "0":
        cmd.extend(["-dt", str(args.dt)])
    if args.max_clusters and int(args.max_clusters) > 0:
        cmd.extend(["-wcl", str(args.max_clusters)])
    cmd.extend(["-cl", "clusters_multi.pdb"])
    
    run_gmx(cmd, input_str=f"{grp}\nCOMPLEX\n")
    subprocess.run(["schrodinger", "run", "python3", f"{Path(__file__).parent}/split_cluster_pdb.py", "--input", "clusters_multi.pdb", "--sys_name", args.sys_name, "--ligand", args.ligand, "--max_clusters", str(args.max_clusters), "--out_receipt", "cluster_receipt.tsv"], check=True)

def cmd_trjconv(args):
    tgt = args.xtcs[0]
    if len(args.xtcs) > 1:
        run_gmx(["gmx", "trjcat", "-f"] + args.xtcs + ["-o", "scratch_concat.xtc", "-cat"])
        tgt = "scratch_concat.xtc"
        
    run_gmx(["gmx", "make_ndx", "-f", args.tpr, "-o", "index.ndx"], input_str="q\n")
    subprocess.run([sys.executable, f"{Path(__file__).parent}/gromacs_ndx_fixer.py", "--ndx", "index.ndx", "--action", "ensure_complex"], check=True)
    run_gmx(["gmx", "trjconv", "-s", args.tpr, "-f", tgt, "-o", "pbc_whole.xtc", "-pbc", "whole"], input_str="0\n")
    run_gmx(["gmx", "trjconv", "-s", args.tpr, "-f", "pbc_whole.xtc", "-n", "index.ndx", "-o", "pbc_centered.xtc", "-center", "-pbc", "mol", "-ur", "compact"], input_str=f"{args.pbc_center}\n{args.pbc_extract}\n")
    run_gmx(["gmx", "trjconv", "-s", args.tpr, "-f", "pbc_centered.xtc", "-n", "index.ndx", "-o", args.out_xtc, "-fit", "rot+trans"], input_str=f"Protein\n{args.pbc_extract}\n")
    run_gmx(["gmx", "trjconv", "-s", args.tpr, "-f", args.out_xtc, "-n", "index.ndx", "-o", args.out_gro, "-dump", "0"], input_str=f"{args.pbc_extract}\n")
    
    for tmp in ["pbc_whole.xtc", "pbc_centered.xtc", "scratch_concat.xtc"]:
        if Path(tmp).exists(): Path(tmp).unlink()

def cmd_equilibration(args):
    if not Path("protocol.txt").exists():
        print("FATAL: protocol.txt playlist not found")
        sys.exit(1)
        
    prev_gro = args.input_gro
    prev_cpt = ""
    with open("protocol.txt", "r") as f:
        for line in f:
            mdp = line.strip()
            if not mdp: continue
            stage = Path(mdp).stem
            cmd = ["gmx", "grompp", "-f", mdp, "-o", f"{stage}.tpr", "-c", prev_gro, "-r", prev_gro, "-p", args.topol, "-n", args.ndx, "-maxwarn", "2"]
            if prev_cpt and "em" not in stage and "minim" not in stage:
                cmd.extend(["-t", prev_cpt])
            run_gmx(cmd)
            
            run_cmd = ["gmx", "mdrun", "-deffnm", stage, "-nb", "gpu", "-ntmpi", "1", "-ntomp", str(args.cpus)]
            if "em" not in stage and "minim" not in stage:
                run_cmd.extend(["-pme", "gpu"])
            run_gmx(run_cmd)
            
            prev_gro = f"{stage}.gro"
            prev_cpt = f"{stage}.cpt"
            
    shutil.copy2(prev_gro, "eq_final.gro")
    shutil.copy2(prev_cpt, "eq_final.cpt")
    shutil.copy2(f"{stage}.tpr", "eq_final.tpr")

def cmd_run_chunk(args):
    for p in args.posre_files:
        pp = Path(p)
        if pp.is_symlink():
            real = pp.resolve()
            pp.unlink()
            shutil.copy2(real, pp)
            
    steps = int(round((float(args.chunk_ns) * 1000.0) / 0.002))
    ext_ps = int(round(float(args.chunk_ns) * 1000.0))
    nstxout = max(1, min(5000, steps // 10))
    nstenergy = max(1, min(1000, steps // 10))
    
    if int(args.chunk_idx) == 1:
        mdp_path = Path("prod.mdp")
        if not mdp_path.exists():
            mdp_path = Path("md.mdp")
        if not mdp_path.exists():
            print("FATAL: Production template (prod.mdp) was not found in staged protocol directory!", file=sys.stderr)
            sys.exit(1)
        mdp = mdp_path.read_text()
        mdp = re.sub(r"nsteps\s*=.*", f"nsteps = {steps}", mdp)
        mdp = re.sub(r"nstxout-compressed\s*=.*", f"nstxout-compressed = {nstxout}", mdp)
        mdp = re.sub(r"nstenergy\s*=.*", f"nstenergy = {nstenergy}", mdp)
        Path("md_run.mdp").write_text(mdp)
        
        run_gmx(["gmx", "grompp", "-f", "md_run.mdp", "-o", "out.tpr", "-c", args.prev_gro, "-r", args.prev_gro, "-t", args.prev_cpt, "-p", args.topol, "-n", args.ndx, "-maxwarn", "3"])
        run_gmx(["gmx", "mdrun", "-deffnm", f"md_{args.chunk_idx}", "-s", "out.tpr", "-nb", "gpu", "-pme", "gpu", "-ntmpi", "1", "-ntomp", str(args.cpus)])
        shutil.copy2("out.tpr", "sys_ref.tpr")
    else:
        run_gmx(["gmx", "convert-tpr", "-s", args.prev_tpr, "-extend", str(ext_ps), "-o", "out.tpr"])
        run_gmx(["gmx", "mdrun", "-deffnm", f"md_{args.chunk_idx}", "-s", "out.tpr", "-cpi", args.prev_cpt, "-noappend", "-nb", "gpu", "-pme", "gpu", "-ntmpi", "1", "-ntomp", str(args.cpus)])
        shutil.copy2("out.tpr", "sys_ref.tpr")
        
    for ext in ["gro", "cpt", "xtc"]:
        for f in Path(".").glob(f"md_{args.chunk_idx}.part*.{ext}"):
            f.rename(f"md_{args.chunk_idx}.{ext}")
            
    shutil.copy2(f"md_{args.chunk_idx}.gro", "out.gro")
    shutil.copy2(f"md_{args.chunk_idx}.cpt", "out.cpt")

def cmd_mmpbsa(args):
    import MDAnalysis as mda
    u = mda.Universe(args.ref_pdb, args.clean_xtc)
    n_frames = len(u.trajectory)
    dt_ns = (u.trajectory[1].time - u.trajectory[0].time) / 1000.0 if n_frames > 1 else 0.01
    start_frame = min(n_frames, max(1, int(round(float(args.start_ns) / dt_ns)) + 1))
    
    pb_str = f"&pb\nmemopt=1,\nmthick={args.membrane_thickness},\nmctrdz=0.0,\nemem=7.0,\nindi=4.0,\nistrng={args.salt_conc},\nradiopt=0,\nporetype=1,\nfillratio=1.25,\nipb=1,\ninp=2,\nsasopt=0,\nsolvopt=2,\nbcopt=10,\nnfocus=1,\nlinit=1000,\neneopt=1,\ncutfd=7.0,\ncutnb=99.0,\nmaxarcdot=15000\n/\n" if args.membrane else f"&gb\nigb=2,\nsaltcon={args.salt_conc}\n/\n"
    
    with open("mmpbsa_pb.in", "w") as f:
        f.write(f"&general\nsys_name=\"Protein_Ligand\",\nstartframe={start_frame},\ninterval={args.step_frames},\ntemperature={args.temp}\n/\n{pb_str}")
        
    mpirun = ["mpirun", "--oversubscribe", "-np", str(args.cpus), "gmx_MMPBSA", "MPI"] if int(args.cpus) > 1 else ["gmx_MMPBSA"]
    run_gmx(mpirun + ["-O", "-i", "mmpbsa_pb.in", "-cs", args.ref_pdb, "-ci", args.ndx, "-cg", "Protein", args.ligand, "-ct", args.clean_xtc, "-cp", args.topol, "-nogui", "-o", "thermodynamics_mmpbsa_pb.dat", "-eo", "thermodynamics_mmpbsa_pb.csv"])
    
    if args.decomp:
        with open("mmpbsa_gb.in", "w") as f:
            f.write(f"&general\nsys_name=\"Protein_Ligand\",\nstartframe={start_frame},\ninterval={args.step_frames},\ntemperature={args.temp}\n/\n&gb\nigb=2,\nsaltcon={args.salt_conc}\n/\n&decomp\nidecomp=2,\ndec_verbose=0\n/\n")
        run_gmx(mpirun + ["-O", "-prefix", "_GMXMMPBSA_GB_", "-i", "mmpbsa_gb.in", "-cs", args.ref_pdb, "-ci", args.ndx, "-cg", "Protein", args.ligand, "-ct", args.clean_xtc, "-cp", args.topol, "-nogui", "-o", "thermodynamics_mmpbsa_gb.dat", "-deo", "decomposition_mmpbsa_gb.csv"])

def cmd_trjcat(args):
    if len(args.xtcs) > 1:
        run_gmx(["gmx", "trjcat", "-f"] + args.xtcs + ["-cat", "-o", args.out_xtc])
    else:
        shutil.copy2(args.xtcs[0], args.out_xtc)

def cmd_ndx_posre(args):
    run_gmx(["gmx", "make_ndx", "-f", args.input_gro, "-o", "index.ndx"], input_str="q\n")
    subprocess.run([sys.executable, f"{Path(__file__).parent}/gromacs_ndx_fixer.py", "--ndx", "index.ndx", "--action", "add_groups"], check=True)
    subprocess.run([sys.executable, f"{Path(__file__).parent}/gromacs_ndx_fixer.py", "--ndx", "index.ndx", "--top", args.topol, "--forces", args.forces, "--action", "generate_posre"], check=True)

def main():
    parser = argparse.ArgumentParser()
    subparsers = parser.add_subparsers(dest="subcommand", required=True)

    p_cluster = subparsers.add_parser("cluster")
    p_cluster.add_argument("--xtcs", nargs="+", required=True)
    p_cluster.add_argument("--ref_pdb", required=True)
    p_cluster.add_argument("--cutoff", required=True)
    p_cluster.add_argument("--method", required=True)
    p_cluster.add_argument("--group", required=True)
    p_cluster.add_argument("--ligand", required=True)
    p_cluster.add_argument("--max_clusters", required=True)
    p_cluster.add_argument("--sys_name", required=True)
    p_cluster.add_argument("--dt")

    p_trjconv = subparsers.add_parser("trjconv")
    p_trjconv.add_argument("--xtcs", nargs="+", required=True)
    p_trjconv.add_argument("--tpr", required=True)
    p_trjconv.add_argument("--pbc_center", required=True)
    p_trjconv.add_argument("--pbc_extract", required=True)
    p_trjconv.add_argument("--out_xtc", required=True)
    p_trjconv.add_argument("--out_gro", required=True)

    p_eq = subparsers.add_parser("equilibration")
    p_eq.add_argument("--input_gro", required=True)
    p_eq.add_argument("--topol", required=True)
    p_eq.add_argument("--ndx", required=True)
    p_eq.add_argument("--cpus", required=True)
    p_eq.add_argument("--posre_files", nargs="+", required=True)

    p_chunk = subparsers.add_parser("run_chunk")
    p_chunk.add_argument("--prev_gro", required=True)
    p_chunk.add_argument("--prev_cpt", required=True)
    p_chunk.add_argument("--prev_tpr", required=True)
    p_chunk.add_argument("--topol", required=True)
    p_chunk.add_argument("--ndx", required=True)
    p_chunk.add_argument("--chunk_ns", required=True)
    p_chunk.add_argument("--chunk_idx", required=True)
    p_chunk.add_argument("--cpus", required=True)
    p_chunk.add_argument("--posre_files", nargs="+", required=True)

    p_mmpbsa = subparsers.add_parser("mmpbsa")
    p_mmpbsa.add_argument("--clean_xtc", required=True)
    p_mmpbsa.add_argument("--ref_pdb", required=True)
    p_mmpbsa.add_argument("--topol", required=True)
    p_mmpbsa.add_argument("--ndx", required=True)
    p_mmpbsa.add_argument("--ligand", required=True)
    p_mmpbsa.add_argument("--start_ns", required=True)
    p_mmpbsa.add_argument("--step_frames", required=True)
    p_mmpbsa.add_argument("--temp", required=True)
    p_mmpbsa.add_argument("--salt_conc", required=True)
    p_mmpbsa.add_argument("--cpus", required=True)
    p_mmpbsa.add_argument("--membrane", action="store_true")
    p_mmpbsa.add_argument("--membrane_thickness", default="28.0")
    p_mmpbsa.add_argument("--decomp", action="store_true")

    p_cat = subparsers.add_parser("trjcat")
    p_cat.add_argument("--xtcs", nargs="+", required=True)
    p_cat.add_argument("--out_xtc", required=True)

    p_ndx = subparsers.add_parser("ndx_posre")
    p_ndx.add_argument("--input_gro", required=True)
    p_ndx.add_argument("--topol", required=True)
    p_ndx.add_argument("--forces", required=True)

    args = parser.parse_args()
    if args.subcommand == "cluster": cmd_cluster(args)
    elif args.subcommand == "trjconv": cmd_trjconv(args)
    elif args.subcommand == "equilibration": cmd_equilibration(args)
    elif args.subcommand == "run_chunk": cmd_run_chunk(args)
    elif args.subcommand == "mmpbsa": cmd_mmpbsa(args)
    elif args.subcommand == "trjcat": cmd_trjcat(args)
    elif args.subcommand == "ndx_posre": cmd_ndx_posre(args)

if __name__ == "__main__":
    main()
