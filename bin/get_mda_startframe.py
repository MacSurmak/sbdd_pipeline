#!/usr/bin/env python3
"""
bin/get_mda_startframe.py
Calculates GROMACS MM-PBSA startframe based on MDAnalysis trajectory metadata.
"""
import argparse
import MDAnalysis as mda


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--pdb", required=True)
    parser.add_argument("--xtc", required=True)
    parser.add_argument("--start_ns", type=float, required=True)
    args = parser.parse_args()

    u = mda.Universe(args.pdb, args.xtc)
    n_frames = len(u.trajectory)
    dt_ns = (
        (u.trajectory[1].time - u.trajectory[0].time) / 1000.0 if n_frames > 1 else 0.01
    )
    sf = min(n_frames, max(1, int(round(args.start_ns / dt_ns)) + 1))
    print(sf)


if __name__ == "__main__":
    main()
