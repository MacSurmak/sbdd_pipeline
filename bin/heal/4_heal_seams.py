#!/usr/bin/env python3
import argparse
import json
import sys
from openmm.app import PDBFile, ForceField, HBonds, NoCutoff, Simulation
from openmm import CustomExternalForce, LangevinMiddleIntegrator
from openmm.unit import kelvin, picoseconds, kilojoules_per_mole, nanometer, nanometers, femtoseconds
from pdbfixer import PDBFixer

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--stitched", required=True)
    parser.add_argument("--map", required=True)
    parser.add_argument("--out", required=True)
    args = parser.parse_args()

    with open(args.map, 'r') as f:
        mapping = json.load(f)

    relax_residues = set()
    for g in mapping["internal_gaps"]:
        start_relax = g["af2_start"] - 2
        end_relax = g["af2_end"] + 2
        for i in range(start_relax, end_relax + 1):
            relax_residues.add(i)

    fixer = PDBFixer(filename=args.stitched)
    fixer.missingResidues = {}
    fixer.findMissingAtoms()
    fixer.addMissingAtoms()
    fixer.addMissingHydrogens(7.4)

    forcefield = ForceField('amber14-all.xml', 'implicit/obc2.xml')
    system = forcefield.createSystem(fixer.topology, nonbondedMethod=NoCutoff, constraints=HBonds)
    
    force = CustomExternalForce("k*((x-x0)^2+(y-y0)^2+(z-z0)^2)")
    force.addGlobalParameter("k", 10000.0 * kilojoules_per_mole/nanometer**2)
    force.addPerParticleParameter("x0")
    force.addPerParticleParameter("y0")
    force.addPerParticleParameter("z0")

    restrained = 0
    for atom in fixer.topology.atoms():
        try:
            res_id = int(atom.residue.id)
        except ValueError:
            res_id = -1
            
        if not (res_id in relax_residues):
            if atom.element.symbol != 'H':
                pos = fixer.positions[atom.index]
                force.addParticle(atom.index, (pos[0], pos[1], pos[2]))
                restrained += 1

    system.addForce(force)

    integrator = LangevinMiddleIntegrator(300*kelvin, 1/picoseconds, 2*femtoseconds)
    simulation = Simulation(fixer.topology, system, integrator)
    simulation.context.setPositions(fixer.positions)

    simulation.minimizeEnergy(maxIterations=2000)

    state = simulation.context.getState(getPositions=True)
    with open(args.out, 'w') as f:
        PDBFile.writeFile(simulation.topology, state.getPositions(), f, keepIds=True)
    
    print(f"OpenMM Relaxation Complete. Restrained {restrained} heavy atoms.")

if __name__ == "__main__":
    main()
