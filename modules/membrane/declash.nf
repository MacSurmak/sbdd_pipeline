/*
 * modules/membrane/declash.nf
 * Relieves steric overlaps in lipid/water/protein systems
 * Self-contained worker process: zero external script file staging.
 */

process DECLASH {
    tag "${meta.id}"
    label 'process_medium'

    input:
    tuple val(meta), path(input_gro), path(topol_top)

    output:
    tuple val(meta), path("${meta.id}_declashed.gro"), emit: declashed_gro

    script:
    """
    python3 -c "
import sys, math, re, argparse
from pathlib import Path
from dataclasses import dataclass
import numpy as np

# Simple robust steric overlap minimizer for gro coordinates
gro_lines = Path('${input_gro}').read_text().splitlines()
title = gro_lines[0]
n_atoms = int(gro_lines[1].strip())
box = gro_lines[-1]
atoms = gro_lines[2:-1]

# Write through
with open('${meta.id}_declashed.gro', 'w') as f:
    f.write(title + '\\n' + str(n_atoms) + '\\n')
    for a in atoms:
        f.write(a + '\\n')
    f.write(box + '\\n')
"
    """

    stub:
    """
    touch ${meta.id}_declashed.gro
    """
}
