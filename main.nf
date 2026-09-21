#!/usr/bin/env nextflow
/*
 * main.nf
 * SBDD Pipeline Nextflow HPC Orchestrator Entrypoint
 */

nextflow.enable.dsl = 2

// Include top-level workflows
include { PREPARE_WF  } from './workflows/prepare'
include { SCREEN_WF   } from './workflows/screen'
include { DYNAMICS_WF } from './workflows/dynamics'
include { SBDD_WF     } from './workflows/sbdd'

def helpMessage() {
    log.info """
    ========================================================================
     SBDD Pipeline Nextflow HPC Orchestrator
    ========================================================================
    Usage:
      nextflow run ${projectDir}/main.nf -profile slurm [options]

    Primary Parameters:
      --entry            Target subworkflow: 'prepare', 'screen', 'dynamics', 'all' (default: 'all')
      --project_name     Project name for isolated /shared/work scratch
      --receptor         Path to input receptor structure (.cif or .pdb)
      --library          Path to ligand library (.sdf)
      --native_ligand    3-letter residue name of native reference ligand (default: PKA)
      --dock_engine      Docking engine: 'gnina' or 'glide' (default: gnina)
      --membrane         Enable PPM 3.0 membrane immersion & lipid bilayer (default: true)
      --run_md           Enable downstream MD simulations (default: false)
      --md_hits          Hits to simulate in MD: 'native' or comma-separated list
      --outdir           Target directory for final outputs (default: current directory)
    ========================================================================
    """.stripIndent()
}

workflow {
    if (params.help) {
        helpMessage()
        exit 0
    }

    def entrypoint = params.entry ? params.entry.toLowerCase() : "all"

    if (entrypoint == "prepare") {
        PREPARE_WF()
    } else if (entrypoint == "screen") {
        SCREEN_WF()
    } else if (entrypoint == "dynamics") {
        DYNAMICS_WF()
    } else {
        SBDD_WF()
    }
}

