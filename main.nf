#!/usr/bin/env nextflow
/*
 * main.nf
 * SBDD Pipeline Nextflow HPC Orchestrator Entrypoint
 */

nextflow.enable.dsl = 2

// Include top-level workflows
include { PREPARE_WF  } from './workflows/prepare'
include { SCREEN_WF   } from './workflows/screen'
include { RESCORE_WF  } from './workflows/rescore'
include { MD_SETUP_WF } from './workflows/md_setup'
include { DYNAMICS_WF } from './workflows/dynamics'
include { ANALYZE_WF  } from './workflows/analyze'
include { HEAL_WF     } from './workflows/heal'

def helpMessage() {
    log.info """
    ========================================================================
     SBDD Pipeline Nextflow HPC Orchestrator
    ========================================================================
    Usage:
      nextflow run ${projectDir}/main.nf -profile slurm [options]

    Primary Parameters:
      --entry            Target subworkflow: 'heal', 'prepare', 'screen', 'rescore', 'dynamics'
      --project_name     Project name for isolated /shared/work scratch
      --receptor         Path to input receptor structure (.cif or .pdb)
      --library          Path to ligand library (.sdf)
      --native_ligand    3-letter residue name of native reference ligand
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
    if (entrypoint == "heal") {
        HEAL_WF()
    } else if (entrypoint == "prepare") {
        PREPARE_WF()
    } else if (entrypoint == "screen" || entrypoint == "dock") {
        SCREEN_WF()
    } else if (entrypoint == "rescore") {
        RESCORE_WF()
    } else if (entrypoint == "md-setup") {
        MD_SETUP_WF()
    } else if (entrypoint == "dynamics") {
        DYNAMICS_WF()
    } else if (entrypoint == "md-analyze") {
        ANALYZE_WF()
    } else {
        if (entrypoint != null && entrypoint != "null") {
            log.error "Unknown entrypoint: ${entrypoint}"
        }
        helpMessage()
        exit 1
    }
}

