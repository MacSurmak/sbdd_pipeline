/*
 * workflows/heal.nf
 * Orchestrates AI-driven Loop Healing (ColabFold + OpenMM)
 */

include { HEAL_STRUCTURE } from '../modules/ai/inpaint'

workflow HEAL_WF {
    if (!params.heal_plan) { error "FATAL: Parameter --heal_plan (JSON) is required." }
    
    def planFile = file(params.heal_plan)
    def planData = new groovy.json.JsonSlurperClassic().parseText(planFile.text)

    ch_tasks = Channel.fromList(planData).map { task ->
        def raw_pdb = file(task.path_raw)
        if (!raw_pdb.exists()) error "FATAL: Raw PDB not found: ${task.path_raw}"

        def meta = [
            id: task.frame_id,
            target: task.target_name,
            fasta_seq: task.fasta_seq,
            project: params.project_name
        ]
        [ meta, raw_pdb ]
    }

    HEAL_STRUCTURE(ch_tasks)

    HEAL_STRUCTURE.out.healed_pdb.map { meta, pdb ->
        "${meta.id}\t${params.outdir}/healed/${meta.id}/${pdb.name}"
    }.collectFile(name: 'heal_receipt.tsv', storeDir: params.outdir, newLine: true)
}
