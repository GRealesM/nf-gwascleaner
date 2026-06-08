#!/usr/bin/env nextflow
nextflow.enable.dsl=2

process INSPECTOR {
    input:
    tuple val(study_id), path(gwas_file)
    path manifest

    output:
    tuple val(study_id), path("meta.json"), path(gwas_file)

    script:
        """
        inspect_gwas.py \
            --input "${gwas_file}" \
            --study_id "${study_id}" \
            --manifest "${manifest}" \
            --output meta.json
        """
}

process STANDARDISE {
    tag "${meta.study_id}"

    input:
    tuple val(meta), path(json_file), path(gwas_file)

    output:
    tuple val(meta), path("${meta.study_id}_std.tsv.gz")

    script:
    """
    standardise_gwas.py \
        --input ${gwas_file} \
        --meta ${json_file} \
        --output ${meta.study_id}_std.tsv.gz
    """
}

process LIFTOVER {
    tag "${meta.study_id}"
    
    // Magic happens here: Nextflow auto-provisions liftOver dynamically
    conda 'ucsc-liftover=482'

    input:
    tuple val(meta), path(gwas_file), path(chain_file)

    output:
    tuple val(meta), path("${meta.study_id}_hg38.tsv.gz"), path("${meta.study_id}_liftover.log")

    script:
    """
    # Note: We now call 'liftOver' directly instead of pointing to a binary file path
    liftover_gwas.py \
        --input ${gwas_file} \
        --chain ${chain_file} \
        --liftover liftOver \
        --output ${meta.study_id}_hg38.tsv.gz \
        --log ${meta.study_id}_liftover.log
    """
}

process SORT_GWAS {
    tag "${meta.study_id}"
    
    // This will be the ONLY folder that saves our final, polished files
    publishDir "${params.outdir}/Final_GWAS", mode: 'copy'

    input:
    tuple val(meta), path(gwas_file)

    output:
    tuple val(meta), path("${meta.study_id}_final.tsv.gz")

    script:
    """
    # 1. Extract the header safely
    zcat ${gwas_file} | head -n 1 > header.tsv
    
    # 2. Extract data, sort it by CHR (Col 1) and BP (Col 2), and save
    # Note: LC_ALL=C disables language locales, making the sort up to 10x faster!
    # -k1,1V sorts Col 1 (CHR) by version (1..22, X, Y)
    # -k2,2n sorts Col 2 (BP) numerically
    zcat ${gwas_file} | tail -n +2 | LC_ALL=C sort -k1,1V -k2,2n > sorted_data.tsv
    
    # 3. Recombine and compress
    cat header.tsv sorted_data.tsv | gzip -c > ${meta.study_id}_final.tsv.gz
    """
}

// Main workflow logic
workflow {
    raw_gwas_ch = Channel
        .fromPath("${params.gwas_dir}/*.gz")
        .map { file -> 
            def auto_id = file.name.replaceAll(/(?i)(\.tsv|\.txt|\.csv)?\.gz$/, "")
            return tuple(auto_id, file) 
        }
        
    manifest_file = file(params.manifest)
    inspected_ch = INSPECTOR(raw_gwas_ch, manifest_file)

    parsed_meta_ch = inspected_ch
        .map { study_id, json_file, gwas_file ->
            def meta = new groovy.json.JsonSlurper().parse(json_file)
            return tuple(meta, json_file, gwas_file)
        }
        .branch { meta, json_file, gwas_file ->
            failed: meta.fatal_error == true
            passed: meta.fatal_error == false
        }

    // --- REPORTING 1: The QC Failures Master File ---
    parsed_meta_ch.failed
        .map { meta, json, file -> 
            "${meta.study_id}\tFAILED\t${meta.error_message}\n" 
        }
        .collectFile(
            name: 'pipeline_qc_failures.tsv', 
            storeDir: "${params.outdir}/reports", 
            seed: "STUDY_ID\tSTATUS\tREASON\n"
        )

    // (We keep the console warning so you still see it live)
    parsed_meta_ch.failed.view { meta, json, file -> 
        "WARNING: ${meta.study_id} failed QC: ${meta.error_message}" 
    }

    // Standardise passing files
    standardised_ch = STANDARDISE(parsed_meta_ch.passed)

    // Branch for Liftover
    build_routing_ch = standardised_ch
        .branch { meta, std_file ->
            needs_lift_18: meta.build == 'hg18'
            needs_lift_19: meta.build == 'hg19'
            ready_hg38: meta.build == 'hg38'
            unknown: meta.build == 'unknown'
        }

    chain_18 = file(params.chain_hg18)
    chain_19 = file(params.chain_hg19)
    
    liftover_inputs_ch = build_routing_ch.needs_lift_18
        .map { meta, std_file -> tuple(meta, std_file, chain_18) }
        .mix(
            build_routing_ch.needs_lift_19
            .map { meta, std_file -> tuple(meta, std_file, chain_19) }
        )

    // Run liftover
    lifted_ch = LIFTOVER(liftover_inputs_ch)

    // --- REPORTING 2: The Liftover Metrics Master File ---
    lifted_ch
        .map { meta, final_file, log_file -> 
            // We read the tiny text file using .text and append it to our row
            "${meta.study_id}\t${log_file.text}\n" 
        }
        .collectFile(
            name: 'liftover_drop_metrics.tsv', 
            storeDir: "${params.outdir}/reports", 
            seed: "STUDY_ID\tLIFTOVER_NOTES\n"
        )

    // To send this to SORT_GWAS, we map the channel to strip away the log_file, 
    // ensuring it matches the 2-item structure of the ready_hg38 channel.
    lifted_for_sort_ch = lifted_ch.map { meta, final_file, log_file -> tuple(meta, final_file) }

    // Mix and Sort
    final_unsorted_ch = lifted_for_sort_ch.mix(build_routing_ch.ready_hg38)
    sorted_ch = SORT_GWAS(final_unsorted_ch)

    // View final completion
    sorted_ch.view { meta, final_file ->
        "SUCCESS: Processing complete -> ${final_file.name}"
    }
    
}
