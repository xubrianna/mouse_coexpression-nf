#!/usr/bin/env nextflow
nextflow.enable.dsl=2

/*
 * Mouse Coexpression Analysis Pipeline
 * Complete pipeline: filtering → gene selection → CPM normalization → coexpression computation
 */

// Print parameter summary
log.info """\
    MOUSE COEXPRESSION PIPELINE
    ===========================
    input_dir      : ${params.input_dir}
    outdir         : ${params.outdir}
    gene_frequency : ${params.gene_frequency_threshold}
    gene_strategy  : ${params.gene_selection_strategy}
    correlation    : ${params.correlation_type}
    """
    .stripIndent()

/*
 * Process 1: Generate gene lists across all cell types and studies
 */
process GENERATE_GENE_LISTS {
    publishDir "${params.outdir}/01_gene_lists", mode: 'copy'
    memory '8 GB'
    
    conda "${params.conda_env}"
    
    input:
    path "filtered/*"
    
    output:
    path "genes-lists-*.json", emit: gene_lists
    
    script:
    """
    python ${projectDir}/bin/generateGeneList.py \\
        --input_dir filtered \\
        --output genes-lists-${params.gene_frequency_threshold}.json \\
        --frequency_threshold ${params.gene_frequency_threshold}
    """
}

/*
 * Process 2: Select genes and apply CPM normalization
 */
process SELECT_GENES_AND_CPM {
    tag "$h5ad.baseName"
    publishDir "${params.outdir}/02_cpm_normalized", mode: 'copy'
    
    conda "${params.conda_env}"
    
    input:
    path h5ad
    path gene_lists
    
    output:
    path "cpm/${h5ad.baseName}.h5ad", emit: cpm_h5ad
    
    script:
    """
    python ${projectDir}/bin/selectGenesCPM.py \\
        --input ${h5ad} \\
        --output cpm/${h5ad.baseName}.h5ad \\
        --gene_lists ${gene_lists} \\
        --strategy ${params.gene_selection_strategy} \\
        --cpm_normalize \\
        --no_log_transform
    """
}

/*
 * Process 3: Compute coexpression for each h5ad file
 */
process COMPUTE_COEXPRESSION {
    tag "$h5ad.baseName"
    publishDir "${params.outdir}/03_coexpression", mode: 'copy'
    
    conda "${params.conda_env}"
    
    input:
    path h5ad
    
    output:
    path "agg-mtx/*.h5ad", optional: true
    path "edge-list/upTriang/*.csv", optional: true
    path "edge-list/all-matrix/*.csv", optional: true
    
    script:
    def replace_nans_arg = params.replace_nans ? '--replace_nans' : ''
    def min_cells_arg = params.min_cells ? "--min_cells_threshold ${params.min_cells}" : ''
    """
    export OMP_NUM_THREADS=${task.cpus}
    export OPENBLAS_NUM_THREADS=${task.cpus}
    export MKL_NUM_THREADS=${task.cpus}
    export NUMEXPR_NUM_THREADS=${task.cpus}

    python ${projectDir}/bin/computeCorrelation.py \\
        --input_file ${h5ad} \\
        --output_dir . \\
        --correlation_type ${params.correlation_type} \\
        ${replace_nans_arg} \\
        ${min_cells_arg}
    """
}

/*
 * Main workflow
 */
workflow {
    // Step 1: Generate gene lists from all filtered files
    h5ad_files_ch = Channel.fromPath("${params.input_dir}/*.h5ad")
    all_filtered = h5ad_files_ch.collect()
    gene_lists = GENERATE_GENE_LISTS(all_filtered)
    
    // Recreate channel for next step (since it was consumed by collect())
    h5ad_files_for_cpm = Channel.fromPath("${params.input_dir}/*.h5ad")
    
    // Step 2: Select genes and apply CPM normalization
    cpm_normalized = SELECT_GENES_AND_CPM(
        h5ad_files_for_cpm,
        gene_lists.gene_lists
    )
    
    // Step 3: Compute coexpression matrices (on CPM normalized output)
    COMPUTE_COEXPRESSION(cpm_normalized.cpm_h5ad)
}

workflow.onComplete {
    log.info """\
        Pipeline completed!
        ==================
        Status    : ${workflow.success ? 'SUCCESS' : 'FAILED'}
        Duration  : ${workflow.duration}
        Results   : ${params.outdir}
        
        Output Structure:
        - 01_gene_lists/     : Gene intersection lists (JSON)
        - 02_cpm_normalized/ : CPM normalized h5ad files
        - 03_coexpression/   : Coexpression matrices and edge lists
        """
        .stripIndent()
}
