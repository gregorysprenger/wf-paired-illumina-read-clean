#!/usr/bin/env nextflow


/*
==============================================================================
                              wf-paired-illumina-read-clean                              
==============================================================================
usage: nextflow run main.nf [--help]
----------------------------------------------------------------------------
*/

def helpMessage() {
    log.info"""
    =========================================
     wf-paired-illumina-read-clean v${version}
    =========================================
    Usage:
    nextflow run -profile <docker|singularity> main.nf --inpath <input directory> --outpath <directory for results>

    Run with test data:
    nextflow run main.nf -profile test,<docker|singularity> --outpath results
    
    Input/output options:
      --inpath             Path to input data directory containing FastA assemblies. Recognized extensions are:  fa, fasta, fas, fna, fsa, fa.gz, fasta.gz, fas.gz, fna.gz, fsa.gz.
      --outpath            The output directory where the results will be saved.
    Analysis options:
      --curated_input      Whether or not input is a curated genome directory. If true, will assemue all genomes are similar enough to return sensible results. Options are: true (default), false.
      --recombination      Use a program to classify SNPs as due to recombination. Options are: gubbins, cfml, both.
      --tree_method        Program used to infer trees (in ParSNP and optionally again after masking positions due to recombination). Options are: fasttree (default), raxml.
      --max_partition_size Max partition size (in bases, limits ParSNP memory usage). Note: results can change slightly depending on this value. Default is: 15000000.
      --bigdata            Whether or not to use more compute resources. Options are true, false (default).
      --max_memory         Specify memory limit on your machine/infrastructure, e.g. '128.GB'. Useful to ensure workflow doesn't request too many resources.
      --max_time           Specify time limit for each process, e.g. '240.h'. Useful to ensure workflow doesn't request too many resources.
      --max_cpus           Specify CPU limit on your machine/infrastructure, e.g. 16. Useful to ensure workflow doesn't request too many resources.
      --min_ggr_size       Minimum filesize to verify that ParSNP produced a .ggr file. Default is: '100k'.
      --min_xmfa_size      Minimum filesize to verify that ParSNP produced a .xmfa file. Default is: '100k'.
    Profile options:
      -profile singularity Use Singularity images to run the workflow. Will pull and convert Docker images from Dockerhub if not locally available.
      -profile docker      Use Docker images to run the workflow. Will pull images from Dockerhub if not locally available.
      -profile conda       TODO: this is not implemented yet.
    Other options:
      -resume              Re-start a workflow using cached results. May not behave as expected with containerization profiles docker or singularity.
      -stub                Use example output files for any process with an uncommented stub block. For debugging/testing purposes.
      -name                Name for the pipeline run. If not specified, Nextflow will automatically generate a random mnemonic
    """.stripIndent()
}

version = "1.0.0"
nextflow.enable.dsl=2

if (params.help) {
    helpMessage()
    exit 0
}

if (params.version){
    println "VERSION: $version"
    exit 0
}

/*
========================================================================================
    VALIDATE INPUTS
========================================================================================
*/

if (!params.inpath) {
    System.err.println "ERROR: parameter inpath must be specified"
    exit 1
}
File inpathFileObj = new File(params.inpath)
if (!inpathFileObj.exists()){
    System.err.println "ERROR: $params.inpath doesn't exist"
    exit 1
}

File outpathFileObj = new File(params.outpath)
if (outpathFileObj.exists()){
    // Per the config file, outpath stores log & trace files so it is created before this point
    // Check that outpath only contains a trace file created this hour
    dayAndHour = new java.util.Date().format('yyyy-MM-dd HH')
    outFiles = outpathFileObj.list()
    if (!(outFiles[0] ==~ /trace.($dayAndHour):\d\d:\d\d.txt/ && outFiles.size() == 1)) {
        // If it contains an older trace file or other files, warn the user
        System.out.println "WARNING: $params.outpath already exists. Output files will be overwritten."
    }
} else {
    outpathFileObj.mkdirs()
}

File logpathFileObj = new File(params.logpath)
if (logpathFileObj.exists()){
    System.out.println "WARNING: $params.logpath already exists. Log files will be overwritten."
} else {
    logpathFileObj.mkdirs()
}

// Print parameters used
log.info """
    =====================================
    wf-paired-illumina-read-clean $version
    =====================================
    inpath:             ${params.inpath}
    outpath:            ${params.outpath}
    logpath:            ${params.logpath}
    workDir:            ${workflow.workDir}
    =====================================
    """
    .stripIndent()

/*
========================================================================================
                 Import local custom modules and subworkflows                 
========================================================================================
*/

include { INFILE_HANDLING } from "./modules/local/infile_handling.nf"
include { REMOVE_PHIX } from "./modules/local/remove_phix.nf"
include { TRIMMOMATIC } from "./modules/local/trimmomatic.nf"
include { EXTRACT_SINGLETONS } from "./modules/local/extract_singleton.nf"
include { KRAKEN_ONE; KRAKEN_TWO; } from "./modules/local/kraken.nf"
include { SPADES } from "./modules/local/spades.nf"
include { FILTER_CONTIGS } from "./modules/local/filter_contigs.nf"
include { CLEAN_READS } from "./modules/local/clean_reads.nf"
include { CLEANED_COVERAGE } from "./modules/local/cleaned_coverage.nf"
include { MLST } from "./modules/local/mlst.nf"
include { ANNOTATE } from "./modules/local/annotate.nf"
include { EXTRACT_RECORDS } from "./modules/local/extract_records.nf"
include { BARRNAP } from "./modules/local/barrnap.nf"
include { BLAST } from "./modules/local/blast.nf"
include { FILTER_BLAST } from "./modules/local/filter_blast.nf"
include { QA } from "./modules/local/qa.nf"
include { GENOME_COVERAGE } from "./modules/local/genome_coverage.nf"

/*
========================================================================================
                   Import nf-core modules and subworkflows                    
========================================================================================
*/

// None

/*
========================================================================================
                            Run the main workflow                             
========================================================================================
*/

workflow {
    input_ch = Channel.fromPath(params.inpath, checkIfExists: true)
    output_ch = Channel.fromPath(params.outpath)
    phix_ch = Channel.fromPath('files/PhiX_NC_001422.1.fasta', checkIfExists: true)
    adapters_ch = Channel.fromPath('files/adapters_Nextera_NEB_TruSeq_NuGEN_ThruPLEX.fas', checkIfExists: true)
    filter_contigs_ch = Channel.fromPath('files/filter.contigs.py', checkIfExists: true)
    extract_record_ch = Channel.fromPath('files/extract.record.from.genbank.py', checkIfExists: true)
    filter_blast_ch = Channel.fromPath('files/filter.blast.py', checkIfExists: true)

    INFILE_HANDLING (
        input_ch
    )

    REMOVE_PHIX (
        phix_ch,
        INFILE_HANDLING.out.R1,
        INFILE_HANDLING.out.R2
    )

    TRIMMOMATIC (
        adapters_ch,
        REMOVE_PHIX.out.noPhiX_R1,
        REMOVE_PHIX.out.noPhiX_R2,
        output_ch
    )

    EXTRACT_SINGLETONS (
        INFILE_HANDLING.out.R1,
        TRIMMOMATIC.out.R1_paired,
        TRIMMOMATIC.out.R2_paired
    )

    KRAKEN_ONE (
        EXTRACT_SINGLETONS.out.R1_paired_gz,
        EXTRACT_SINGLETONS.out.R2_paired_gz,
        EXTRACT_SINGLETONS.out.single_gz
    )

    KRAKEN_TWO (
        EXTRACT_SINGLETONS.out.R1_paired_gz,
        EXTRACT_SINGLETONS.out.R2_paired_gz,
        EXTRACT_SINGLETONS.out.single_gz
    )

    SPADES (
        EXTRACT_SINGLETONS.out.R1_paired_gz,
        EXTRACT_SINGLETONS.out.R2_paired_gz,
        EXTRACT_SINGLETONS.out.single_gz,
        output_ch
    )

    FILTER_CONTIGS (
        filter_contigs_ch,
        SPADES.out.contigs,
        EXTRACT_SINGLETONS.out.R1_paired_gz,
        output_ch
    )

    CLEAN_READS (
        FILTER_CONTIGS.out.uncorrected_contigs,
        EXTRACT_SINGLETONS.out.R1_paired_gz,
        EXTRACT_SINGLETONS.out.R2_paired_gz,
        EXTRACT_SINGLETONS.out.single_gz,
        output_ch
    )

    CLEANED_COVERAGE (
        CLEAN_READS.out.single_bam,
        CLEAN_READS.out.paired_bam
    )

    MLST (
        CLEAN_READS.out.base_fna
    )

    ANNOTATE (
        CLEAN_READS.out.base_fna
    )

    EXTRACT_RECORDS (
        extract_record_ch,
        ANNOTATE.out.annotation
    )

    BARRNAP (
        EXTRACT_RECORDS.out.extracted_rna,
        CLEAN_READS.out.base_fna,
        ANNOTATE.out.annotation
    )

    BLAST (
        BARRNAP.out.extracted_base,
        CLEAN_READS.out.base_fna
    )

    FILTER_BLAST (
        filter_blast_ch,
        BLAST.out.blast_tsv,
        CLEAN_READS.out.base_fna,
        output_ch
    )

    QA (
        CLEAN_READS.out.base_fna,
        EXTRACT_SINGLETONS.out.R1_paired_gz,
        EXTRACT_SINGLETONS.out.R2_paired_gz,
        EXTRACT_SINGLETONS.out.single_gz
    )

    GENOME_COVERAGE (
        CLEANED_COVERAGE.out.summary_stats,
        QA.out.summary_assemblies,
        QA.out.summary_bases
    )
}

/*
========================================================================================
                        Completion e-mail and summary                         
========================================================================================
*/

workflow.onComplete {
    log.info """
                |=====================================
                |Pipeline Execution Summary
                |=====================================
                |Workflow Version : ${version}
                |Nextflow Version : ${nextflow.version}
                |Command Line     : ${workflow.commandLine}
                |Resumed          : ${workflow.resume}
                |Completed At     : ${workflow.complete}
                |Duration         : ${workflow.duration}
                |Success          : ${workflow.success}
                |Exit Code        : ${workflow.exitStatus}
                |Launch Dir       : ${workflow.launchDir}
                |=====================================
             """.stripMargin()
}

workflow.onError {
    def err_msg = """
                     |=====================================
                     |Error summary
                     |=====================================
                     |Completed at : ${workflow.complete}
                     |exit status  : ${workflow.exitStatus}
                     |workDir      : ${workflow.workDir}
                     |Error Report :
                     |${workflow.errorReport ?: '-'}
                     |=====================================
                  """.stripMargin()
    log.info err_msg
}


