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
    nextflow run main.nf -profile test,<docker|singularity>
    
    Input/output options:
      --inpath             Path to input data directory containing FastQ assemblies. Recognized extensions are:  fastq.gz, fq.gz.
      --outpath            The output directory where the results will be saved.
    Analysis options:
      --size               Specify file size that is used to verify minimum file sizes in all processes in BYTES, e.g. 1000.
      --kraken1_db         Specify path to database for Kraken1. Default database is Mini Kraken.
      --kraken2_db         Specify path to database for Kraken2. Default database is Mini Kraken.
      --blast_db           Specify path to 16S ribosomal database for BLAST. Default database is NCBI's 16S ribosomal database.
      --bigdata            Whether or not to use more compute resources. Options are true, false (default).
      --max_memory         Specify memory limit on your machine/infrastructure, e.g. '128.GB'. Useful to ensure workflow doesn't request too many resources.
      --max_time           Specify time limit for each process, e.g. '240.h'. Useful to ensure workflow doesn't request too many resources.
      --max_cpus           Specify CPU limit on your machine/infrastructure, e.g. 16. Useful to ensure workflow doesn't request too many resources.
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
    dayAndHour = new java.util.Date().format('yyyy-MM-dd_HH-mm-ss')
    outFiles = outpathFileObj.list()
    if (!(outFiles[0] ==~ /trace.($dayAndHour).txt/ && outFiles.size() == 1)) {
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

    // SETUP: Define input, output, and dependency channels
    input_ch = Channel.fromFilePairs(params.inpath+'/*R{1,2}*.{fastq,fq}.gz', checkIfExists: true)
    ch_versions = Channel.empty()

    // PROCESS: Read files from input directory, validate and stage input files
    INFILE_HANDLING (
        input_ch
    )

    ch_versions = ch_versions.mix(INFILE_HANDLING.out.versions)

    // PROCESS: Run bbduk to remove PhiX reads
    REMOVE_PHIX (
        INFILE_HANDLING.out.input,
        INFILE_HANDLING.out.base,
        INFILE_HANDLING.out.size
    )

    ch_versions = ch_versions.mix(REMOVE_PHIX.out.versions)

    // PROCESS: Run trimmomatic to clip adapters and do quality trimming
    TRIMMOMATIC (
        REMOVE_PHIX.out.noPhiX_R1,
        REMOVE_PHIX.out.noPhiX_R2,
        INFILE_HANDLING.out.base,
        INFILE_HANDLING.out.size
    )

    ch_versions = ch_versions.mix(TRIMMOMATIC.out.versions)

    // PROCESS: Run flash to merge overlapping sister reads into singleton reads
    EXTRACT_SINGLETONS (
        INFILE_HANDLING.out.input,
        TRIMMOMATIC.out.R1_paired,
        TRIMMOMATIC.out.R2_paired,
        INFILE_HANDLING.out.base,
        INFILE_HANDLING.out.size
    )

    ch_versions = ch_versions.mix(EXTRACT_SINGLETONS.out.versions)

    // PROCESS: Run kraken1 on paired reads
    KRAKEN_ONE (
        EXTRACT_SINGLETONS.out.R1_paired_gz,
        EXTRACT_SINGLETONS.out.R2_paired_gz,
        EXTRACT_SINGLETONS.out.single_gz,
        INFILE_HANDLING.out.base
    )

    ch_versions = ch_versions.mix(KRAKEN_ONE.out.versions)

    // PROCESS: Run kraken2 on paired reads
    KRAKEN_TWO (
        EXTRACT_SINGLETONS.out.R1_paired_gz,
        EXTRACT_SINGLETONS.out.R2_paired_gz,
        EXTRACT_SINGLETONS.out.single_gz,
        INFILE_HANDLING.out.base
    )

    ch_versions = ch_versions.mix(KRAKEN_TWO.out.versions)
    
    // PATTERN: Collate method version information
    ch_versions.collectFile(name: 'software_versions.yml', storeDir: params.logpath)
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


