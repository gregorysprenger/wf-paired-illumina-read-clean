process SPADES {

    publishDir "${params.outpath}/asm",
        mode: "${params.publish_dir_mode}",
        pattern: "*"
    publishDir "${params.process_log_dir}",
        mode: "${params.publish_dir_mode}",
        pattern: ".command.*",
        saveAs: { filename -> "${task.process}${filename}"}

    container "staphb/spades@sha256:e9c50ffb4b6f0ce4d3c504dd0ce1cb3381ae942ff4d5bac24dc78119b3bfd0dd"
    

    input:
        path R1_paired_gz
        path R2_paired_gz
        path single_gz
        path outpath

    output:
        path "spades/"
        path "spades/contigs.fasta", emit: contigs


    shell:
    '''
    source bash_functions.sh

    # Get basename of input file
    base=$(basename "!{R1_paired_gz}" | cut -d _ -f 1 | sed 's/[-.,]//g')

    # Assemble with SPAdes
    failed=0
    NSLOTS=$(cat /sys/devices/system/cpu/present | cut -d '-' -f2)

    echo "INFO: Starting SPades"

    while [[ ! -f tmp/contigs.fasta ]] && [ ${failed} -lt 2 ]; do
        RAMSIZE_TOT=$(free --giga | grep '^Mem:' | awk '{print int($2*0.8)}')
        echo "INFO: RAMSIZE = ${RAMSIZE_TOT}"
        if [ ${failed} -gt 0 ]; then
            echo "ERROR: assembly file not produced by SPAdes for ${base}" >&2
            mv -f tmp/spades.log \
            tmp/"${failed}"of3-asm-attempt-failed.spades.log 2> /dev/null
            echo "INFO: SPAdes failure ${failed}; retrying assembly for ${base}" >&2
            spades.py --restart-from last -o tmp -t ${NSLOTS} >&2
        else
            spades.py --pe1-1 !{R1_paired_gz}\
            --pe1-2 !{R2_paired_gz}\
            --pe1-s !{single_gz}\
            --memory "${RAMSIZE_TOT}"\
            -o tmp --phred-offset 33\
            -t ${NSLOTS} --only-assembler >&2
        fi
        failed=$((${failed}+1))
    done

    echo "INFO: SPAdes finished"

    verify_file_minimum_size "tmp/contigs.fasta" 'SPAdes output assembly' '1M'
    if grep -E -q 'N{60}' "tmp/contigs.fasta"; then
        # avoid this again: https://github.com/ablab/spades/issues/273
        echo "ERROR: contigs.fasta contains 60+ Ns" >&2
        exit 1
    fi

    echo "INFO: Verify files"
    if [[ $(find -L !{outpath}/asm/${base}.fna -type f -size +2M 2> /dev/null) ]] && \
    [ -s !{outpath}/asm/${base}.InDels-corrected.cnt.txt ] && \
    [ -s !{outpath}/asm/${base}.SNPs-corrected.cnt.txt ] && \
    $(grep -P -q "^${base}\t" !{outpath}/qa/Summary.Illumina.CleanedReads-AlnStats.tab) && \
    $(grep -P -q "!{outpath}/asm/${base}.fna\t" !{outpath}/qa/Summary.MLST.tab) && \
    [[ $(find -L !{outpath}/annot/${base}.gbk -type f -size +3M 2> /dev/null) ]]; then
        echo "INFO: found polished assembly for ${base}" >&2
        exit 0
    fi

    gzip tmp/spades.log\
    tmp/params.txt

    mkdir spades
    mv tmp/spades.log.gz tmp/params.txt.gz tmp/contigs.fasta tmp/warnings.log tmp/assembly_graph_with_scaffolds.gfa spades/
    
    '''
}