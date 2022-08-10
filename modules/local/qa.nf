process QA {
    publishDir "${params.outpath}/qa",
        mode: "${params.publish_dir_mode}",
        pattern: "*.tab"
    publishDir "${params.process_log_dir}",
        mode: "${params.publish_dir_mode}",
        pattern: ".command.*",
        saveAs: { filename -> "${task.process}${filename}"}
    
    label "process_low"
    
    container "snads/quast@sha256:c8147a279feafbc88bafeeda3817ff32d43db87d31dd0978df1cd2f8022d324c"

    input:
        path base_fna
        path R1_paired_gz
        path R2_paired_gz
        path single_gz

    output:
        path "Summary.Assemblies.tab", emit: summary_assemblies
        path "Summary.Illumina.CleanedReads-Bases.tab", emit: summary_bases
        path ".command.out"
        path ".command.err"

    shell:
    '''
    # Get basename of input file
    base=$(basename "!{base_fna}" | cut -d . -f 1 | sed 's/[-.,]//g')

    NSLOTS=$(cat /sys/devices/system/cpu/present | cut -d '-' -f2)
    echo "INFO: Number of threads found: !{task.cpus}"

    quast.py --output-dir quast --min-contig 100 --threads !{task.cpus} \
    --no-html --gene-finding --gene-thresholds 300 --contig-thresholds 500,1000 \
    --ambiguity-usage one --strict-NA --silent "!{base_fna}" >&2

    mv -f quast/transposed_report.tsv Summary.Assemblies.tab

    # Count nucleotides per read set
    echo -n '' > Summary.Illumina.CleanedReads-Bases.tab
    for (( i=0; i<3; i+=3 )); do
        R1=$(basename "!{R1_paired_gz}" _R1.paired.fq.gz)
        R2=$(basename "!{R2_paired_gz}" _R2.paired.fq.gz)
        single=$(basename "!{single_gz}" _single.fq.gz)

        # Verify each set of reads groups properly
        nr_uniq_str=$(echo -e "${R1}\n${R2}\n${single}" | sort -u | wc -l)
        if [ "${nr_uniq_str}" -ne 1 ]; then
            echo "ERROR: improperly grouped ${R1} ${R2} ${single}" >&2
            exit 1
        fi
        echo -ne "${R1}\t" >> Summary.Illumina.CleanedReads-Bases.tab
        zcat "!{R1_paired_gz}" "!{R2_paired_gz}" "!{single_gz}" | \
        awk 'BEGIN{SUM=0} {if(NR%4==2){SUM+=length($0)}} END{print SUM}' \
        >> Summary.Illumina.CleanedReads-Bases.tab
    done

    '''
}