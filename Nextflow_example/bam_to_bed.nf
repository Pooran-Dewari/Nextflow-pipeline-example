#!/usr/bin/env nextflow
nextflow.enable.dsl=2

params.input = "bams/*.bam"
params.blacklist_file = file("Final_blacklist_atlantic_salmon_devmap_bodymap_combined.bed")


process BAM_TO_BED {
    publishDir 'results/beds'

    input:
    path read

    output:
    file ("${read}.bed")

    script:
    """
    
    bedtools bamtobed -i $read  | awk 'OFS="\t" {print \$1, \$2, \$3,"",\$5, \$6}' > ${read}.bed
    
    """
}


process FILTER_BLACKLIST {
    publishDir 'results/filtered_beds'

    input:
    path read
    file blacklist_file
    output:
    file ("${read}.filtered.bed")

    script:
    """
    bedtools subtract -a $read -b $blacklist_file > ${read}.filtered.bed

    """

}

workflow {
    input_ch = Channel.fromPath(params.input)
    blacklist_file = params.blacklist_file
    BAM_TO_BED(input_ch)
    FILTER_BLACKLIST(BAM_TO_BED.out.flatten(), blacklist_file)
}

