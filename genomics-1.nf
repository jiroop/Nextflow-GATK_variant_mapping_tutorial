#!/usr/bin/env nextflow


// Pipeline parameters
 
params.outdir = "./results"
params.reads_bam = "${projectDir}/data/sample_bams.txt"

params.reference = "${projectDir}/data/ref/ref.fasta"
params.reference_index = "${projectDir}/data/ref/ref.fasta.fai"
params.reference_dict = "${projectDir}/data/ref/ref.dict"
params.intervals = "${projectDir}/data/ref/intervals.bed"


// Primary input


// Generate BAM index file
process SAMTOOLS_INDEX {

    container 'quay.io/biocontainers/samtools:1.20--h50ea8bc_1'

    publishDir params.outdir, mode: 'symlink'

    input:
        path input_bam

    output:
        tuple path(input_bam), path("${input_bam}.bai") 

    script:
    """
    samtools index "${input_bam}"
    """

}

// Call variants with GATK HaplotypeCaller

process GATK_HAPLOTYPECALLER {

    container 'quay.io/biocontainers/gatk4:4.6.2.0--py310hdfd78af_1'

    publishDir params.outdir, mode: 'symlink'

    input:
        tuple path (input_bam), path(input_bam_index)
        path ref_fasta
        path ref_index
        path ref_dict
        path interval_list

    output:
        path "${input_bam}.vcf"     , emit: vcf
        path "${input_bam}.vcf.idx" , emit: vcf_index

    script:
    """
    gatk HaplotypeCaller \
    -R ${ref_fasta} \
    -I ${input_bam} \
    -O ${input_bam}.vcf \
    -L ${interval_list}
    """

}
    
workflow {

    // Create input channel

    reads_ch = Channel
        .fromPath(params.reads_bam)
        .view()
        .splitText()
        .map { it.trim() }
        .map { file("${projectDir}/${it}") }
        .view()
    

    // Load the file paths for accessory files

    ref_file = file(params.reference)
    ref_index_file = file(params.reference_index)
    ref_dict_file = file(params.reference_dict)
    intervals_file = file(params.intervals)


    // Create index file for input BAM file

    SAMTOOLS_INDEX(reads_ch)

    GATK_HAPLOTYPECALLER(SAMTOOLS_INDEX.out, ref_file, ref_index_file, ref_dict_file, intervals_file)

}
