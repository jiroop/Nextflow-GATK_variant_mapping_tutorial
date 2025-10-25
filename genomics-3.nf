#!/usr/bin/env nextflow

// Pipeline parameters
 
params.outdir = "./results"
params.reads_bam = "${projectDir}/data/sample_bams.txt"

params.reference = "${projectDir}/data/ref/ref.fasta"
params.reference_index = "${projectDir}/data/ref/ref.fasta.fai" // for speeding up access to the reference file
params.reference_dict = "${projectDir}/data/ref/ref.dict" // for GATK to know metadata about each chromosome in the reference
params.intervals = "${projectDir}/data/ref/intervals.bed" // target regions of interest for variant calling
params.cohort_name = "family_trio"


// Include modules

include { SAMTOOLS_INDEX } from './modules/samtools/index/main.nf'
include { GATK_HAPLOTYPECALLER } from './modules/gatk/haplotypecaller/main.nf'
include { GATK_JOINTGENOTYPING } from './modules/gatk/jointgenotyping/main.nf'


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

    all_gvcfs_ch = GATK_HAPLOTYPECALLER.out.vcf.collect() //   collects all outputs from all parallel runs of this process, using the emitted label
    all_idxs_ch = GATK_HAPLOTYPECALLER.out.idx.collect()


    GATK_JOINTGENOTYPING(all_gvcfs_ch, all_idxs_ch, intervals_file, params.cohort_name, ref_file, ref_index_file, ref_dict_file) // create the joint genotyped VCF

}
