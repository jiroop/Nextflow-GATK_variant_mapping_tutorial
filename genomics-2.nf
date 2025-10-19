#!/usr/bin/env nextflow

// Pipeline parameters
 
params.outdir = "./results"
params.reads_bam = "${projectDir}/data/sample_bams.txt"

params.reference = "${projectDir}/data/ref/ref.fasta"
params.reference_index = "${projectDir}/data/ref/ref.fasta.fai" // for speeding up access to the reference file
params.reference_dict = "${projectDir}/data/ref/ref.dict" // for GATK to know metadata about each chromosome in the reference
params.intervals = "${projectDir}/data/ref/intervals.bed" // target regions of interest for variant calling
params.cohort_name = "family_trio"


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
        path "${input_bam}.g.vcf"     , emit: vcf
        path "${input_bam}.g.vcf.idx" , emit: idx // this is the index file for the g.vcf that is automatically created by GATK

    script:
    """
    gatk HaplotypeCaller \
    -R ${ref_fasta} \
    -I ${input_bam} \
    -O ${input_bam}.g.vcf \
    -L ${interval_list} \
    -ERC GVCF 

    #.g.vcf for genomic vcf
    # The -ERC GVCF switches the haplotype caller to GVCF mode for doing genomic VCF calling for joint genotyping analysis 

    """
}
    

// Combine GVCFs into GenomicsDB datastore format

process GATK_JOINTGENOTYPING {
    
    container 'quay.io/biocontainers/gatk4:4.6.2.0--py310hdfd78af_1'

    publishDir params.outdir, mode: 'symlink'

    input:
        path all_gvcfs
        path all_idxs
        path interval_list
        val cohort_name
        path ref_fasta
        path ref_index
        path ref_dict


    output:
        path "${cohort_name}_joint.vcf"         , emit: vcf    
        path "${cohort_name}_joint.vcf.idx"     , emit: idx


    script:
    def gvcfs_line = all_gvcfs.collect { gvcf -> "-V ${gvcf}" }.join(' ') // create a string with each gvcf individually listed in format "-V {gvcf1} -V {gvcf2}..."
    """
    gatk GenomicsDBImport \
        ${gvcfs_line} \
        -L ${interval_list} \
        --genomicsdb-workspace-path ${cohort_name}_gdb


    gatk GenotypeGVCFs \
        -R ${ref_fasta} \
        -V gendb://${cohort_name}_gdb \
        -L ${interval_list} \
        -O ${cohort_name}_joint.vcf

        # "gendb://" is a special gatk syntax to indicate that the input is a GenomicsDB datastore
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

    all_gvcfs_ch = GATK_HAPLOTYPECALLER.out.vcf.collect() //   collects all outputs from all parallel runs of this process, using the emitted label
    all_idxs_ch = GATK_HAPLOTYPECALLER.out.idx.collect()


    GATK_JOINTGENOTYPING(all_gvcfs_ch, all_idxs_ch, intervals_file, params.cohort_name, ref_file, ref_index_file, ref_dict_file) // create the joint genotyped VCF

}
