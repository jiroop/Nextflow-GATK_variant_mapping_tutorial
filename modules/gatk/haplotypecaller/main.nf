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