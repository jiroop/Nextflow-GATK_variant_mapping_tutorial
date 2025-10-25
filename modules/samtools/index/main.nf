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