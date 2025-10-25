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