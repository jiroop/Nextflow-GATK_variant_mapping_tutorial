Nextflow Genomics Tutorial

A simple Nextflow pipeline for learning genomics workflow management with variant calling. This follows the tutorial at: https://training.nextflow.io/basic_training/

Description

This pipeline demonstrates core Nextflow concepts:
* Docker containers for bioinformatics tools
* Processing multiple samples in parallel
* Working with genomics file formats (BAM, VCF)
* Channel operations and file handling

What It Does

1. Index BAM files: Creates index files (.bai) for efficient BAM file access
2. Call variants: Identifies genetic variants using GATK HaplotypeCaller
3. Output VCF files: Produces variant call files for downstream analysis

Requirements

* Nextflow 25.04+
* Docker Desktop (Mac/Windows) or Docker Engine (Linux)

Usage

Run with default parameters:

nextflow run genomics-1.nf

Run with custom output directory:

nextflow run genomics-1.nf --outdir my_results

Resume a previous run:

nextflow run genomics-1.nf -resume

Input Files

Edit sample_bams.txt to specify your BAM files (one absolute path per line):

/full/path/to/sample1.bam
/full/path/to/sample2.bam

Output

Results are saved to the results/ directory:
* .bam.bai - BAM index files
* .vcf - Variant call files
* .vcf.idx - VCF index files

Clean Up

Remove cached files before a fresh run:

rm -rf work/ .nextflow*