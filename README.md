# Nextflow Genomics Tutorial

A simple Nextflow pipeline for learning genomics workflow management with variant calling. This follows the tutorial at: https://training.nextflow.io/basic_training/

## Description

This pipeline demonstrates core Nextflow concepts:
* Docker containers for bioinformatics tools
* Processing multiple samples in parallel
* Working with genomics file formats (BAM, VCF)
* Channel operations and file handling

## What It Does

1. **Index BAM files**: Creates index files (.bai) for efficient BAM file access
2. **Call variants**: Identifies genetic variants using GATK HaplotypeCaller
3. **Output VCF files**: Produces variant call files for downstream analysis

## Part 1: Per-Sample Variant Calling

Builds a two-step pipeline that indexes BAM files and calls variants independently for each sample.

### `genomics-1.nf`

A linear workflow that processes multiple samples in parallel:
- **Step 1**: Generate BAM index files using Samtools
- **Step 2**: Call variants with GATK HaplotypeCaller using the indexed BAMs

**Inputs**:
- **Sample BAMs**: List of aligned sequencing files (specified in `sample_bams.txt`)
- **Reference genome**: FASTA file with index (.fai) and dictionary (.dict)
- **Intervals**: BED file defining genomic regions to analyze

**Key features**:
- Uses tuples to keep BAM files paired with their index files through the pipeline
- Reads sample list from `sample_bams.txt` for easy batch processing
- Runs each sample independently in parallel
- Handles genomics file format requirements (reference files, indexes, intervals)

**What you learn**: Process chaining, tuple usage for keeping related files together, parallel execution, and working with bioinformatics file conventions.

## Requirements

* Nextflow 25.04+
* Docker Desktop (Mac/Windows) or Docker Engine (Linux)

## Usage

Run with default parameters:
```bash
nextflow run genomics-1.nf
```

Run with custom output directory:
```bash
nextflow run genomics-1.nf --outdir my_results
```

Resume a previous run:
```bash
nextflow run genomics-1.nf -resume
```

## Input Files

Edit `sample_bams.txt` to specify your BAM files (one absolute path per line):
```
/full/path/to/sample1.bam
/full/path/to/sample2.bam
```

## Output

Results are saved to the `results_genomics/` directory:
* `.bam.bai` - BAM index files
* `.vcf` - Variant call files
* `.vcf.idx` - VCF index files

## Clean Up

Remove cached files before a fresh run:
```bash
rm -rf work/ .nextflow*
```