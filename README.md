# Nextflow Genomics Tutorial
A simple Nextflow pipeline for learning genomics workflow management with variant calling. This follows the tutorial at: https://training.nextflow.io/basic_training/

## Description
This pipeline demonstrates core Nextflow concepts:
* Docker containers for bioinformatics tools
* Processing multiple samples in parallel
* Working with genomics file formats (BAM, VCF, GVCF)
* Channel operations and file handling
* Joint variant calling across multiple samples

## What It Does
1. **Index BAM files**: Creates index files (.bai) for efficient BAM file access
2. **Generate GVCFs**: Creates genomic VCF files for each sample containing information about all positions
3. **Combine samples**: Merges all GVCFs into a GenomicsDB datastore
4. **Joint genotyping**: Performs cohort-level variant calling to identify variants across all samples simultaneously

## Part 1: Per-Sample Variant Calling
Builds a two-step pipeline that indexes BAM files and calls variants independently for each sample.

### `genomics-1.nf`
A linear workflow that processes multiple samples in parallel:
- **Step 1**: Generate BAM index files using Samtools
- **Step 2**: Call variants with GATK HaplotypeCaller using the indexed BAMs

**Key features**:
- Uses tuples to keep BAM files paired with their index files through the pipeline
- Reads sample list from `sample_bams.txt` for easy batch processing
- Runs each sample independently in parallel
- Handles genomics file format requirements (reference files, indexes, intervals)

**What you learn**: Process chaining, tuple usage for keeping related files together, parallel execution, and working with bioinformatics file conventions.

## Part 2: Joint Variant Calling on a Cohort
Extends the pipeline to perform joint genotyping across multiple samples, providing more accurate variant calls and genotype information for all samples together.

### `genomics-2.nf`
A workflow that combines individual sample data for cohort-level analysis:
- **Step 1**: Generate BAM index files using Samtools
- **Step 2**: Create GVCFs with GATK HaplotypeCaller in GVCF mode (`-ERC GVCF`)
- **Step 3**: Combine all GVCFs into a GenomicsDB datastore
- **Step 4**: Run joint genotyping to produce cohort-level variant calls

**Key features**:
- Uses `collect()` operator to gather outputs from all samples
- Demonstrates Groovy string manipulation for dynamic command construction
- Runs multiple GATK tools in a single process
- Produces a single VCF with genotype information for all samples

**What you learn**: Channel operators (`collect()`), Groovy closures for data transformation, handling directory outputs, constructing complex command lines, and the difference between per-sample and joint variant calling.

## Requirements
* Nextflow 24.10+
* Docker Desktop (Mac/Windows) or Docker Engine (Linux)

## Usage
Run per-sample variant calling:
```bash
nextflow run genomics-1.nf
```

Run joint variant calling:
```bash
nextflow run genomics-2.nf
```

Run with custom output directory:
```bash
nextflow run genomics-2.nf --outdir my_results
```

Resume a previous run:
```bash
nextflow run genomics-2.nf -resume
```

Specify cohort name for joint calling:
```bash
nextflow run genomics-2.nf --cohort_name my_family
```

## Input Files
Edit `sample_bams.txt` to specify your BAM files (one path per line):
```
data/bam/reads_mother.bam
data/bam/reads_father.bam
data/bam/reads_son.bam
```

## Output

### Part 1 Output (`genomics-1.nf`)
Results in `results_genomics/` directory:
* `.bam.bai` - BAM index files
* `.vcf` - Individual variant call files per sample
* `.vcf.idx` - VCF index files

### Part 2 Output (`genomics-2.nf`)
Results in `results_genomics/` directory:
* `.bam.bai` - BAM index files
* `.g.vcf` - Genomic VCF files per sample
* `.g.vcf.idx` - GVCF index files
* `[cohort_name].joint.vcf` - Final joint-called VCF with genotypes for all samples
* `[cohort_name].joint.vcf.idx` - Joint VCF index file

## Understanding the Output

### Individual VCF (Part 1)
Contains variants for a single sample with genotype information for that sample only.

### Joint-Called VCF (Part 2)
Contains variants across all samples with columns for each sample showing:
- **GT**: Genotype (0/0=homozygous reference, 0/1=heterozygous, 1/1=homozygous alternate)
- **AD**: Allelic depths for reference and alternate alleles
- **DP**: Total read depth at the position
- **GQ**: Genotype quality score
- **PL**: Phred-scaled likelihoods for each possible genotype

## Clean Up
Remove cached files before a fresh run:
```bash
rm -rf work/ .nextflow*
```

Remove results:
```bash
rm -rf results_genomics/
```