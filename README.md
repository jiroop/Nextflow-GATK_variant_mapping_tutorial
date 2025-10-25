# Nextflow Genomics Tutorial
A simple Nextflow pipeline for learning genomics workflow management with variant calling. This follows the tutorial at: https://training.nextflow.io/basic_training/

## Description
This pipeline demonstrates core Nextflow concepts:
* Docker containers for bioinformatics tools
* Processing multiple samples in parallel
* Working with genomics file formats (BAM, VCF, GVCF)
* Channel operations and file handling
* Joint variant calling across multiple samples
* Modular pipeline design with reusable components
* Comprehensive testing with nf-test

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

## Part 3: Moving Code into Modules
Refactors the monolithic workflow into reusable, maintainable modules following best practices for pipeline organization.

### `genomics-3.nf`
A modularized version of the joint variant calling workflow:
- Extracts process definitions into separate module files
- Organizes modules in a structured directory hierarchy
- Uses `include` statements to import modules into the workflow

**Module structure**:
```
modules/
├── gatk
│   ├── haplotypecaller
│   │   └── main.nf
│   └── jointgenotyping
│       └── main.nf
└── samtools
    └── index
        └── main.nf
```

**Key features**:
- Each process lives in its own `main.nf` file within a toolkit/command directory structure
- Modules are self-contained and reusable across different workflows
- Pipeline remains fully functional with `-resume` capability intact
- Sets foundation for adding module-level tests

**What you learn**: Code organization, module creation and importing, directory structure conventions, and how modularization improves maintainability without affecting functionality.

## Part 4: Adding Tests with nf-test
Implements comprehensive testing using the nf-test framework to ensure pipeline reliability and correctness.

### `genomics-4.nf`
The same modularized workflow with a complete test suite:
- Module-level tests for each process
- Workflow-level tests for end-to-end validation
- Multiple test strategies for different scenarios

**Test structure**:
```
modules/
├── gatk
│   ├── haplotypecaller
│   │   ├── main.nf
│   │   └── tests/
│   │       └── main.nf.test
│   └── jointgenotyping
│       ├── main.nf
│       └── tests/
│           └── main.nf.test
└── samtools
    └── index
        ├── main.nf
        └── tests/
            └── main.nf.test

tests/
├── genomics-4.nf.test          # Workflow-level test
└── nextflow.config              # Test parameters
```

**Key features**:
- **Snapshot testing**: Verifies outputs match expected results across runs
- **Content assertions**: Checks specific variant calls and VCF contents
- **Setup blocks**: Uses previously generated results as test inputs
- **Co-located tests**: Each module has its tests in an adjacent `tests/` directory
- **Centralized test config**: Workflow parameters defined in `tests/nextflow.config`
- **Run all tests at once**: Single command executes entire test suite

**Testing strategies demonstrated**:
1. **SAMTOOLS_INDEX tests**: Snapshot-based verification of BAM index generation
2. **GATK_HAPLOTYPECALLER tests**: Content assertions checking specific VCF variants
3. **GATK_JOINTGENOTYPING tests**: Validation of multi-sample joint genotyping
4. **Workflow test**: End-to-end pipeline execution verification

**What you learn**: Test generation, snapshot vs. content assertions, test organization, using pre-generated data as test fixtures, workflow-level testing, and automated quality assurance for bioinformatics pipelines.

## Requirements
* Nextflow 24.10+
* Docker Desktop (Mac/Windows) or Docker Engine (Linux)
* nf-test 0.9+ (for running tests)

## Usage

### Run the pipeline
Run per-sample variant calling:
```bash
nextflow run genomics-1.nf
```

Run joint variant calling:
```bash
nextflow run genomics-2.nf
```

Run modularized version:
```bash
nextflow run genomics-3.nf
```

Run with custom output directory:
```bash
nextflow run genomics-4.nf --outdir my_results
```

Resume a previous run:
```bash
nextflow run genomics-4.nf -resume
```

Specify cohort name for joint calling:
```bash
nextflow run genomics-4.nf --cohort_name my_family
```

### Run tests

Initialize nf-test (first time only):
```bash
nf-test init
```

Run a specific module test:
```bash
nf-test test modules/samtools/index/tests/main.nf.test
```

Run a specific workflow test:
```bash
nf-test test tests/genomics-4.nf.test
```

Run all tests:
```bash
nf-test test
```

Update snapshots after expected changes:
```bash
nf-test test --update-snapshot
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

### Part 2-4 Output (`genomics-2.nf`, `genomics-3.nf`, `genomics-4.nf`)
Results in `results_genomics/` directory:
* `.bam.bai` - BAM index files
* `.g.vcf` - Genomic VCF files per sample
* `.g.vcf.idx` - GVCF index files
* `[cohort_name].joint.vcf` - Final joint-called VCF with genotypes for all samples
* `[cohort_name].joint.vcf.idx` - Joint VCF index file

### Test Output
Test results and snapshots in `.nf-test/` directory (automatically managed by nf-test)

## Understanding the Output

### Individual VCF (Part 1)
Contains variants for a single sample with genotype information for that sample only.

### Joint-Called VCF (Parts 2-4)
Contains variants across all samples with columns for each sample showing:
- **GT**: Genotype (0/0=homozygous reference, 0/1=heterozygous, 1/1=homozygous alternate)
- **AD**: Allelic depths for reference and alternate alleles
- **DP**: Total read depth at the position
- **GQ**: Genotype quality score
- **PL**: Phred-scaled likelihoods for each possible genotype

## Project Structure
```
nextflow_genomics_tutorial/
├── genomics-1.nf              # Part 1: Basic pipeline
├── genomics-2.nf              # Part 2: Joint calling
├── genomics-3.nf              # Part 3: Modularized
├── genomics-4.nf              # Part 4: With tests
├── data/
│   ├── bam/                   # Input BAM files
│   └── ref/                   # Reference genome files
├── modules/                   # Reusable process modules
│   ├── gatk/
│   │   ├── haplotypecaller/
│   │   │   ├── main.nf
│   │   │   └── tests/
│   │   └── jointgenotyping/
│   │       ├── main.nf
│   │       └── tests/
│   └── samtools/
│       └── index/
│           ├── main.nf
│           └── tests/
├── tests/
│   ├── genomics-4.nf.test     # Workflow test
│   └── nextflow.config         # Test parameters
├── nf-test.config              # nf-test configuration
├── sample_bams.txt             # List of input BAM files
└── README.md
```

## Clean Up
Remove cached files before a fresh run:
```bash
rm -rf work/ .nextflow*
```

Remove results:
```bash
rm -rf results_genomics/
```

Remove test artifacts:
```bash
rm -rf .nf-test/ tests/results/
```

## Learning Path
1. **Part 1**: Learn basic Nextflow syntax, processes, and parallel execution
2. **Part 2**: Master channel operators and complex data flows
3. **Part 3**: Apply software engineering best practices with modularization
4. **Part 4**: Implement testing for reliable, maintainable pipelines

Each part builds on the previous, creating a production-ready genomics pipeline with proper code organization and quality assurance.