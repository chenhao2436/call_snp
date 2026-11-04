# call_snp

**A reproducible, GPU-accelerated workflow for short-read SNP discovery: FASTQ → BAM → GVCF → genome-wide joint-called VCF.**

[![Shell](https://img.shields.io/badge/shell-bash-4EAA25)](https://www.gnu.org/software/bash/)
[![Workflow](https://img.shields.io/badge/workflow-GATK%20%7C%20Parabricks-blue)](https://github.com/chenhao2436/call_snp)
[![Status](https://img.shields.io/badge/status-alpha-orange)](https://github.com/chenhao2436/call_snp)

**call_snp** is a command-line workflow for processing paired-end sequencing data and performing multi-sample joint genotyping. It currently contains two Bash scripts:

1. `fq2gvcf_sample.sh`: align paired-end FASTQ reads and generate one GVCF per sample using NVIDIA Parabricks.
2. `gvcf2vcf_genome.sh`: import sample GVCFs into GenomicsDB, run GATK joint genotyping by reference contig, and gather the shards into a genome-wide VCF.

> **Development status: alpha / not yet validated as a turnkey release.** The repository is being prepared for reproducible use. The current FASTQ-to-GVCF script has issues that must be fixed before production use (see [Known issues and validation status](#known-issues-and-validation-status)). Do not launch a large production run using the current script without reviewing and correcting these issues.

## Contents

- [Workflow overview](#workflow-overview)
- [Software requirements](#software-requirements)
- [Input files](#input-files)
- [Installation](#installation)
- [Step 1: FASTQ to GVCF](#step-1-fastq-to-gvcf)
- [Step 2: GVCF to genome-wide VCF](#step-2-gvcf-to-genome-wide-vcf)
- [Output structure](#output-structure)
- [Resource planning](#resource-planning)
- [Reproducibility and quality control](#reproducibility-and-quality-control)
- [Known issues and validation status](#known-issues-and-validation-status)
- [Troubleshooting](#troubleshooting)
- [Roadmap](#roadmap)
- [Citation and contact](#citation-and-contact)

## Workflow overview

```text
Paired-end FASTQ files
        |
        v
Reference preparation
(BWA index + FASTA .fai + GATK sequence dictionary)
        |
        v
Parabricks fq2bam
        |
        v
Coordinate-sorted BAM
        |
        v
Parabricks HaplotypeCaller (--gvcf)
        |
        v
One GVCF per sample + index
        |
        v
GVCF sample map (sample<TAB>absolute_path)
        |
        v
GATK GenomicsDBImport (per contig)
        |
        v
GATK GenotypeGVCFs (per contig)
        |
        v
GatherVcfs
        |
        v
Genome-wide joint-called VCF + index
```

The two stages can be run separately. This is useful when read alignment and GVCF generation are performed on a GPU server, while joint genotyping is performed on a separate compute server.

## Software requirements

### Stage 1: FASTQ to GVCF

- Linux server.
- Bash.
- Docker Engine.
- NVIDIA GPU and a working NVIDIA Container Toolkit / Docker GPU runtime.
- NVIDIA Clara Parabricks container. The current script refers to `nvcr.io/nvidia/clara/clara-parabricks:4.1.1-1`.
- BWA, for reference indexing.
- SAMtools, for FASTA indexing (the current script's indexing command needs correction; see below).
- GATK, for creating the reference sequence dictionary.

### Stage 2: GVCF to VCF

- Linux server.
- Bash 4 or later (the script uses associative arrays and `mapfile`).
- GATK 4 with `GenomicsDBImport`, `GenotypeGVCFs`, `GatherVcfs`, and `IndexFeatureFile`.
- Java runtime compatible with the installed GATK release.
- Standard command-line utilities: `awk`, `sed`, `tee`, `realpath`, `wc`, `date`, and `mkdir`.

The Conda file in this repository provides the host-side command-line environment. Parabricks is distributed and run through Docker, not installed by pip or Conda.

## Input files

### Reference genome

Provide a reference FASTA, for example:

```text
reference/
├── reference.fa
├── reference.fa.fai       # created by samtools faidx
├── reference.dict        # created by gatk CreateSequenceDictionary
└── BWA index files       # created by bwa index
```

All samples in one joint-calling run must use the **same reference assembly and contig naming convention**. Do not mix references or assemblies in one sample map.

### Paired-end FASTQ files

The current Stage 1 script expects filenames in this pattern:

```text
{sample}_1_clean.fq.gz
{sample}_2_clean.fq.gz
```

The sample identifier must not contain whitespace or tabs. The two files must be paired reads from the same sample.

### Sample list for Stage 1

One sample ID per line, for example:

```text
Parent_A
Parent_B
Pool_Long
Pool_Round
```

Blank lines and comment handling are not currently implemented in Stage 1. Use a clean list containing only sample IDs.

### Sample map for Stage 2

Stage 2 requires a tab-delimited, two-column file:

```text
Parent_A    /absolute/path/to/Parent_A.g.vcf.gz
Parent_B    /absolute/path/to/Parent_B.g.vcf.gz
Pool_Long   /absolute/path/to/Pool_Long.g.vcf.gz
Pool_Round  /absolute/path/to/Pool_Round.g.vcf.gz
```

Use a literal TAB between the sample name and the path, not spaces. Each path must point to an existing, readable, indexed GVCF. Sample names must be unique.

Example generation command (run in the directory containing the GVCFs):

```bash
for f in *.g.vcf.gz; do
    sample=$(basename "$f" .g.vcf.gz)
    printf '%s\t%s\n' "$sample" "$(realpath "$f")"
done > gvcf_sample.list
```

Review the resulting list before starting joint genotyping.

## Installation

Clone the repository:

```bash
git clone https://github.com/chenhao2436/call_snp.git
cd call_snp
```

Create the host-side Conda environment:

```bash
conda env create -f environment.yml
conda activate call-snp
```

Check the main tools:

```bash
bash --version
bwa 2>&1 | head
samtools --version
gatk --version
docker --version
```

For GPU-enabled Docker, verify that the host and container runtime can see the GPU:

```bash
nvidia-smi
docker run --rm --gpus all nvcr.io/nvidia/clara/clara-parabricks:4.1.1-1 nvidia-smi
```

The Parabricks image may require NVIDIA NGC access or prior image authentication. Follow NVIDIA's current container access instructions.

## Step 1: FASTQ to GVCF

### 1. Prepare directories and configuration

The current script contains hard-coded server paths near the top of `fq2gvcf_sample.sh`. Edit these paths to match your environment:

- Reference FASTA.
- FASTQ input directory.
- BAM output directory.
- GVCF output directory.
- Sample list.

The script currently reads the sample list from its **first positional argument** (`$1`). Pass the sample-list path explicitly.

### 2. Prepare reference indexes

The required reference resources are:

```bash
bwa index /path/to/reference.fa
samtools faidx /path/to/reference.fa
gatk CreateSequenceDictionary -R /path/to/reference.fa
```

The current Stage 1 script attempts to run `samtools index` on the FASTA. This is not the correct SAMtools command for creating a FASTA index. Correct this before execution.

### 3. Run the workflow

After correcting the issues listed below and updating all paths:

```bash
bash fq2gvcf_sample.sh /path/to/sample.list
```

The intended per-sample operations are:

- Parabricks `fq2bam`: paired-end FASTQ → aligned BAM.
- Parabricks `haplotypecaller --gvcf`: BAM → single-sample GVCF.

Expected outputs:

```text
bam/
├── Parent_A.bam
└── Parent_B.bam

gvcf/
├── Parent_A.g.vcf.gz
├── Parent_B.g.vcf.gz
└── gvcf_sample.list
```

The exact output locations are determined by the paths configured in the script.

**Important:** The current script's second Docker command mounts the GVCF directory at `/3` but also uses `/3/${sample}.bam` as its BAM input path. The BAM input should refer to the BAM directory mount. This must be fixed before Stage 1 can run correctly.

## Step 2: GVCF to genome-wide VCF

### 1. Prepare and validate GVCFs

Every GVCF must be complete, readable, generated against the same reference, and indexed. For compressed GVCFs, a `.tbi`, `.idx`, or `.csi` index is expected by the current script.

If an index is missing, create it with GATK:

```bash
gatk IndexFeatureFile -I sample.g.vcf.gz
```

### 2. Build the sample map

Create a two-column, tab-delimited sample map as described in [Input files](#input-files). Use absolute paths to avoid working-directory ambiguity.

### 3. Run joint genotyping

The script accepts three positional arguments:

```bash
bash gvcf2vcf_genome.sh [gvcf_list] [reference.fa] [output_directory]
```

Example:

```bash
bash gvcf2vcf_genome.sh \
    /data/project/gvcf/gvcf_sample.list \
    /data/reference/reference.fa \
    /data/project/joint_vcf
```

If arguments are omitted, the script falls back to paths hard-coded in its configuration section. For portable use, pass all three arguments explicitly.

### 4. Configure compute resources

The following variables are configured near the top of `gvcf2vcf_genome.sh`:

| Variable | Meaning | Current default |
|---|---|---:|
| `MAX_JOBS` | Number of contigs processed concurrently | 3 |
| `JAVA_XMX_IMPORT` | Maximum heap for GenomicsDBImport | 20g |
| `JAVA_XMX_GENOTYPE` | Maximum heap for GenotypeGVCFs | 28g |
| `JAVA_XMX_GATHER` | Maximum heap for GatherVcfs | 12g |
| `READER_THREADS` | GVCF reader threads | 4 |
| `INCLUDE_NON_VARIANT_SITES` | Include non-variant sites | no |

Adjust these values to the available memory, CPU, storage throughput, and sample count. The defaults are not universal recommendations.

### 5. Output

The final file name is constructed from `qunti_name` and `REF_NAME` in the script. The current defaults are:

```text
F2_3_PP64.genome.joint.vcf.gz
```

The script also creates logs, per-contig VCF shards, status flags, temporary directories, a sample map, and interval metadata under the output directory.

## Output structure

The Stage 2 output is organized approximately as follows:

```text
output_directory/
├── logs/
│   ├── *.run_*.log
│   └── *.log
├── status/
│   ├── *.ok
│   └── *.fail
├── per_interval/
│   └── <interval>/
│       └── <interval>.joint.vcf.gz
├── genomicsdb/
├── tmp/
├── *.sample_map
├── *.intervals.tsv
├── *.shards.list
└── <population>_<reference>.genome.joint.vcf.gz
```

The script uses per-contig processing and status files to support skipping completed shards on a subsequent run. Review the log and status directory before treating a run as successful.

## Resource planning

Joint genotyping can be memory- and I/O-intensive. Before running a large cohort:

- Estimate the total size of all GVCFs and available scratch space.
- Confirm that the output filesystem has enough free space.
- Start with a low `MAX_JOBS` value and increase only after monitoring memory and I/O.
- Keep temporary GenomicsDB workspaces on fast local storage when possible.
- Avoid running multiple large joint-genotyping jobs on the same storage system without resource coordination.

The script attempts to use `/dev/shm` for scratch and falls back to the output directory's `tmp` folder. Ensure the selected scratch filesystem has sufficient capacity.

## Reproducibility and quality control

For each analysis, record:

- Git commit or release tag of `call_snp`.
- Reference genome name, version, and checksum.
- Parabricks container tag and digest.
- GATK, BWA, SAMtools, Java, Docker, and driver versions.
- Sample-list and sample-map files.
- Key resource settings and the complete run log.

Before downstream analysis, verify:

1. Every expected sample has a non-empty GVCF.
2. Every GVCF has a valid index.
3. The sample map contains the expected unique sample IDs.
4. The final VCF is non-empty, readable, and indexed.
5. The VCF header contains the expected sample columns.
6. The contig names and order are consistent with the reference.
7. The number of samples and key variant statistics are biologically plausible.

This repository currently performs variant calling and joint genotyping. It does **not** yet implement variant filtering, SNP-only extraction, BSA analysis, GWAS, functional annotation, or expression analysis.

## Known issues and validation status

The following items were identified by static inspection of the current scripts. They should be addressed before calling the workflow production-ready.

| Priority | File | Issue | Required action |
|---|---|---|---|
| Critical | `fq2gvcf_sample.sh` | The standalone text `基础配置` is not a comment and may be executed as a command. | Change it to a comment. |
| Critical | `fq2gvcf_sample.sh` | `samtools index` is used on a FASTA file. | Use `samtools faidx`. |
| Critical | `fq2gvcf_sample.sh` | The HaplotypeCaller Docker command points its BAM input to the GVCF mount. | Mount the BAM directory separately and use the correct container path. |
| High | `fq2gvcf_sample.sh` | Most paths and the Parabricks image are hard-coded. | Add validated command-line options or a configuration file. |
| High | `fq2gvcf_sample.sh` | The script does not robustly validate FASTQ pairs, sample IDs, tool availability, or output files. | Add preflight checks, strict shell options, logging, and per-sample failure handling. |
| High | `gvcf2vcf_genome.sh` | Population and reference labels are hard-coded. | Expose them as command-line options or configuration parameters. |
| High | `gvcf2vcf_genome.sh` | The final `mv` command is unquoted and redundant because the destination is already the configured final path. | Remove it or replace it with a validated, quoted operation. |
| Medium | Both | No automated test dataset or end-to-end test is included. | Add a small test dataset and documented expected outputs. |

**Validation statement:** This README is based on static inspection of the repository scripts. No end-to-end run, benchmark, or biological validation is claimed.

## Troubleshooting

### Docker cannot access the GPU

Check the host GPU driver and Docker GPU runtime:

```bash
nvidia-smi
docker run --rm --gpus all <parabricks-image> nvidia-smi
```

### GATK reports a missing reference index or dictionary

Create the FASTA index and sequence dictionary:

```bash
samtools faidx reference.fa
gatk CreateSequenceDictionary -R reference.fa
```

### GATK reports a missing GVCF index

```bash
gatk IndexFeatureFile -I sample.g.vcf.gz
```

### Joint genotyping fails for one contig

Inspect the corresponding file under `logs/` and the matching `.fail` status flag. Confirm that every GVCF is readable, indexed, and based on the same reference. Check scratch space and Java memory settings before retrying.

### The final VCF is missing

The script only gathers shards after all contigs have completed successfully. Check the master log, per-contig logs, and `status/` directory. Do not treat the existence of a partial shard as a completed run.

## Roadmap

The repository's intended development direction is:

- [x] FASTQ-to-GVCF workflow (initial implementation; requires fixes and validation).
- [x] GVCF-to-genome-wide-VCF workflow (initial implementation; requires validation).
- [ ] Replace hard-coded paths with a unified configuration interface.
- [ ] Add preflight validation and structured logs to both stages.
- [ ] Add a small reproducible test dataset and automated shell tests.
- [ ] Add release tags and versioned environment documentation.
- [ ] Add BSA visualization and interactive plotting modules.
- [ ] Add transcriptome gene-ID conversion and plotting utilities.
- [ ] Add variant-expression integration and functional annotation modules.

The downstream modules above are planned work, not features currently included in this repository.

## Citation and contact

If you use this workflow, please cite the software repository and the underlying tools used in your analysis:

- Repository: https://github.com/chenhao2436/call_snp
- NVIDIA Parabricks: https://docs.nvidia.com/clara/parabricks/
- GATK: https://gatk.broadinstitute.org/
- BWA: https://github.com/lh3/bwa
- SAMtools: https://www.htslib.org/

For bugs and feature requests, please open a GitHub Issue with the command used, relevant log excerpts, software versions, and a minimal description of the input layout. Do not upload human genomic data, credentials, or other sensitive data to public issues.
