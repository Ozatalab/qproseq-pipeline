# qPRO-seq Shell Pipeline

This repository provides a universal, SLURM-compatible shell script for processing qPRO-seq data with optional Drosophila spike-in normalization.

---

## 🧪 Features

- Adapter trimming with Cutadapt
- First nucleotide trimming from Read 1
- Alignment to both **mm10** (mouse) and **dm6** (Drosophila spike-in)
- Filtering, sorting, indexing with SAMtools
- Spike-in read counting
- Coverage track generation with deepTools (`bamCoverage`)

---

## 🚀 Quick Start

Edit the variables at the top of `qproseq_pipeline.sh`:

```bash
input_dir="/your/path/to/fastq"
output_dir="/your/path/to/output"
bowtie2_mm10="/path/to/mm10/index/genome"
bowtie2_dm6="/path/to/dm6/index/genome"

