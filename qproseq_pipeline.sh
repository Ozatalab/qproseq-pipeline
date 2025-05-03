#!/bin/bash -l

#SBATCH -A naiss2024-22-1454
#SBATCH -J qproseq_pipeline
#SBATCH -p main
#SBATCH -t 24:00:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH -e slurm-%j.err
#SBATCH -o slurm-%j.out
#SBATCH --mail-type=ALL

set -euo pipefail

# ---------------- USER CONFIGURATION ----------------

input_dir="/path/to/fastq"
output_dir="/path/to/output"
bowtie2_mm10="/path/to/mm10/index/genome"
bowtie2_dm6="/path/to/dm6/index/genome"

adapter_r1="TGGAATTCTCGGGTGCCAAGG"
adapter_r2="GATCGTCGGACTGTAGAACTCTGAAC"

# ---------------- SETUP ----------------

mkdir -p "$output_dir"/{cut,trimmed,aligned_mm10,aligned_dm6,log}

# ---------------- MAIN LOOP ----------------

for file_r1 in "$input_dir"/*_1.fastq.gz; do
  file_r2="${file_r1/_1.fastq.gz/_2.fastq.gz}"
  sample_name=$(basename "$file_r1" _1.fastq.gz)

  echo "🔹 Processing $sample_name"

  cut_r1="$output_dir/cut/${sample_name}_cut_1.fastq.gz"
  cut_r2="$output_dir/cut/${sample_name}_cut_2.fastq.gz"
  trimmed_r1="$output_dir/trimmed/${sample_name}_cut_1_trimmed.fastq"
  log_file="$output_dir/log/${sample_name}.log"
  out_dir_mm10="$output_dir/aligned_mm10/${sample_name}"
  out_dir_dm6="$output_dir/aligned_dm6/${sample_name}"
  mkdir -p "$out_dir_mm10" "$out_dir_dm6"

  # Step 1: Adapter trimming
  cutadapt -a "$adapter_r1" -A "$adapter_r2" -z -e 0.1 --minimum-length=10 \
    -o "$cut_r1" -p "$cut_r2" "$file_r1" "$file_r2" > "$log_file" 2>&1

  # Step 2: Remove first base from R1
  cutadapt -u -1 -o "$trimmed_r1" <(zcat "$cut_r1") >> "$log_file" 2>&1

  # Step 3: Align to mm10
  bowtie2 -p 16 -x "$bowtie2_mm10" \
    -1 "$trimmed_r1" -2 <(zcat "$cut_r2") \
    --local --very-sensitive --no-unal --no-mixed --no-discordant -I 10 -X 700 \
    -S "$out_dir_mm10/${sample_name}.sam" >> "$log_file" 2>&1

  # Step 4: Align to dm6 (spike-in)
  bowtie2 -p 16 -x "$bowtie2_dm6" \
    -1 "$trimmed_r1" -2 <(zcat "$cut_r2") \
    --local --very-sensitive --no-unal --no-mixed --no-discordant -I 10 -X 700 \
    -S "$out_dir_dm6/${sample_name}.sam" >> "$log_file" 2>&1

  # Step 5: Process mm10 alignment
  samtools view -F 4 -b "$out_dir_mm10/${sample_name}.sam" > "$out_dir_mm10/${sample_name}.bam"
  samtools sort "$out_dir_mm10/${sample_name}.bam" -o "$out_dir_mm10/${sample_name}_sorted.bam"
  samtools index "$out_dir_mm10/${sample_name}_sorted.bam"
  samtools view -b -q 20 -o "$out_dir_mm10/${sample_name}_uniquely_mapped.bam" "$out_dir_mm10/${sample_name}_sorted.bam"
  samtools index "$out_dir_mm10/${sample_name}_uniquely_mapped.bam"
  samtools flagstat "$out_dir_mm10/${sample_name}_uniquely_mapped.bam" > "$out_dir_mm10/${sample_name}.stat"

  # Optional: Spike-in read count (dm6)
  samtools view -c -F 4 "$out_dir_dm6/${sample_name}.sam" > "$out_dir_dm6/${sample_name}_spikein_count.txt"

  # Step 6: Generate BigWig
  bamCoverage -b "$out_dir_mm10/${sample_name}_uniquely_mapped.bam" \
    -o "$out_dir_mm10/${sample_name}.bw" \
    --normalizeUsing CPM

  echo "✅ $sample_name complete"
done

echo "🎉 All samples finished with mm10 + dm6 alignment"
