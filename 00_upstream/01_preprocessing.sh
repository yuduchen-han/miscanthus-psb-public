#!/usr/bin/env bash
# 科学目的：从 paired MIG-seq FASTQ 制备有参 Stacks 调用所需 BAM/catalog。
# 输入：raw paired FASTQ、MIGadapter.fasta、M. sinensis v7.0 reference、早期 383 人 popmap。
# 输出：合并 FASTQ、FastQC/MultiQC、排序 BAM 及索引、Stacks catalog。
# 流程：去接头/引物及质量修剪 → 拼接 gzip 流（非重叠组装）→ QC → BWA-MEM → ref_map.pl/gstacks。
# UNRESOLVED：v3 FASTQ/BAM 的完整命令级执行日志未保存；本脚本是历史方法的可执行 protocol，非原执行文件。
set -euo pipefail

WORK_DIR="${WORK_DIR:-$PWD/upstream_work}"
RAW_DIR="${RAW_DIR:-$PWD/raw_reads}"
REF="${REF:-$PWD/reference/Msinensis_497_v7.0.fa}"
ADAPTER="${ADAPTER:-$PWD/reference/MIGadapter.fasta}"
EARLY_POPMAP="${EARLY_POPMAP:?Set EARLY_POPMAP to the 383-sample early Stacks popmap; do not use final_359_samples.tsv}"
THREADS=8
TRIM_DIR="$WORK_DIR/03_merged_final_v3"
QC_DIR="$WORK_DIR/04_fastqc_reports"
BAM_DIR="$WORK_DIR/05_aligned_v3"
STACKS_DIR="$WORK_DIR/06_stacks_output_v3"
TMP_DIR="$WORK_DIR/trim_tmp"

[[ -s "$EARLY_POPMAP" && -s "$REF" && -s "$ADAPTER" ]] || { echo 'Missing early popmap, reference or adapter' >&2; exit 1; }
mkdir -p "$TRIM_DIR" "$QC_DIR" "$BAM_DIR" "$STACKS_DIR" "$TMP_DIR"
export _JAVA_OPTIONS="-XX:+UseSerialGC -Xmx2g"

for r1 in "$RAW_DIR"/*_L1_1.fq.gz; do
  [[ -f "$r1" ]] || { echo "No paired R1 FASTQ in $RAW_DIR" >&2; exit 1; }
  sample="$(basename "$r1" _L1_1.fq.gz)"
  r2="$RAW_DIR/${sample}_L1_2.fq.gz"
  [[ -f "$r2" ]] || { echo "Missing R2: $r2" >&2; exit 1; }

  trimmomatic PE -threads "$THREADS" -phred33 "$r1" "$r2" \
    "$TMP_DIR/p1.fq.gz" "$TMP_DIR/u1.fq.gz" "$TMP_DIR/p2.fq.gz" "$TMP_DIR/u2.fq.gz" \
    "ILLUMINACLIP:$ADAPTER:2:30:10" MINLEN:36
  trimmomatic PE -threads "$THREADS" -phred33 "$TMP_DIR/p1.fq.gz" "$TMP_DIR/p2.fq.gz" \
    "$TMP_DIR/tp1.fq.gz" "$TMP_DIR/up1.fq.gz" "$TMP_DIR/tp2.fq.gz" "$TMP_DIR/up2.fq.gz" \
    HEADCROP:17 SLIDINGWINDOW:4:15 CROP:79 MINLEN:79
  trimmomatic SE -threads "$THREADS" -phred33 "$TMP_DIR/u1.fq.gz" "$TMP_DIR/tu1.fq.gz" \
    HEADCROP:17 SLIDINGWINDOW:4:15 CROP:79 MINLEN:79
  trimmomatic SE -threads "$THREADS" -phred33 "$TMP_DIR/u2.fq.gz" "$TMP_DIR/tu2.fq.gz" \
    HEADCROP:17 SLIDINGWINDOW:4:15 CROP:79 MINLEN:79

  cat "$TMP_DIR/tp1.fq.gz" "$TMP_DIR/tp2.fq.gz" "$TMP_DIR/tu1.fq.gz" "$TMP_DIR/tu2.fq.gz" \
    > "$TRIM_DIR/$sample.fq.gz"
  fastqc -t "$THREADS" -q -o "$QC_DIR" "$TRIM_DIR/$sample.fq.gz"
  bwa mem -t "$THREADS" -R "@RG\tID:${sample}\tSM:${sample}" "$REF" "$TRIM_DIR/$sample.fq.gz" |
    samtools view -@ "$THREADS" -Sb - |
    samtools sort -@ "$THREADS" -o "$BAM_DIR/$sample.bam"
  samtools index "$BAM_DIR/$sample.bam"
  rm -f "$TMP_DIR"/*.fq.gz
done
multiqc "$QC_DIR" -o "$QC_DIR"
rmdir "$TMP_DIR"

# ref_map.pl 调用 gstacks；早期 popmap 限定 383 个 BAM，而非最终 359 人。
ref_map.pl -T "$THREADS" --samples "$BAM_DIR" --popmap "$EARLY_POPMAP" -o "$STACKS_DIR"
