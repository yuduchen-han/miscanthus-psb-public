#!/usr/bin/env bash
# 科学目的：从早期 Stacks catalog 和已确认的 359 人 popmap 生成论文最终 8625 SNP VCF。
# 输入：06_stacks_output_v3 catalog、同目录 final_359_samples.tsv。
# 输出：359 人/9121 SNP Stacks VCF、PLINK prune.in/out、359 人/8625 SNP final VCF。
# 流程：populations → VCF 转 PLINK → LD 筛选 → 从原始 9121 SNP 集排除 prune.out → 保留最终 VCF。
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
WORK_DIR="${WORK_DIR:-$PWD/upstream_work}"
STACKS_DIR="${STACKS_DIR:-$WORK_DIR/06_stacks_output_v3}"
FINAL_DIR="$WORK_DIR/final_359"
POPMAP="$SCRIPT_DIR/final_359_samples.tsv"
PREFIX="$FINAL_DIR/populations_ld_final359"
THREADS=8

[[ -s "$POPMAP" && -d "$STACKS_DIR" ]] || { echo 'Missing final popmap or Stacks catalog' >&2; exit 1; }
[[ "$(wc -l < "$POPMAP" | tr -d ' ')" -eq 359 ]] || { echo 'Expected 359 popmap rows' >&2; exit 1; }
[[ "$(awk -F '\t' '{print $2}' "$POPMAP" | sort -u | wc -l | tr -d ' ')" -eq 14 ]] || { echo 'Expected 14 regions' >&2; exit 1; }
mkdir -p "$FINAL_DIR"

populations -P "$STACKS_DIR" -M "$POPMAP" -t "$THREADS" \
  -R 0.80 --min-mac 3 --max-obs-het 0.60 --write-single-snp --vcf -O "$FINAL_DIR"

count_snps() { awk '!/^#/ {n++} END {print n+0}' "$1"; }
count_samples() { awk -F '\t' '$1=="#CHROM" {print NF-9; exit}' "$1"; }
[[ "$(count_samples "$FINAL_DIR/populations.snps.vcf")" -eq 359 &&
   "$(count_snps "$FINAL_DIR/populations.snps.vcf")" -eq 9121 ]] || { echo 'Expected 359 samples and 9121 pre-LD SNPs' >&2; exit 1; }

plink --vcf "$FINAL_DIR/populations.snps.vcf" --allow-extra-chr --make-bed --out "$PREFIX"
plink --bfile "$PREFIX" --allow-extra-chr --indep-pairwise 50 5 0.4 --out "${PREFIX}.pruned"
[[ "$(wc -l < "${PREFIX}.pruned.prune.in" | tr -d ' ')" -eq 8625 &&
   "$(wc -l < "${PREFIX}.pruned.prune.out" | tr -d ' ')" -eq 496 ]] || { echo 'Expected 8625 retained and 496 removed SNPs' >&2; exit 1; }

# 与已确认的 final PLINK producer 一致：从 9121 SNP bfile 中排除 prune.out，而非另行提取 prune.in。
plink --bfile "$PREFIX" --allow-extra-chr --exclude "${PREFIX}.pruned.prune.out" \
  --recode vcf --out "${PREFIX}.pruned"
[[ "$(count_samples "${PREFIX}.pruned.vcf")" -eq 359 &&
   "$(count_snps "${PREFIX}.pruned.vcf")" -eq 8625 ]] || { echo 'Expected 359 samples and 8625 final SNPs' >&2; exit 1; }
# 已核实历史 final VCF 与 PLINK --recode vcf 文件逐字节相同；不改写染色体、SNP ID 或样本名。
cp "${PREFIX}.pruned.vcf" "$FINAL_DIR/populations_final_359_8625.vcf"
