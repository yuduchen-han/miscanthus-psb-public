#!/usr/bin/env bash
# 科学目的：用 final 359 × 14-region popmap 计算 Table 1 的 Stacks 全位点统计。
# 输入：STACKS_INPUT_DIR（原分析 Stacks 输出）、已确认的 final_359_samples.tsv。
# 输出：table1_final359_stats/；流程：populations 保留原过滤参数。
set -e
: "${STACKS_INPUT_DIR:?Set STACKS_INPUT_DIR to the original Stacks input directory}"
mkdir -p 01_Data/vcf_260901_fin/table1_final359_stats
populations \
  -P "$STACKS_INPUT_DIR" \
  -M 03_Public_Code/00_upstream/final_359_samples.tsv \
  -O 01_Data/vcf_260901_fin/table1_final359_stats \
  -t 8 -R 0.80 --min-mac 3 --max-obs-het 0.60 --write-single-snp
