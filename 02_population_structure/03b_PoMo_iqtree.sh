#!/usr/bin/env bash
set -e

# 03b_PoMo_iqtree.sh
#
# 科学目的：
# 基于46个Pop_L3的all-sites PoMo COUNTSFILE构建population-level ML tree。
#
# 输入：
# 04_Results_Archive/03_ML_Tree/PoMo_population/populations.cf
#
# 模型：
# GTR+P：PoMo模型
# N=9：IQ-TREE默认virtual population size
# weighted binomial sampling：IQ-TREE默认
# 1000次Ultrafast Bootstrap评估节点支持度

IQTREE="01_Data/iqtree-2.4.0-macOS/bin/iqtree2"
INPUT="04_Results_Archive/03_ML_Tree/PoMo_population/populations.cf"
PREFIX="04_Results_Archive/03_ML_Tree/PoMo_population/population_PoMo"

"$IQTREE" \
  -s "$INPUT" \
  -m GTR+P \
  -bb 1000 \
  -nt AUTO \
  -pre "$PREFIX"
