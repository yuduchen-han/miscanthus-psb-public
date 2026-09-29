# Upstream workflow（PSB，359 individuals × 8625 SNPs）

从新的空工作目录运行两个脚本；脚本位于 `03_Public_Code/00_upstream/`，默认把输出写入当前目录的 `upstream_work/`。仅发布小型脚本与最终 popmap；**不在 Git 中存放 FASTQ、BAM、Stacks catalog 或 VCF**。本文档是可复现 protocol，不是历史命令日志的替代品。依赖 Trimmomatic、FastQC、MultiQC、BWA、SAMtools、Stacks（历史 v2.68）、PLINK 1.9。

## CONFIRMED：最终 filtering provenance

`final_359_samples.tsv` 是已确认的 `popmap_table1_final359.txt` 的逐字节副本（无表头，两列 `SampleID<TAB>Region_L2`；359 人、14 组）；来源：外置历史工作目录 `scripts/03_final_v4_260902/popmap_table1_final359.txt`，**不是从下游 Genlight 反向生成**。

最终执行证据位于该目录的 `pop359.sh`、`LD359.sh` 和 `final_359_260902/` 下的 Stacks/PLINK 日志及筛选清单：

```
359-sample popmap → Stacks populations (-R 0.80 --min-mac 3 --max-obs-het 0.60 --write-single-snp --vcf)
→ 359 individuals × 9121 SNPs → PLINK --allow-extra-chr --indep-pairwise 50 5 0.4
→ prune.in 8625 / prune.out 496 → --exclude prune.out --recode vcf → 359 × 8625 final VCF
```

`02_final_snp_filtering.sh` 沿用实证 producer 的三个 PLINK 步骤，不以另一个等价算法替换。最终 `populations_final_359_8625.vcf` 是 `populations_ld_final359.pruned.vcf` 的原样副本；历史两文件及仓库 `01_Data/vcf_260902_8625_fin/populations_final_359_8625.vcf` 的 SHA-1 一致。独立 `--indep-pairwise` 步骤的日志未留存；该命令参数来自 `LD359.sh`，而保留/剔除数目有实际 `.prune.in/.prune.out` 和后续 PLINK 日志佐证。

## Historical preprocessing provenance boundary — UNRESOLVED execution

`01_preprocessing.sh` 从历史 v3 候选方法整理：paired raw FASTQ → Trimmomatic adapter/primer 与质量修剪（`ILLUMINACLIP:MIGadapter.fasta:2:30:10 MINLEN:36`；随后 `HEADCROP:17 SLIDINGWINDOW:4:15 CROP:79 MINLEN:79`，未配对 reads 同样修剪）→ 按历史方式串接四个 gzip read streams（**不是**重叠序列组装）→ FastQC/MultiQC → BWA-MEM 比对 *M. sinensis* v7.0 → sorted/indexed BAM → `ref_map.pl`/`gstacks`。早期 Stacks 必须使用独立的 **383 人 popmap**，不是此目录的最终 359 人 popmap。

上述方法见历史 `00_rawcode/99-migseq_pipeline.sh` 等候选脚本；`06_stacks_output_v3/ref_map.log` 和 `gstacks.log` 确认 `ref_map.pl`/`gstacks` 使用 383 个 BAM，但对应 v3 trimming、merging、BWA 的完整 execution logs 未全部留存。**不能声称本脚本是历史 final files 的实际 producer，也不能将这些步骤的 command-level historical provenance 标为 CONFIRMED。** 未加入旧版 VCFtools MAF 0.05 或 PLINK `50 10 0.1` 流程。

## 输入与运行（仅 protocol，未运行）

原始 paired reads 应从公共 raw-read archive 获取；参考基因组、adapter、早期 383 人 popmap、或用于单独复现 filtering 的 Stacks catalog 应从相应外部数据发布获取。**具体 archive accession / 下载地址及 v3 原始执行参数归属未核实，不能由本目录推断。** 预期输入示例：

```
raw_reads/<SampleID>_L1_1.fq.gz  +  <SampleID>_L1_2.fq.gz
reference/Msinensis_497_v7.0.fa  +  reference/MIGadapter.fasta
early_popmap.tsv                 # 383 人；由外部发布提供
upstream_work/06_stacks_output_v3/  # 完成 01 后的 catalog；单独运行 02 时须外部提供
```

```bash
EARLY_POPMAP="$PWD/early_popmap.tsv" RAW_DIR="$PWD/raw_reads" \
  REF="$PWD/reference/Msinensis_497_v7.0.fa" ADAPTER="$PWD/reference/MIGadapter.fasta" \
  WORK_DIR="$PWD/upstream_work" bash path/to/00_upstream/01_preprocessing.sh
WORK_DIR="$PWD/upstream_work" bash path/to/00_upstream/02_final_snp_filtering.sh
```

`WORK_DIR` 必须是有足够空间的**新目录**，勿指向历史正式输出。若 catalog 单独获取，可用 `STACKS_DIR` 指向其路径。最终 VCF 路径：`$WORK_DIR/final_359/populations_final_359_8625.vcf`。本目录未提供大型数据，也未执行这两个脚本或证明它们能逐位点复现历史 BAM/catalog。

## 已核实的样本历史

384 个 raw R1 文件 → 383 个早期 Stacks BAM inputs（早期 popmap 排除原名 `MIsin-FUJ-FUJ02S04`，但包含 `MIsin-FUJ-FUJ02S04-2`）→ 376 人 analysis set（`PROS_Pipeline/popmap376.txt`；TW 3、OIT 4 已在此之前移除）→ 剔除 BOS4 8、FUJ4 8、IBR01S08 1 → **359 人**。`FUJ02S04-2` 留在最终 359 人 VCF/Genlight。`docs/FACTS.md` 中“375 人”及将 `FUJ02S04-2` 排除的旧表述与正式输出冲突，**不作为此 final provenance 的依据**。从 383/376 人 popmap 到最终 359 人 popmap 的首次生成命令仍未找到。
