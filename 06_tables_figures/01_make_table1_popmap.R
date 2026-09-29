# 科学目的：用最终 359 样本的 Region_L2 生成 Stacks Table 1 popmap。
# 输入：02_final_genlight.rds；输出：01_Data/vcf_260901_fin/popmap_table1_final359.txt。
# 流程：按 SampleID 对齐 metadata → 映射 14 个 region → 无表头 TSV。

library(adegenet)

# Load data
gl_final <- readRDS("00_rds_modules/02_final_genlight.rds")
stopifnot(nInd(gl_final)==359, nLoc(gl_final)==8625)
meta <- gl_final@other$meta
sample_ids <- indNames(gl_final)

# Prepare data: preserve the original fallback when SampleID is absent in metadata.
if ("SampleID" %in% names(meta)) {
  idx <- match(sample_ids, meta$SampleID)
  stopifnot(!anyNA(idx))
  meta <- meta[idx, , drop=FALSE]
  stopifnot(identical(as.character(meta$SampleID), as.character(sample_ids)))
} else {
  stopifnot(nrow(meta)==length(sample_ids))
  meta$SampleID <- sample_ids
}
popmap <- data.frame(SampleID=sample_ids, Region=as.character(meta$Region_L2))
stopifnot(!anyNA(popmap$Region), !any(popmap$Region==""),
          nrow(popmap)==359, length(unique(popmap$Region))==14,
          !any(grepl("-(TW|OIT)-", popmap$SampleID)))

# Save results
out_dir <- "01_Data/vcf_260901_fin"
dir.create(out_dir, recursive=TRUE, showWarnings=FALSE)
write.table(popmap, file.path(out_dir, "popmap_table1_final359.txt"),
            sep="\t", quote=FALSE, row.names=FALSE, col.names=FALSE)
