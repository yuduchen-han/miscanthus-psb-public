# 科学目的：以全部 8625 SNP 的 Nei GST 为背景比较固定的 21 个 consensus SNP。
# 输入：最终 Genlight、11_Consensus_RF_pcadapt_RDA.rds；输出：genome-wide GST、21 SNP 注释及稿件表。
# 流程：Species 分组 → Genlight 剂量转 Genind → 一次 Gst_Nei → 背景分位数 → 注释连接。

library(adegenet)
library(mmod)
library(readr)
library(dplyr)

# Load data
gl_final <- readRDS("00_rds_modules/02_final_genlight.rds")
consensus <- readRDS("00_rds_modules/11_Consensus_RF_pcadapt_RDA.rds")
snp_21 <- consensus$triple_snps
stopifnot(nInd(gl_final)==359L, nLoc(gl_final)==8625L,
          length(snp_21)==21L, all(snp_21 %in% locNames(gl_final)),
          identical(as.character(gl_final@other$meta$SampleID), as.character(indNames(gl_final))))

# Prepare data: preserve 0/0, 0/1, 1/1 and missing NA in the codominant Genind.
pop(gl_final) <- factor(gl_final@other$meta$Species)
geno_code <- apply(as.matrix(gl_final), 2, function(x) {
  out <- rep(NA_character_, length(x))
  out[x==0] <- "0/0"
  out[x==1] <- "0/1"
  out[x==2] <- "1/1"
  out
})
gi_all <- adegenet::df2genind(as.data.frame(geno_code), sep="/", ploidy=2,
                                type="codom", pop=pop(gl_final), NA.char="NA")
stopifnot(nInd(gi_all)==359L, nLoc(gi_all)==8625L)

# Analysis: per-locus Nei GST is computed once for the entire SNP universe.
gst_per_locus <- mmod::Gst_Nei(gi_all)$per.locus
gst_all <- tibble(SNP_ID=names(gst_per_locus), GST_Nei=as.numeric(gst_per_locus))
stopifnot(nrow(gst_all)==8625L)
background <- gst_all$GST_Nei[!is.na(gst_all$GST_Nei)]
q90 <- quantile(background, probs=0.90, names=FALSE)
q95 <- quantile(background, probs=0.95, names=FALSE)
gst_21 <- gst_all[match(snp_21, gst_all$SNP_ID), , drop=FALSE]
stopifnot(identical(gst_21$SNP_ID, snp_21))
result_full <- consensus$triple_annotation_full %>% left_join(gst_21, by="SNP_ID")
stopifnot(sum(!is.na(result_full$GST_Nei))==21L)
result_manuscript <- result_full %>%
  select(SNP_ID, Chr, SNP_Pos, GST_Nei, Gene_ID, Hit_Type, Best.hit.arabi.defline)
summary_gst <- data.frame(
  Statistic=c("Genome mean", "Genome median", "Genome 90th percentile",
              "Genome 95th percentile", "Candidate mean", "Candidate median",
              "Candidates >= genome 90th percentile", "Candidates >= genome 95th percentile"),
  Value=c(mean(background), median(background), q90, q95,
          mean(gst_21$GST_Nei, na.rm=TRUE), median(gst_21$GST_Nei, na.rm=TRUE),
          sum(gst_21$GST_Nei>=q90, na.rm=TRUE),
          sum(gst_21$GST_Nei>=q95, na.rm=TRUE)))

# Save results
out_dir <- "04_Results_Archive/14c_21SNP_GST_mmod"
dir.create(out_dir, recursive=TRUE, showWarnings=FALSE)
write_csv(gst_all, file.path(out_dir, "Table_All8625SNP_GST_mmod.csv"))
write_csv(result_full, file.path(out_dir, "Table_21SNP_GST_Annotation_Full.csv"))
write_csv(result_manuscript, file.path(out_dir, "Table_21SNP_GST_For_Manuscript.csv"))
write_csv(summary_gst, file.path(out_dir, "Summary_21SNP_vs_Genome_GST_mmod.csv"))
