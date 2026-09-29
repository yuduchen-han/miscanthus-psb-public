# 科学目的：对最终 359 × 8625 SNP 数据执行 pcadapt，获得 q<0.01 且 MAF>=0.10 的稳健候选。
# 输入：00_rds_modules/02_final_genlight.rds；输出：全位点统计、FDR/robust 表及分析 RDS。
# 流程：核对数据 → 0/1/2 转 LFMM（缺失=9）→ K=20 检查、K=3 分析 → BH → 保存。

library(adegenet)
library(pcadapt)
library(dplyr)

dir_out <- "04_Results_Archive/09_pcadapt"
dir.create(dir_out, recursive = TRUE, showWarnings = FALSE)
scree_K <- 20
optimal_K <- 3
pcadapt_min_maf <- 0.05
fdr_threshold <- 0.01
robust_maf_threshold <- 0.10

# Load data. Sample and locus order determine the identity of every p-value.
gl_main <- readRDS("00_rds_modules/02_final_genlight.rds")
stopifnot(nInd(gl_main) == 359, nLoc(gl_main) == 8625)
meta <- gl_main@other$meta
stopifnot(identical(as.character(meta$SampleID), as.character(indNames(gl_main))),
          !anyDuplicated(locNames(gl_main)))

# Prepare data: pcadapt LFMM uses 9 for missing diploid genotypes.
geno <- as.matrix(gl_main)
stopifnot(all(is.na(geno) | geno %in% 0:2))
geno[is.na(geno)] <- 9
lfmm_path <- tempfile(pattern = "pcadapt_final359_", fileext = ".lfmm")
write.table(geno, lfmm_path, quote = FALSE, row.names = FALSE, col.names = FALSE, sep = " ")
pc_input <- read.pcadapt(lfmm_path, type = "lfmm")

# Analysis: K=20 scree diagnosis documents the established K=3 selection.
pc_test <- pcadapt(pc_input, K = scree_K, min.maf = pcadapt_min_maf)
pc_final <- pcadapt(pc_input, K = optimal_K, min.maf = pcadapt_min_maf)
unlink(lfmm_path)
pvalid <- !is.na(pc_final$pvalues)
qvalues <- rep(NA_real_, length(pc_final$pvalues))
qvalues[pvalid] <- p.adjust(pc_final$pvalues[pvalid], method = "fdr")
pcadapt_statistics <- data.frame(SNP_ID = locNames(gl_main), Pvalue = pc_final$pvalues,
                                 Qvalue = qvalues, MAF = pc_final$maf)
outliers_fdr <- pcadapt_statistics %>% filter(!is.na(Qvalue), Qvalue < fdr_threshold)
outliers_robust <- outliers_fdr %>% filter(!is.na(MAF), MAF >= robust_maf_threshold)

# Save results: robust table is the input to the three-method consensus script.
write.csv(pcadapt_statistics, file.path(dir_out, "Table_pcadapt_All_SNP_Statistics.csv"),
          row.names = FALSE)
write.csv(outliers_fdr, file.path(dir_out, "Table_pcadapt_Outliers_FDR.csv"), row.names = FALSE)
write.csv(outliers_robust, file.path(dir_out, "Table_pcadapt_Outliers_Robust.csv"),
          row.names = FALSE)
saveRDS(list(settings = list(scree_K = scree_K, optimal_K = optimal_K,
                             min_maf = pcadapt_min_maf, fdr_method = "BH",
                             fdr_threshold = fdr_threshold,
                             robust_maf_threshold = robust_maf_threshold),
             scree_singular_values = pc_test$singular.values, pcadapt_final = pc_final,
             qvalues = qvalues, robust_snp_ids = outliers_robust$SNP_ID),
        "00_rds_modules/09_pcadapt_results.rds")
