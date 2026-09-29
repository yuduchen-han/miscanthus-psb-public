# 科学目的：以固定的 21 consensus SNP 描述两个 variety 的 RF 区分和 PCA 模式。
# 输入：最终 Genlight、consensus RDS；输出：importance、OOB 与 Figure 5 所需 PCA 数据。
# 流程：固定 SNP 提取 → 均值插补 → 3 次 ranger → importance 汇总 → 同 SNP 的 PCA。
# 本脚本是事后 characterization；绝不重新选择 candidate SNP。

library(adegenet)
library(ranger)
library(dplyr)

dir_out <- "04_Results_Archive/08c_RF_Variety_21SNP"
dir.create(dir_out, recursive=TRUE, showWarnings=FALSE)
num_trees <- 5000
rf_repeats <- 3
rf_threads <- 8
mtry_fraction <- 1/3

# Load data
gl_rf <- readRDS("00_rds_modules/02_final_genlight.rds")
consensus <- readRDS("00_rds_modules/11_Consensus_RF_pcadapt_RDA.rds")
snp_21 <- consensus$triple_snps
stopifnot(nInd(gl_rf)==359L, nLoc(gl_rf)==8625L, length(snp_21)==21L,
          all(snp_21 %in% locNames(gl_rf)),
          identical(as.character(gl_rf@other$meta$SampleID), as.character(indNames(gl_rf))))

# Prepare data: no locus selection other than the upstream fixed 21 SNPs.
meta <- gl_rf@other$meta
y_rf <- factor(meta$Species, levels=c("Miscanthus sinensis var. sinensis",
                                     "Miscanthus sinensis var. condensatus"))
gl_21 <- gl_rf[, snp_21]
gt <- as.matrix(gl_21)
rownames(gt) <- indNames(gl_21)
colnames(gt) <- locNames(gl_21)
snp_means <- colMeans(gt, na.rm=TRUE)
for (i in seq_len(ncol(gt))) {
  na_idx <- is.na(gt[, i])
  if (any(na_idx)) gt[na_idx, i] <- snp_means[i]
}
mtry_21 <- max(1, round(ncol(gt)*mtry_fraction))
stopifnot(ncol(gt)==21, mtry_21==7)

# Analysis: original seeds and permutation importance, averaged across three runs.
oob <- numeric(rf_repeats)
importance_list <- vector("list", rf_repeats)
for (r in seq_len(rf_repeats)) {
  fit <- ranger(x=gt, y=y_rf, num.trees=num_trees, mtry=mtry_21,
                importance="permutation", classification=TRUE,
                seed=500000+r, num.threads=rf_threads)
  oob[r] <- fit$prediction.error
  importance_list[[r]] <- data.frame(SNP_ID=names(fit$variable.importance),
                                      Importance=as.numeric(fit$variable.importance), Run=r)
}
importance_mean <- bind_rows(importance_list) %>% group_by(SNP_ID) %>%
  summarise(Mean_Importance=mean(Importance), SD_Importance=sd(Importance),
            .groups="drop") %>% arrange(desc(Mean_Importance))
pca_21 <- prcomp(gt, center=TRUE, scale.=TRUE)
var_exp <- pca_21$sdev^2/sum(pca_21$sdev^2)*100
pca_data <- data.frame(SampleID=rownames(gt), PC1=pca_21$x[,1], PC2=pca_21$x[,2],
                       Variety=y_rf, Habitat=meta$Habitat_new)

# Save results: 08d reads diagnostic_pca$plot_data and $variance_explained.
write.csv(importance_mean, file.path(dir_out, "Table_RF_Variety_21SNP_Importance.csv"), row.names=FALSE)
write.csv(pca_data, file.path(dir_out, "Table_RF_Variety_21SNP_PCA_Data.csv"), row.names=FALSE)
write.csv(data.frame(PC=seq_along(var_exp), Variance_Percent=var_exp),
          file.path(dir_out, "Table_RF_Variety_21SNP_PCA_Variance.csv"), row.names=FALSE)
saveRDS(list(parameters=list(num_trees=num_trees, mtry=mtry_21,
                             mtry_fraction=mtry_fraction, repeats=rf_repeats,
                             importance="permutation", threads=rf_threads),
             oob_errors=oob, mean_oob_error=mean(oob), mean_importance=importance_mean,
             diagnostic_pca=list(variance_explained=var_exp, plot_data=pca_data)),
        "00_rds_modules/08c_RF_Variety_21SNP.rds")
