# 科学目的：诊断 RF-Habitat 的 mtry / ntree；不筛选候选 SNP。
# 输入：最终 359 × 8625 Genlight（Habitat_new）；输出：各参数 OOB、importance 可重复性表和 PDF。
# 流程：逐 SNP 均值插补 → 每组参数独立 3 次 ranger → OOB 与 permutation importance 相关性。
# mtry 与 ntree 两阶段分别运行；selected_mtry=862 是既有诊断后人工选择，非自动最优值。

library(ranger)
library(dplyr)
library(ggplot2)
library(adegenet)

dir_out <- "04_Results_Archive/08a_RF_Habitat_Parameter_Diagnostic"
dir.create(dir_out, recursive=TRUE, showWarnings=FALSE)
run_stage <- "ntree" # 改为 "mtry" 可生成 Stage A 的诊断；两阶段不同时运行。
mtry_diagnostic_ntree <- 2000
n_repeats <- 3
selected_mtry <- 862

# Load data
gl_main <- readRDS("00_rds_modules/02_final_genlight.rds")
stopifnot(nInd(gl_main)==359, nLoc(gl_main)==8625)
y_rf <- factor(gl_main@other$meta$Habitat_new, levels=c("Inland", "Seashore", "Island"))

# Prepare data: 与历史诊断一致的逐 SNP 均值插补。
gt <- as.matrix(gl_main)
snp_means <- colMeans(gt, na.rm=TRUE)
for (i in seq_len(ncol(gt))) {
  missing_idx <- is.na(gt[, i])
  if (any(missing_idx)) gt[missing_idx, i] <- snp_means[i]
}

# Analysis: 每组参数固定 3 个独立 seed；两两 Pearson correlation 衡量 importance 可重复性。
run_repeated_rf <- function(gt, y, num_trees, mtry_value, n_repeats) {
  oob <- numeric(n_repeats)
  importance_list <- vector("list", n_repeats)
  for (r in seq_len(n_repeats)) {
    rf_fit <- ranger(x=gt, y=y, num.trees=num_trees, mtry=mtry_value,
                     importance="permutation", classification=TRUE,
                     seed=202600+r, num.threads=8)
    oob[r] <- rf_fit$prediction.error
    importance_list[[r]] <- rf_fit$variable.importance
  }
  importance_matrix <- do.call(cbind, importance_list)
  importance_cor <- cor(importance_matrix, method="pearson")
  pairwise_cor <- importance_cor[upper.tri(importance_cor)]
  list(mean_oob=mean(oob), sd_oob=sd(oob),
       mean_importance_cor=mean(pairwise_cor), min_importance_cor=min(pairwise_cor))
}

if (run_stage == "mtry") {
  p <- ncol(gt)
  mtry_values <- unique(round(c(sqrt(p), 2*sqrt(p), 0.1*p, 0.2*p)))
  mtry_df <- bind_rows(lapply(mtry_values, function(mtry_value) {
    result <- run_repeated_rf(gt, y_rf, mtry_diagnostic_ntree, mtry_value, n_repeats)
    data.frame(Mtry=mtry_value, Num_Trees=mtry_diagnostic_ntree,
               Mean_OOB_Error=result$mean_oob, SD_OOB_Error=result$sd_oob,
               Mean_Importance_Correlation=result$mean_importance_cor,
               Min_Importance_Correlation=result$min_importance_cor)
  }))
  write.csv(mtry_df, file.path(dir_out, "Table_RF_mtry_Diagnostic.csv"), row.names=FALSE)
  p_mtry_oob <- ggplot(mtry_df, aes(Mtry, Mean_OOB_Error)) +
    geom_line() + geom_point(size=3) +
    geom_errorbar(aes(ymin=Mean_OOB_Error-SD_OOB_Error,
                      ymax=Mean_OOB_Error+SD_OOB_Error), width=30) +
    theme_classic(base_size=14) + labs(x="mtry", y="Mean OOB prediction error")
  ggsave(file.path(dir_out, "Figure_RF_mtry_OOB.pdf"), p_mtry_oob, width=7, height=5)
  p_mtry_cor <- ggplot(mtry_df, aes(Mtry, Mean_Importance_Correlation)) +
    geom_line() + geom_point(size=3) + theme_classic(base_size=14) +
    labs(x="mtry", y="Mean correlation of permutation importance")
  ggsave(file.path(dir_out, "Figure_RF_mtry_Importance_Repeatability.pdf"),
         p_mtry_cor, width=7, height=5)
}

if (run_stage == "ntree") {
  if (is.na(selected_mtry)) stop("Run Stage A first and set selected_mtry.")
  ntree_values <- c(500, 1000, 2000, 5000, 10000)
  ntree_df <- bind_rows(lapply(ntree_values, function(num_trees) {
    result <- run_repeated_rf(gt, y_rf, num_trees, selected_mtry, n_repeats)
    data.frame(Num_Trees=num_trees, Mtry=selected_mtry,
               Mean_OOB_Error=result$mean_oob, SD_OOB_Error=result$sd_oob,
               Mean_Importance_Correlation=result$mean_importance_cor,
               Min_Importance_Correlation=result$min_importance_cor)
  }))
  write.csv(ntree_df, file.path(dir_out, "Table_RF_ntree_Diagnostic.csv"), row.names=FALSE)
  p_ntree_oob <- ggplot(ntree_df, aes(Num_Trees, Mean_OOB_Error)) +
    geom_line() + geom_point(size=3) +
    geom_errorbar(aes(ymin=Mean_OOB_Error-SD_OOB_Error,
                      ymax=Mean_OOB_Error+SD_OOB_Error), width=100) +
    theme_classic(base_size=14) + labs(x="Number of trees", y="Mean OOB prediction error")
  ggsave(file.path(dir_out, "Figure_RF_ntree_OOB.pdf"), p_ntree_oob, width=7, height=5)
  p_ntree_cor <- ggplot(ntree_df, aes(Num_Trees, Mean_Importance_Correlation)) +
    geom_line() + geom_point(size=3) + theme_classic(base_size=14) +
    labs(x="Number of trees", y="Mean correlation of permutation importance")
  ggsave(file.path(dir_out, "Figure_RF_ntree_Importance_Repeatability.pdf"),
         p_ntree_cor, width=7, height=5)
}
