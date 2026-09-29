# 科学目的：用个体基因型与海岸距离/Bio 环境变量的 RDA 寻找环境关联候选 SNP。
# 输入：最终 Genlight (359 × 8625) + publication sample_metadata.csv；输出：模型检验、图和候选表。
# 流程：按 SampleID 对齐 → 基因型均值插补 → 环境变量 NA/VIF 过滤 → RDA → permutation → 3 SD 候选。

library(adegenet)
library(dplyr)
library(vegan)
library(ggplot2)
library(ggrepel)
library(usdm)

dir_out <- "04_Results_Archive/10_RDA_HabitatNew"
dir.create(dir_out, recursive = TRUE, showWarnings = FALSE)
vif_threshold <- 10
n_perm <- 999
parallel_cores <- 1
outlier_sd_threshold <- 3

# Load data: never join genotype and environment by implicit row position.
gl_main <- readRDS("00_rds_modules/02_final_genlight.rds")
meta_gen <- gl_main@other$meta
meta_env <- read.csv("03_Public_Code/data/metadata/sample_metadata.csv", stringsAsFactors = FALSE)
stopifnot(nInd(gl_main) == 359, nLoc(gl_main) == 8625,
          !anyDuplicated(meta_gen$SampleID), !anyDuplicated(meta_env$SampleID),
          identical(as.character(meta_gen$SampleID), as.character(indNames(gl_main))))
idx_env <- match(indNames(gl_main), meta_env$SampleID)
stopifnot(!anyNA(idx_env))
meta_env_final <- meta_env[idx_env, , drop = FALSE]
stopifnot(identical(as.character(meta_env_final$SampleID), as.character(indNames(gl_main))))

# Prepare genotype: discard all-NA loci; impute each remaining locus by its mean;
# then discard zero-variance loci. SNP IDs remain attached to columns.
gt_raw <- as.matrix(gl_main)
rownames(gt_raw) <- indNames(gl_main)
colnames(gt_raw) <- locNames(gl_main)
stopifnot(!anyDuplicated(colnames(gt_raw)))
all_na_snp <- colSums(!is.na(gt_raw)) == 0
gt <- gt_raw[, !all_na_snp, drop = FALSE]
stopifnot(ncol(gt) >= 2)
snp_means <- colMeans(gt, na.rm = TRUE)
for (i in seq_len(ncol(gt))) {
  na_idx <- is.na(gt[, i])
  if (any(na_idx)) gt[na_idx, i] <- snp_means[i]
}
snp_sd <- apply(gt, 2, sd)
zero_var_snp <- is.na(snp_sd) | snp_sd == 0
gt <- gt[, !zero_var_snp, drop = FALSE]
stopifnot(ncol(gt) >= 2)

# Prepare predictors: remove variables containing NA (no environmental imputation),
# apply usdm VIF selection, then remove zero-variance variables and standardize.
bio_cols <- grep("^Bio", names(meta_env_final), value = TRUE)
stopifnot(length(bio_cols) > 0)
env_all <- meta_env_final %>% select(distance_to_coast_km, all_of(bio_cols))
env_all <- as.data.frame(lapply(env_all, function(x) as.numeric(as.character(x))))
rownames(env_all) <- meta_env_final$SampleID
env_complete <- env_all[, colSums(is.na(env_all)) == 0, drop = FALSE]
stopifnot(ncol(env_complete) >= 2)
vif_res <- usdm::vifstep(env_complete, th = vif_threshold)
env_core <- usdm::exclude(env_complete, vif_res)
env_sd <- apply(env_core, 2, sd, na.rm = TRUE)
env_core <- env_core[, !(is.na(env_sd) | env_sd == 0), drop = FALSE]
stopifnot(ncol(env_core) >= 2)
env_scaled <- scale(env_core)
stopifnot(identical(rownames(gt), rownames(env_scaled)))

# Analysis: full model and its permutation test; species scores define candidates.
rda_model <- vegan::rda(gt ~ ., data = as.data.frame(env_scaled))
set.seed(2026)
signif_full <- vegan::anova.cca(rda_model, permutations = n_perm, parallel = parallel_cores)
R2 <- vegan::RsquareAdj(rda_model)
snp_scores <- vegan::scores(rda_model, choices = c(1, 2), display = "species")
df_snps <- data.frame(snp_scores, SNP_ID = rownames(snp_scores))
lim1 <- outlier_sd_threshold * sd(df_snps$RDA1, na.rm = TRUE)
lim2 <- outlier_sd_threshold * sd(df_snps$RDA2, na.rm = TRUE)
df_snps$Outlier_RDA1 <- abs(df_snps$RDA1) > lim1
df_snps$Outlier_RDA2 <- abs(df_snps$RDA2) > lim2
df_snps$is_outlier <- df_snps$Outlier_RDA1 | df_snps$Outlier_RDA2
rda_outliers <- df_snps$SNP_ID[df_snps$is_outlier]

# Manuscript RDA triplot: site scores, predictor arrows, variety shapes and habitat colors.
rda_scores <- vegan::scores(rda_model, choices = c(1, 2), display = c("sites", "bp"))
df_sites <- data.frame(rda_scores$sites, SampleID = rownames(rda_scores$sites),
                       Species = meta_gen$Species, Habitat_new = meta_gen$Habitat_new)
df_env <- data.frame(rda_scores$biplot, Variable = rownames(rda_scores$biplot))
eig <- rda_model$CCA$eig
var_RDA1 <- if (length(eig) >= 2 && sum(eig) > 0) round(eig[1] / sum(eig) * 100, 1) else NA
var_RDA2 <- if (length(eig) >= 2 && sum(eig) > 0) round(eig[2] / sum(eig) * 100, 1) else NA
habitat_colors <- c(Inland = "#d95f02", Seashore = "#1b9e77", Island = "#7570b3")
species_shape <- gl_main@other$palettes$shape
stopifnot(!is.null(species_shape))
arrow_mul <- 3
p_rda <- ggplot() +
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray80") +
  geom_vline(xintercept = 0, linetype = "dashed", color = "gray80") +
  geom_point(data = df_sites, aes(RDA1, RDA2, color = Habitat_new, shape = Species),
             size = 3, alpha = 0.75) +
  geom_segment(data = df_env, aes(x = 0, y = 0, xend = RDA1 * arrow_mul,
                                  yend = RDA2 * arrow_mul),
               arrow = arrow(length = unit(0.2, "cm")), linewidth = 1) +
  geom_text_repel(data = df_env, aes(RDA1 * arrow_mul, RDA2 * arrow_mul, label = Variable),
                  fontface = "bold", size = 4, box.padding = 0.5, max.overlaps = 50) +
  scale_color_manual(values = habitat_colors, name = "Habitat") +
  scale_shape_manual(values = species_shape, name = "Type") +
  theme_classic(base_size = 14) +
  labs(title = "RDA: Genotype-environment association",
       subtitle = sprintf("Full model P = %.4f | Adjusted R² = %.2f%% | Variables = %d",
                          signif_full$`Pr(>F)`[1], R2$adj.r.squared * 100, ncol(env_core)),
       x = paste0("RDA1 (", var_RDA1, "%)"), y = paste0("RDA2 (", var_RDA2, "%)")) +
  theme(legend.position = "bottom")

# Save results: the pure SNP-ID list is read by the three-method consensus script.
write.csv(as.data.frame(signif_full), file.path(dir_out, "Table_RDA_Full_Model_Permutation_Test.csv"))
write.csv(data.frame(Selected_Variable = names(env_core)),
          file.path(dir_out, "Table_RDA_Selected_Environmental_Variables_VIF.csv"), row.names = FALSE)
write.csv(df_snps, file.path(dir_out, "Table_RDA_All_SNP_Scores.csv"), row.names = FALSE)
write.csv(df_snps[df_snps$is_outlier, , drop = FALSE],
          file.path(dir_out, "Table_RDA_Candidate_SNPs.csv"), row.names = FALSE)
write.csv(data.frame(SNP_ID = rda_outliers),
          file.path(dir_out, "Table_RDA_Candidate_SNPs_Pure_List.csv"), row.names = FALSE)
ggsave(file.path(dir_out, "Figure_RDA_Triplot.pdf"), p_rda, width = 9, height = 7)
saveRDS(list(parameters = list(vif_threshold = vif_threshold, permutations = n_perm,
                                outlier_sd_threshold = outlier_sd_threshold),
             selected_env_variables = names(env_core), excluded_env_variables_vif = vif_res@excluded,
             rda_model = rda_model, signif_full = signif_full, R2 = R2,
             loading_thresholds = c(RDA1 = lim1, RDA2 = lim2), candidate_snp_ids = rda_outliers),
        "00_rds_modules/10_RDA_candidate_SNPs.rds")
