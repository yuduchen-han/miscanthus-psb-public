# 科学目的：展示最终 359 × 8625 数据的全基因组 PCA。
# 输入：00_rds_modules/02_final_genlight.rds；输出：PCA 模型、Figure PDF/PNG 与图源 CSV。
# 流程：核对样本/metadata → Habitat/variety 因子 → glPca → 轴方差 → 绘图保存。

library(adegenet)
library(ggplot2)

# Load data
gl_main <- readRDS("00_rds_modules/02_final_genlight.rds")
stopifnot(nInd(gl_main) == 359, nLoc(gl_main) == 8625,
          identical(as.character(gl_main@other$meta$SampleID), as.character(indNames(gl_main))))

# Prepare data
gl_main@other$meta$Habitat_new <- factor(gl_main@other$meta$Habitat_new,
                                            levels = c("Inland", "Seashore", "Island"))
gl_main@other$meta$Species <- factor(gl_main@other$meta$Species,
  levels = c("Miscanthus sinensis var. sinensis", "Miscanthus sinensis var. condensatus"))
species_shape <- gl_main@other$palettes$shape
habitat_colors <- c(Inland="#d95f02", Seashore="#1b9e77", Island="#7570b3")

# Analysis
pca_all <- glPca(gl_main, nf = 3, parallel = TRUE, n.cores = 4)
eig_global <- round(100 * pca_all$eig / sum(pca_all$eig), 1)
df_pca_global <- cbind(as.data.frame(pca_all$scores), gl_main@other$meta)
variety_labels <- c(
  "Miscanthus sinensis var. sinensis" =
    expression(italic("Miscanthus sinensis") ~ var. ~ italic("sinensis")),
  "Miscanthus sinensis var. condensatus" =
    expression(italic("Miscanthus sinensis") ~ var. ~ italic("condensatus")))
p_pca_all <- ggplot(df_pca_global, aes(PC1, PC2)) +
  geom_hline(yintercept=0, linetype="dashed", color="gray80") +
  geom_vline(xintercept=0, linetype="dashed", color="gray80") +
  geom_point(aes(color=Habitat_new, shape=Species), size=3.5, alpha=0.8) +
  scale_color_manual(values=habitat_colors, name="Habitat") +
  scale_shape_manual(values=species_shape, labels=variety_labels, name="Variety") +
  labs(x=paste0("PC1 (", eig_global[1], "%)"),
       y=paste0("PC2 (", eig_global[2], "%)")) +
  theme_bw(base_size=14)

# Save results
saveRDS(pca_all, "00_rds_modules/02_pca_all.rds")
dir.create("04_Results_Archive/inter_results", recursive=TRUE, showWarnings=FALSE)
write.csv(df_pca_global, "04_Results_Archive/inter_results/df_pca_global.csv", row.names=FALSE)
dir.create("04_Results_Archive/01_PCA", recursive=TRUE, showWarnings=FALSE)
ggsave("04_Results_Archive/01_PCA/pca_all_habitat.pdf", p_pca_all, width=10, height=8)
ggsave("04_Results_Archive/01_PCA/pca_all_habitat.png", p_pca_all, width=10, height=8, dpi=600)
