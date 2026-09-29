# 科学目的：将固定 21 consensus SNP 的 PCA 结果绘成最终 Figure 5。
# 输入：08c_RF_Variety_21SNP.rds；输出：最终 PDF/PNG。
# 流程：读取已计算的坐标/方差 → Habitat 着色、variety 形状 → 保存。

library(ggplot2)

# Load data
res <- readRDS("00_rds_modules/08c_RF_Variety_21SNP.rds")
pca_data <- res$diagnostic_pca$plot_data
var_exp <- res$diagnostic_pca$variance_explained
stopifnot(nrow(pca_data)==359, length(var_exp)>=2)

# Prepare data / Analysis
habitat_colors <- c(Inland="#d95f02", Seashore="#1b9e77", Island="#7570b3")
variety_labels <- c(
  "Miscanthus sinensis var. sinensis" =
    expression(italic("Miscanthus sinensis") ~ var. ~ italic("sinensis")),
  "Miscanthus sinensis var. condensatus" =
    expression(italic("Miscanthus sinensis") ~ var. ~ italic("condensatus")))
p_pca_final <- ggplot(pca_data, aes(PC1, PC2, color=Habitat, shape=Variety)) +
  geom_hline(yintercept=0, linetype="dashed", color="gray80") +
  geom_vline(xintercept=0, linetype="dashed", color="gray80") +
  geom_point(size=3, alpha=0.8) +
  scale_color_manual(values=habitat_colors, breaks=c("Inland","Seashore","Island"),
                     name="Habitat") +
  scale_shape_manual(values=c(17,16), labels=variety_labels, name="Variety") +
  labs(x=paste0("PC1 (", round(var_exp[1],1), "%)"),
       y=paste0("PC2 (", round(var_exp[2],1), "%)")) +
  theme_classic(base_size=13) +
  theme(axis.line=element_blank(),
        panel.border=element_rect(colour="black", fill=NA, linewidth=0.5),
        legend.position="bottom", legend.box="vertical") +
  guides(shape=guide_legend(order=1), color=guide_legend(order=2))

# Save results
dir_out <- "04_Results_Archive/08c_RF_Variety_21SNP"
ggsave(file.path(dir_out, "Figure_RF_Variety_21SNP_Diagnostic_PCA_final.pdf"),
       p_pca_final, width=8, height=6)
ggsave(file.path(dir_out, "Figure_RF_Variety_21SNP_Diagnostic_PCA_final.png"),
       p_pca_final, width=8, height=6, dpi=600)
