# 科学目的：用已计算的 IBD 配对距离与 Mantel 结果绘制 manuscript Figure 6。
# 输入：11_ibd_popl3_results.rds；输出：IBD_PopL3_Mantel_final.pdf/png。
# 流程：读取 pairwise_table / mantel_result → 各类群体对标色 → 保存。

library(ggplot2)

# Load data
res <- readRDS("00_rds_modules/11_ibd_popl3_results.rds")
pairwise_table <- res$pairwise_table
mantel_result <- res$mantel_result
mantel_r <- unname(mantel_result$statistic)
mantel_p <- mantel_result$signif

# Prepare data
mantel_annotation <- if (mantel_p < 0.001) {
  sprintf("Mantel~italic(r) == %.3f*','~~italic(P) < 0.001", mantel_r)
} else {
  sprintf("Mantel~italic(r) == %.3f*','~~italic(P) == %.3f", mantel_r, mantel_p)
}

# Analysis: Figure 6 uses the original all-pair lm smoothing, not a second Mantel.
p <- ggplot(pairwise_table, aes(Geographic_distance_km, Genetic_distance)) +
  geom_point(aes(color=Pair_type), size=2.5, alpha=0.75) +
  geom_smooth(method="lm", formula=y~x, se=TRUE, color="darkgrey",
              fill="lightgrey", linewidth=0.9) +
  scale_color_manual(values=c("Within mainland"="#4C78A8",
                              "Mainland vs. Islands"="#F58518",
                              "Within islands"="#E45756"), name="Comparison type") +
  annotate("text", x=Inf, y=Inf, label=mantel_annotation, parse=TRUE,
           hjust=1.05, vjust=1.2, size=4) +
  labs(x="Geographical distance (km)", y=expression(F[ST]/(1-F[ST]))) +
  theme_classic(base_size=13) +
  theme(axis.title=element_text(face="bold"), axis.line=element_blank(),
        panel.border=element_rect(colour="black", fill=NA, linewidth=0.5),
        legend.position=c(0.82,0.14), legend.background=element_blank(),
        plot.margin=margin(10,15,10,10))

# Save results
dir_out <- "04_Results_Archive/11_IBD_PopL3"
ggsave(file.path(dir_out, "IBD_PopL3_Mantel_final.pdf"), p, width=6.5, height=6)
ggsave(file.path(dir_out, "IBD_PopL3_Mantel_final.png"), p, width=6.5, height=6, dpi=600)
