# 科学目的：从最终 359 × 8625 Genlight 估计 sNMF K=1:8，并展示 K=2:5 ancestry。
# 输入：02_final_genlight.rds；输出：LEA project、全部 cross-entropy、K=2:5 图与图源数据。
# 流程：转 .geno（缺失=9）→ 1000 runs/K → 最低 CE run → ancestry barplots。

library(LEA)
library(adegenet)
library(dplyr)
library(tidyr)
library(ggplot2)
library(ggsci)
library(patchwork)

dir_out <- "04_Results_Archive/04_structure"
dir.create(dir_out, recursive=TRUE, showWarnings=FALSE)
geno_file <- file.path(dir_out, "populations_final_359_8625.geno")

# Load data
gl_main <- readRDS("00_rds_modules/02_final_genlight.rds")
stopifnot(nInd(gl_main)==359, nLoc(gl_main)==8625,
          identical(as.character(gl_main@other$meta$SampleID), as.character(indNames(gl_main))))

# Prepare data: LEA .geno uses SNP rows, individual columns and 9 for missing calls.
mat_geno <- t(as.matrix(gl_main))
mat_geno[is.na(mat_geno)] <- 9
write.table(mat_geno, geno_file, sep="", row.names=FALSE, col.names=FALSE, quote=FALSE)

# Analysis: force a fresh LEA project; keep every K/run's cross-entropy.
project_ld <- snmf(geno_file, K=1:8, entropy=TRUE, repetitions=1000, project="force")
ce_all <- bind_rows(lapply(1:8, function(k) {
  ce <- as.numeric(cross.entropy(project_ld, K=k))
  data.frame(K=k, Run=seq_along(ce), CrossEntropy=ce)
}))
ce_summary <- ce_all %>% group_by(K) %>%
  summarise(Min=min(CrossEntropy, na.rm=TRUE), Mean=mean(CrossEntropy, na.rm=TRUE),
            SD=sd(CrossEntropy, na.rm=TRUE), Best_run=which.min(CrossEntropy), .groups="drop")
# K=3 has the lowest minimum CE in the established final analysis.
region_l1_order <- c("Inland", "Miura Peninsula", "Boso Peninsula",
                     "Izu Peninsula", "Izu islands")
meta <- gl_main@other$meta %>%
  mutate(Region_L1=factor(Region_L1, levels=region_l1_order), Lat=as.numeric(Lat))
stopifnot(!anyNA(meta$Region_L1))
plots <- list()
bar_data <- list()
for (k in 2:5) {
  q <- Q(project_ld, K=k, run=ce_summary$Best_run[ce_summary$K==k])
  stopifnot(nrow(q)==nInd(gl_main))
  df_q <- as.data.frame(q)
  names(df_q) <- paste0("Cluster", seq_len(k))
  df_q$SampleID <- indNames(gl_main)
  df_plot <- inner_join(df_q, meta, by="SampleID") %>%
    arrange(Region_L1, desc(Lat), Pop_L3, SampleID)
  stopifnot(nrow(df_plot)==nInd(gl_main))
  sample_order <- unique(df_plot$SampleID)
  df_plot <- df_plot %>% pivot_longer(starts_with("Cluster"),
                                      names_to="Ancestry", values_to="Proportion") %>%
    mutate(SampleID=factor(SampleID, levels=sample_order))
  bar_data[[as.character(k)]] <- mutate(df_plot, K=k)
  p <- ggplot(df_plot, aes(SampleID, Proportion, fill=Ancestry)) +
    geom_col(width=1) + scale_fill_igv() +
    facet_grid(~Region_L1, scales="free_x", space="free_x", drop=TRUE) +
    labs(y=paste0("K = ", k), x=NULL) + theme_minimal() +
    theme(legend.position="none", axis.text.x=element_blank(),
          axis.ticks.x=element_blank(), panel.grid=element_blank(),
          panel.spacing=unit(0.1, "lines"),
          axis.title.y=element_text(angle=0, vjust=0.5, face="bold", size=12))
  if (k!=2) p <- p + theme(strip.text.x=element_blank(), strip.background=element_blank())
  else p <- p + theme(strip.text.x=element_text(size=6, face="bold"))
  plots[[as.character(k)]] <- p
}

# Save results: the LEA project is also consumed by the separate K=3 map script.
write.csv(ce_all, file.path(dir_out, "CrossEntropy_All_Runs_K1_K8.csv"), row.names=FALSE)
write.csv(ce_summary, file.path(dir_out, "CrossEntropy_Summary_K1_K8.csv"), row.names=FALSE)
write.csv(bind_rows(bar_data), file.path(dir_out, "Ancestry_Barplots_K2_to_K5_Data.csv"),
          row.names=FALSE)
p_ce <- ggplot(ce_summary, aes(K, Min)) + geom_line(linewidth=0.7) +
  geom_point(size=2.8) + scale_x_continuous(breaks=1:8) +
  labs(title="Cross-Entropy for Optimal K (LD-pruned)",
       x="Number of ancestral populations", y="Minimum cross-entropy") +
  theme_classic(base_size=13)
ggsave(file.path(dir_out, "Optimal_K_CrossEntropy_LD.pdf"), p_ce, width=8, height=6)
ggsave(file.path(dir_out, "Optimal_K_CrossEntropy_LD.png"), p_ce, width=8, height=6, dpi=600)
combined <- wrap_plots(plots, ncol=1)
ggsave(file.path(dir_out, "Structure_Combined_K2_to_K5_Final.pdf"), combined, width=16, height=10)
ggsave(file.path(dir_out, "Structure_Combined_K2-K5_Final.png"), combined, width=16, height=10, dpi=600)
