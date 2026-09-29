# 科学目的：绘制最终 359 个体 midpoint-rooted ML tree（仅标记双重显著节点）。
# 输入：03_ml_tree_models.rds、02_final_genlight.rds；输出：manuscript PDF/PNG。
# 流程：匹配树端和 metadata → 解析 SH-aLRT/UFBoot → 按 80/95 阈值绘图。

library(adegenet)
library(ggtree)
library(ggplot2)
library(dplyr)

# Load data
dir_out <- "04_Results_Archive/03_ML_Tree"
tree_midpoint <- readRDS("00_rds_modules/03_ml_tree_models.rds")$tree_midpoint
gl_main <- readRDS("00_rds_modules/02_final_genlight.rds")
meta_df <- gl_main@other$meta
stopifnot(length(tree_midpoint$tip.label)==359,
          setequal(tree_midpoint$tip.label, meta_df$SampleID))

# Prepare data
meta_tree <- meta_df %>%
  filter(SampleID %in% tree_midpoint$tip.label) %>%
  mutate(SampleID=factor(SampleID, levels=tree_midpoint$tip.label)) %>%
  arrange(SampleID) %>% mutate(SampleID=as.character(SampleID))
region_palette <- gl_main@other$palettes$region
species_shape <- gl_main@other$palettes$shape

# Analysis: node labels contain SH-aLRT/UFBoot; mark only >=80 AND >=95.
p <- ggtree(tree_midpoint, layout="fan", open.angle=15, linewidth=0.3) %<+% meta_tree
plot_data <- p$data
plot_data$SH_aLRT <- suppressWarnings(as.numeric(sub("/.*$", "", plot_data$label)))
plot_data$UFboot <- suppressWarnings(as.numeric(sub("^.*/", "", plot_data$label)))
plot_data$is_supported <- !plot_data$isTip & !is.na(plot_data$SH_aLRT) &
  !is.na(plot_data$UFboot) & plot_data$SH_aLRT >= 80 & plot_data$UFboot >= 95
p$data <- plot_data
p_final <- p +
  geom_tippoint(aes(color=Region_L2, shape=Species), size=2.2, alpha=0.9) +
  geom_nodepoint(aes(subset=is_supported), shape=23, fill="gold", color="black",
                 size=1.5, stroke=0.2) +
  scale_color_manual(values=region_palette, name="Region", na.translate=FALSE) +
  scale_shape_manual(values=species_shape, name="Variety",
    labels=c(expression(italic("Miscanthus sinensis") ~ "var." ~ italic("sinensis")),
             expression(italic("Miscanthus sinensis") ~ "var." ~ italic("condensatus"))),
    na.translate=FALSE) +
  theme_void(base_size=14) +
  theme(legend.position="right", legend.title=element_text(size=11),
        legend.text=element_text(size=9),
        plot.margin=margin(10, 60, 10, 10, unit="pt"),
        plot.title=element_text(hjust=0.5, face="bold"),
        plot.subtitle=element_text(hjust=0.5))

# Save results
ggsave(file.path(dir_out, "Tree_01_Global_TrueDist_Final.pdf"), p_final, width=13, height=12)
ggsave(file.path(dir_out, "Tree_01_Global_TrueDist_Final.png"), p_final,
       width=13, height=12, dpi=600)
