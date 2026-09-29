# 科学目的：以 Region_L2 标记 46 个 Pop_L3 的 PoMo population tree。
# 输入：PoMo treefile、publication sample_metadata.csv；输出：midpoint-rooted PDF/PNG。
# 流程：读取树与 population→region 对照 → midpoint rooting → 绘图保存。

library(ape)
library(phangorn)
library(ggtree)
library(dplyr)
library(readr)
library(ggplot2)

# Load data
dir_out <- "04_Results_Archive/03_ML_Tree/PoMo_population"
tree <- read.tree(file.path(dir_out, "population_PoMo.treefile"))
metadata <- read_csv("03_Public_Code/data/metadata/sample_metadata.csv", show_col_types=FALSE)
pop_info <- metadata %>% distinct(Pop_L3, Region_L2)
stopifnot(Ntip(tree)==46, all(tree$tip.label %in% pop_info$Pop_L3))

# Prepare data / Analysis: branch lengths are not ordinary substitutions/site.
tree <- midpoint(tree)
region_palette <- c(AOG="#5050FFFF", BOS="#CE3D32FF", FUJ="#749B58FF",
                    HJ="#466983FF", HKN="#5DB1DDFF", IBR="#802268FF",
                    ISU="#6BD76BFF", KOZ="#D595A7FF", MIY="#E6C200",
                    MUR="#89288FFF", NI="#D2AF81FF", OSH="#87C55FFF",
                    SH="#9CB084FF", TC="#808080")
p <- ggtree(tree) %<+% pop_info +
  geom_tippoint(aes(color=Region_L2), size=2) +
  geom_tiplab(aes(color=Region_L2), size=2.5, offset=0.00008,
              show.legend=FALSE) +
  scale_color_manual(values=region_palette) + theme_tree2() + labs(color="Region")

# Save results: no bootstrap labels are added.
ggsave(file.path(dir_out, "population_PoMo_tree.pdf"), p, width=10, height=12)
ggsave(file.path(dir_out, "population_PoMo_tree.png"), p, width=10, height=12, dpi=300)
