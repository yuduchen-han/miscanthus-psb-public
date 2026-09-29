# 科学目的：从最终 VCF 构建供所有下游分析使用的 Genlight。
# 输入：359 个体、8625 个 LD-pruned SNP 的 VCF；publication sample_metadata.csv。
# 输出：00_rds_modules/02_final_genlight.rds。
# 流程：VCF 转换 → 清理 SampleID → 按样本显式匹配 metadata → 设置 Region_L2 与 palettes。

library(vcfR)
library(adegenet)

# Load data
vcf <- read.vcfR("01_Data/vcf_260902_8625_fin/populations_final_359_8625.vcf",
                 verbose = FALSE)
gl_final <- vcfR2genlight(vcf)
meta_df <- read.csv("03_Public_Code/data/metadata/sample_metadata.csv")

# Prepare data: PLINK 后缀清理与原 workflow 完全一致；不得再次筛选样本/SNP。
indNames(gl_final) <- sub("_.*", "", indNames(gl_final))
stopifnot(nInd(gl_final) == 359, nLoc(gl_final) == 8625,
          !anyDuplicated(indNames(gl_final)), !anyDuplicated(meta_df$SampleID))
meta_aligned <- meta_df[match(indNames(gl_final), meta_df$SampleID), , drop = FALSE]
rownames(meta_aligned) <- NULL
stopifnot(identical(as.character(meta_aligned$SampleID), as.character(indNames(gl_final))),
          !anyNA(meta_aligned$Region_L2))

# Analysis: Region_L2 population assignment and palettes used by downstream figures.
gl_final@other$meta <- meta_aligned
pop(gl_final) <- meta_aligned$Region_L2
region_palette <- c(AOG="#5050FFFF", BOS="#CE3D32FF", FUJ="#749B58FF",
                    HJ="#466983FF", HKN="#5DB1DDFF", IBR="#802268FF",
                    ISU="#6BD76BFF", KOZ="#D595A7FF", MIY="#E6C200",
                    MUR="#89288FFF", NI="#D2AF81FF", OSH="#87C55FFF",
                    SH="#9CB084FF", TC="#808080")
macro_palette <- c("Izu Islands"="#3366FF", "Boso Peninsula"="#FF3333",
                   Inland="#00A86B", "Izu Peninsula"="#9933FF",
                   "Miura Peninsula"="#FF9933", Other="#808080")
species_shape <- c("Miscanthus sinensis var. sinensis"=17,
                   "Miscanthus sinensis var. condensatus"=16)
habitat_palette <- c(Inland="#d95f02", Seashore="#1b9e77", Island="#7570b3")
gl_final@other$palettes <- list(region=region_palette, macro=macro_palette,
                                shape=species_shape, habitat=habitat_palette)

# Save results
saveRDS(gl_final, "00_rds_modules/02_final_genlight.rds")
