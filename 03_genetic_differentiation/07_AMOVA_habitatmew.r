# 科学目的：以 Habitat_new / Pop_L3 为层级估算 AMOVA 并生成 manuscript Table 2。
# 输入：最终 359 × 8625 Genlight；输出：原始统计对象及 Table_AMOVA_Habitat_Final.csv。
# 流程：按 metadata 对齐 → 筛选有效层级 → poppr.amova → 999 次置换 → 原始方差百分比。

library(poppr)
library(ade4)
library(dplyr)

dir_out <- "04_Results_Archive/07_AMOVA_Habitat"
dir.create(dir_out, recursive=TRUE, showWarnings=FALSE)

# Load data
gl_main <- readRDS("00_rds_modules/02_final_genlight.rds")
stopifnot(nInd(gl_main)==359, nLoc(gl_main)==8625)
meta <- gl_main@other$meta
stopifnot(identical(as.character(meta$SampleID), as.character(indNames(gl_main))))

# Prepare data: retain raw filtering, including explicit TW/OIT identification.
keep_core <- !grepl("-(TW|OIT)-", meta$SampleID)
stopifnot(all(keep_core)) # final 359 must not silently include outgroups
meta_core <- meta[keep_core, , drop=FALSE]
meta_core$Habitat <- trimws(as.character(meta_core$Habitat_new))
keep_hab <- !is.na(meta_core$Habitat) & meta_core$Habitat != "" &
  !is.na(meta_core$Pop_L3) & meta_core$Pop_L3 != ""
gl_core <- gl_main[keep_core, ]
gl_core@other$meta <- meta_core
gl_hab <- gl_core[keep_hab, ]
meta_hab <- meta_core[keep_hab, , drop=FALSE]
gl_hab@other$meta <- meta_hab
strata(gl_hab) <- meta_hab[, c("Habitat", "Pop_L3"), drop=FALSE]

# Analysis: negative Sigma components are retained in the total and percentages.
amova_hab <- poppr.amova(gl_hab, ~ Habitat / Pop_L3,
                          clonecorrect=FALSE, within=FALSE)
set.seed(2026)
amova_sig <- randtest(amova_hab, nrepet=999)
var_comp_final <- amova_hab$componentsofcovariance
var_comp_final$Original_Level <- rownames(var_comp_final)
rownames(var_comp_final) <- NULL
var_comp_final <- var_comp_final %>%
  filter(!grepl("Total", Original_Level, ignore.case=TRUE)) %>%
  mutate(Source_of_variation = case_when(
    grepl("Between samples", Original_Level, ignore.case=TRUE) ~ "Among populations within habitats",
    grepl("Between Habitat", Original_Level, ignore.case=TRUE) ~ "Among habitats",
    grepl("Within samples", Original_Level, ignore.case=TRUE) ~ "Within populations",
    TRUE ~ Original_Level))
total_sigma <- sum(var_comp_final$Sigma, na.rm=TRUE)
var_comp_final <- var_comp_final %>% mutate(Percent_variation=Sigma/total_sigma*100)
sig_table <- data.frame(Test=amova_sig$names, P_value=as.numeric(amova_sig$pvalue))
get_p <- function(pattern) {
  x <- sig_table$P_value[grepl(pattern, sig_table$Test, ignore.case=TRUE)]
  if (length(x)==1) x else NA_real_
}
table_amova_final <- var_comp_final %>%
  mutate(P=case_when(
    Source_of_variation=="Among habitats" ~ get_p("between Habitat"),
    Source_of_variation=="Among populations within habitats" ~ get_p("between samples"),
    Source_of_variation=="Within populations" ~ get_p("within samples"),
    TRUE ~ NA_real_)) %>%
  transmute(`Source of variation`=Source_of_variation,
            `Variance component`=round(Sigma, 5),
            `Variation (%)`=round(Percent_variation, 2), P=P)

# Save results
write.csv(table_amova_final, file.path(dir_out, "Table_AMOVA_Habitat_Final.csv"),
          row.names=FALSE, na="")
saveRDS(list(amova=amova_hab, permutation=amova_sig,
             variance_components=var_comp_final, final_table=table_amova_final),
        "00_rds_modules/07_amova_habitat_final.rds")
