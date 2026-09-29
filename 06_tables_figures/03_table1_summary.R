# 科学目的：生成 Table 1：Stacks 全位点 Ho/He/pi/FIS 与 Genlight LD-pruned SNP 的 AR。
# 输入：Stacks All positions summary、最终 359 × 8625 Genlight、regional metadata。
# 输出：Table1_Island_Level_Final.csv；两套 SNP universes 不混用。
# 流程：读取 Stacks 全位点块 → hierfstat rarefied AR → region 注释 → 合并输出。

library(readr)
library(dplyr)
library(adegenet)
library(hierfstat)

dir_out <- "04_Results_Archive/13_Table1_260901"
dir.create(dir_out, recursive=TRUE, showWarnings=FALSE)

# Load data: skip the variant-only section; keep All positions (variant and fixed).
stats_file <- "01_Data/vcf_260901_fin/table1_final359_stats/populations.sumstats_summary.tsv"
lines <- readLines(stats_file)
split_idx <- grep("# All positions \\(variant and fixed\\)", lines)
stopifnot(length(split_idx)==1)
raw_stats <- read_tsv(stats_file, skip=split_idx, show_col_types=FALSE)
table_raw <- raw_stats %>% select(
  Pop_ID=`# Pop ID`, Private, Num_Indv,
  Obs_Het, Obs_Het_StdErr=`StdErr...15`,
  Exp_Het, Exp_Het_StdErr=`StdErr...21`,
  Pi, Pi_StdErr=`StdErr...27`, Fis, Fis_StdErr=`StdErr...30`)
gl_main <- readRDS("00_rds_modules/02_final_genlight.rds")
stopifnot(nInd(gl_main)==359, nLoc(gl_main)==8625,
          identical(as.character(gl_main@other$meta$SampleID), as.character(indNames(gl_main))))

# Prepare data: AR is computed on LD-pruned genotypes, separately from Stacks metrics.
pop_vec <- trimws(as.character(gl_main@other$meta$Region_L2))
valid <- !is.na(pop_vec) & pop_vec!=""
gl_core <- gl_main[valid, ]
pop_vec <- pop_vec[valid]
stopifnot(!any(grepl("-(TW|OIT)-", indNames(gl_core))))
sample_size_table <- data.frame(Pop_ID=names(table(pop_vec)),
                                Sample_Size=as.integer(table(pop_vec)))
region_levels <- sort(unique(pop_vec))
stopifnot(length(region_levels)==14)
gt <- as.matrix(gl_core)
stopifnot(all(is.na(gt) | gt %in% 0:2))
geno_hf <- matrix(NA_integer_, nrow(gt), ncol(gt), dimnames=dimnames(gt))
geno_hf[gt==0] <- 11L
geno_hf[gt==1] <- 12L
geno_hf[gt==2] <- 22L

# Analysis: hierfstat::allelic.richness performs the original rarefaction.
hf <- data.frame(pop=as.integer(factor(pop_vec, levels=region_levels)),
                 geno_hf, check.names=FALSE)
ar_result <- hierfstat::allelic.richness(hf)
AR_matrix <- ar_result$Ar[, as.character(seq_along(region_levels)), drop=FALSE]
ar_table <- data.frame(Pop_ID=region_levels,
                       Allelic_Richness=colMeans(AR_matrix, na.rm=TRUE))
pop_meta <- read.csv("03_Public_Code/data/metadata/population_metadata.csv")
sample_meta <- read.csv("03_Public_Code/data/metadata/sample_metadata.csv")
habitat_map <- sample_meta %>% select(Region_L2, Habitat_new) %>%
  distinct(Region_L2, .keep_all=TRUE)
# Table 1 保留现有的 Region_L2 层级分组；不以个体/Pop_L3 的 Species 替换稿件标签。
pop_info <- pop_meta %>% left_join(habitat_map, by="Region_L2", suffix=c("_pop", "")) %>%
  select(Region_L2, Sampling_area, Species=Table1_Variety, Habitat_new, Lat, Lon,
         distance_to_coast_km, Bio1, Bio12) %>%
  distinct(Region_L2, .keep_all=TRUE) %>%
  mutate(Lat=round(Lat, 2), Lon=round(Lon, 2),
         Temp_Bio1=round(Bio1, 1), Precip_Bio12=round(Bio12, 0),
         distance_to_coast_km=round(distance_to_coast_km, 2))
table_final <- table_raw %>%
  left_join(pop_info, by=c("Pop_ID"="Region_L2")) %>%
  left_join(ar_table, by="Pop_ID") %>%
  left_join(sample_size_table, by="Pop_ID") %>%
  filter(!is.na(Sampling_area)) %>%
  mutate(Species=recode(Species,
                        "Miscanthus sinensis"="Miscanthus sinensis var. sinensis",
                        "Miscanthus condensatus"="Miscanthus sinensis var. condensatus"),
         Species=factor(Species, levels=c("Miscanthus sinensis var. sinensis",
                                          "Miscanthus sinensis var. condensatus")),
         Abbreviation=Pop_ID, N=Sample_Size,
         `Allelic richness`=sprintf("%.4f", Allelic_Richness),
         Ho=sprintf("%.4f ± %.4f", Obs_Het, Obs_Het_StdErr),
         He=sprintf("%.4f ± %.4f", Exp_Het, Exp_Het_StdErr),
         `π`=sprintf("%.5f ± %.5f", Pi, Pi_StdErr),
         Fis_display=ifelse(abs(Fis)<0.0005, 0, Fis),
         Fis=sprintf("%.3f ± %.3f", Fis_display, Fis_StdErr)) %>%
  select(Variety=Species, `Sampling area`=Sampling_area, Abbreviation,
         `Latitude (°N)`=Lat, `Longitude (°E)`=Lon,
         `Annual mean temperature (°C)`=Temp_Bio1,
         `Annual precipitation (mm)`=Precip_Bio12,
         N, `Allelic richness`, Ho, He, `π`, Fis) %>%
  arrange(Variety, `Sampling area`) %>% group_by(Variety) %>%
  mutate(Variety=ifelse(row_number()==1, as.character(Variety), "")) %>% ungroup()

# Save results: CSV is the manuscript table; no Word preview or archive RDS.
write.csv(ar_table, file.path(dir_out, "Table1_Rarefied_AR_hierfstat.csv"), row.names=FALSE)
write.csv(table_final, file.path(dir_out, "Table1_Island_Level_Final.csv"), row.names=FALSE)
