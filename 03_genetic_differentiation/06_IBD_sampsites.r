# 科学目的：以 Pop_L3 为单位检验遗传距离 FST/(1-FST) 与地理距离的 IBD。
# 输入：00_rds_modules/02_final_genlight.rds（最终 359 个体、8625 SNP）。
# 输出：原始 pairwise FST、距离与配对表、Mantel/回归摘要，以及供 Figure 6 使用的 RDS。
# 流程：核对样本 → 筛选采样点 → WC FST → 经纬度质心距离 → Mantel → 保存。

library(adegenet)
library(hierfstat)
library(vegan)
library(geosphere)
library(dplyr)

min_n_per_pop <- 3
n_perm <- 9999
mantel_method <- "pearson"
random_seed <- 2026
outgroup_levels <- c("TW", "OIT", "Outgroup")
dir_out <- "04_Results_Archive/11_IBD_PopL3"
dir.create(dir_out, recursive = TRUE, showWarnings = FALSE)

# Load data: the metadata order must match the genotype order.
gl_all <- readRDS("00_rds_modules/02_final_genlight.rds")
stopifnot(nInd(gl_all) == 359, nLoc(gl_all) == 8625)
meta_all <- gl_all@other$meta
stopifnot(identical(as.character(meta_all$SampleID), as.character(indNames(gl_all))))

# Prepare data: retain the raw script's outgroup, coordinate and population filters.
keep_core <- !( (!is.na(meta_all$Region_L2) & meta_all$Region_L2 %in% outgroup_levels) |
                  grepl("-(TW|OIT)-", meta_all$SampleID) )
gl_core <- gl_all[keep_core, ]
meta_core <- meta_all[keep_core, , drop = FALSE]
meta_core$Pop_L3 <- trimws(as.character(meta_core$Pop_L3))
meta_core$Lat <- suppressWarnings(as.numeric(meta_core$Lat))
meta_core$Lon <- suppressWarnings(as.numeric(meta_core$Lon))
valid <- !is.na(meta_core$Pop_L3) & meta_core$Pop_L3 != "" &
  !is.na(meta_core$Lat) & !is.na(meta_core$Lon)
valid_populations <- names(which(table(meta_core$Pop_L3[valid]) >= min_n_per_pop))
keep <- valid & meta_core$Pop_L3 %in% valid_populations
gl <- gl_core[keep, ]
meta <- meta_core[keep, , drop = FALSE]
gl@other$meta <- meta
stopifnot(identical(as.character(meta$SampleID), as.character(indNames(gl))))
pop(gl) <- droplevels(factor(meta$Pop_L3))
stopifnot(nPop(gl) >= 3)

# Analysis: diploid genlight dosages 0/1/2 → hierfstat 11/12/22; keep missing NA.
geno <- as.matrix(gl)
stopifnot(all(is.na(geno) | geno %in% 0:2))
geno_h <- matrix(NA_integer_, nrow(geno), ncol(geno))
geno_h[geno == 0] <- 11L
geno_h[geno == 1] <- 12L
geno_h[geno == 2] <- 22L
colnames(geno_h) <- if (is.null(locNames(gl))) paste0("Locus_", seq_len(ncol(geno))) else
  make.unique(as.character(locNames(gl)))
pop_factor <- droplevels(factor(meta$Pop_L3))
population_lookup <- data.frame(Population_code = seq_along(levels(pop_factor)),
                                Population = levels(pop_factor))
fst <- as.matrix(pairwise.WCfst(data.frame(pop = as.integer(pop_factor), geno_h,
                                           check.names = FALSE), diploid = TRUE))
labels <- population_lookup$Population[match(as.integer(rownames(fst)),
                                             population_lookup$Population_code)]
stopifnot(!anyNA(labels))
dimnames(fst) <- list(labels, labels)
population_order <- sort(labels)
fst <- fst[population_order, population_order, drop = FALSE]
stopifnot(isTRUE(all.equal(fst, t(fst), check.attributes = FALSE)))
diag(fst) <- 0
fst_distance <- fst
fst_distance[fst_distance < 0] <- 0 # 原始 FST 保留负值；仅距离计算截为零。
stopifnot(!any(fst_distance >= 1, na.rm = TRUE))
genetic_distance_matrix <- fst_distance / (1 - fst_distance)
diag(genetic_distance_matrix) <- 0

population_coordinates <- meta %>%
  group_by(Population = Pop_L3) %>%
  summarise(Longitude = mean(Lon, na.rm = TRUE), Latitude = mean(Lat, na.rm = TRUE),
            Longitude_min = min(Lon, na.rm = TRUE), Longitude_max = max(Lon, na.rm = TRUE),
            Latitude_min = min(Lat, na.rm = TRUE), Latitude_max = max(Lat, na.rm = TRUE),
            N = n(), .groups = "drop")
if (any((population_coordinates$Longitude_max - population_coordinates$Longitude_min) > 0.01 |
        (population_coordinates$Latitude_max - population_coordinates$Latitude_min) > 0.01))
  warning("A population spans more than 0.01 degrees; using its centroid.")
population_coordinates <- population_coordinates[match(population_order,
                                                         population_coordinates$Population), ]
stopifnot(identical(population_coordinates$Population, population_order))
coordinates <- as.matrix(population_coordinates[, c("Longitude", "Latitude")])
geographic_distance_matrix <- distm(coordinates, fun = distHaversine) / 1000
dimnames(geographic_distance_matrix) <- list(population_order, population_order)
diag(geographic_distance_matrix) <- 0

pairs <- which(lower.tri(fst), arr.ind = TRUE)
pairwise_table <- data.frame(Population_1 = population_order[pairs[, "row"]],
                             Population_2 = population_order[pairs[, "col"]],
                             FST_raw = fst[pairs], FST_for_distance = fst_distance[pairs],
                             Genetic_distance = genetic_distance_matrix[pairs],
                             Geographic_distance_km = geographic_distance_matrix[pairs]) %>%
  mutate(Population_pair = paste(Population_1, Population_2, sep = " – "),
         Log_geographic_distance = log10(Geographic_distance_km)) %>%
  filter(is.finite(Genetic_distance), is.finite(Geographic_distance_km),
         Geographic_distance_km > 0)
habitat_lookup <- meta %>% distinct(Pop_L3, Habitat_new)
pairwise_table <- pairwise_table %>%
  left_join(habitat_lookup, by = c("Population_1" = "Pop_L3")) %>%
  rename(Habitat_1 = Habitat_new) %>%
  left_join(habitat_lookup, by = c("Population_2" = "Pop_L3")) %>%
  rename(Habitat_2 = Habitat_new) %>%
  mutate(Pair_type = case_when(
    Habitat_1 == "Island" & Habitat_2 == "Island" ~ "Within islands",
    xor(Habitat_1 == "Island", Habitat_2 == "Island") ~ "Mainland vs. Islands",
    TRUE ~ "Within mainland"))

set.seed(random_seed)
mantel_result <- mantel(as.dist(genetic_distance_matrix), as.dist(geographic_distance_matrix),
                        method = mantel_method, permutations = n_perm, na.rm = TRUE)
lm_ibd <- lm(Genetic_distance ~ Geographic_distance_km, data = pairwise_table)
lm_summary <- summary(lm_ibd)
statistical_summary <- data.frame(
  Analysis = c("Mantel test", "Linear regression"),
  Statistic = c(unname(mantel_result$statistic), unname(coef(lm_ibd)["Geographic_distance_km"])),
  P_value = c(mantel_result$signif,
              coef(lm_summary)["Geographic_distance_km", "Pr(>|t|)"]),
  R_squared = c(NA_real_, lm_summary$r.squared),
  Permutations = c(n_perm, NA_integer_))

# Save results: the RDS fields are read by 06b_IBD_plot_din.r (Figure 6).
write.csv(fst, file.path(dir_out, "Pairwise_FST_Weir_Cockerham.csv"))
write.csv(genetic_distance_matrix, file.path(dir_out, "Genetic_Distance_FST_over_1_minus_FST.csv"))
write.csv(geographic_distance_matrix, file.path(dir_out, "Geographic_Distance_km.csv"))
write.csv(pairwise_table, file.path(dir_out, "IBD_Pairwise_Data_PopL3.csv"), row.names = FALSE)
write.csv(statistical_summary, file.path(dir_out, "IBD_Statistical_Summary.csv"), row.names = FALSE)
saveRDS(list(pairwise_table = pairwise_table, mantel_result = mantel_result),
        "00_rds_modules/11_ibd_popl3_results.rds")
