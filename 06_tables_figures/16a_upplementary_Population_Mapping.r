# 16_1_Supplementary_Population_Mapping.R
# Supplementary table:
# Mapping between 14 regional sampling regions and 46 finer-scale populations
#
# Scientific purpose:
#   - Region_L2: 14 broader sampling regions used in regional-level analyses
#   - Pop_L3: 46 finer-scale sampling populations used in the PoMo analysis
#   - N: number of individuals in each Pop_L3 in the final metadata
#   - Coordinates: sampling coordinates of each Pop_L3

library(here)
library(readr)
library(dplyr)

# 1. Paths
pop_meta <- read.csv(
  here(
    "03_Public_Code",
    "data",
    "metadata",
    "population_metadata.csv"
  )
)

sampling_area_map <- pop_meta %>%
  select(
    Region_L2,
    Sampling_area
  ) %>%
  distinct()

dir_out <- here(
  "04_Results_Archive",
  "16_Supplementary_Population_Mapping"
)

dir.create(
  dir_out,
  recursive = TRUE,
  showWarnings = FALSE
)

# 2. Read final dataset metadata

gl_main <- readRDS(
  here(
    "00_rds_modules",
    "02_final_genlight.rds"
  )
)

stopifnot(adegenet::nInd(gl_main) == 359)
stopifnot(adegenet::nLoc(gl_main) == 8625)

meta <- gl_main@other$meta

cat("[DATA] Individuals:", nrow(meta), "\n")
cat("[DATA] Region_L2:", dplyr::n_distinct(meta$Region_L2), "\n")
cat("[DATA] Pop_L3:", dplyr::n_distinct(meta$Pop_L3), "\n")

stopifnot(nrow(meta) == 359)
stopifnot(n_distinct(meta$Region_L2) == 14)
stopifnot(n_distinct(meta$Pop_L3) == 46)

# 3. Basic checks

# Each Pop_L3 should belong to exactly one Region_L2
region_check <- meta %>%
  distinct(Pop_L3, Region_L2) %>%
  count(Pop_L3, name = "n_regions") %>%
  filter(n_regions != 1)

if (nrow(region_check) > 0) {
  stop("Some Pop_L3 populations are assigned to multiple Region_L2 regions.")
}

# Each Pop_L3 should have one sampling coordinate
coord_check <- meta %>%
  group_by(Pop_L3) %>%
  summarise(
    n_lat = n_distinct(Lat),
    n_lon = n_distinct(Lon),
    .groups = "drop"
  ) %>%
  filter(n_lat != 1 | n_lon != 1)

if (nrow(coord_check) > 0) {
  print(coord_check)
  stop("Some Pop_L3 populations contain multiple sampling coordinates.")
}

# Each Pop_L3 should contain only one variety
variety_check <- meta %>%
  distinct(Pop_L3, Species) %>%
  count(Pop_L3, name = "n_varieties") %>%
  filter(n_varieties != 1)

if (nrow(variety_check) > 0) {
  print(variety_check)
  stop("Some Pop_L3 populations contain multiple varieties.")
}

# 4. Generate supplementary mapping table

table_s_sampling <- meta %>%
  group_by(
    Species,
    Region_L2,
    Pop_L3,
    Lat,
    Lon
  ) %>%
  summarise(
    N = n(),
    .groups = "drop"
  ) %>%
  left_join(
    sampling_area_map,
    by = "Region_L2"
  ) %>%
  mutate(
    Species = factor(
      Species,
      levels = c(
        "Miscanthus sinensis var. sinensis",
        "Miscanthus sinensis var. condensatus"
      )
    )
  ) %>%
  arrange(
    Species,
    Region_L2,
    Pop_L3
  ) %>%
  transmute(
    Variety = as.character(Species),
    `Sampling area` = Sampling_area,
    `Region abbreviation` = Region_L2,
    Population = Pop_L3,
    N = N,
    `Latitude (°N)` = round(Lat, 5),
    `Longitude (°E)` = round(Lon, 5)
  )

# 5. Final checks

stopifnot(nrow(table_s_sampling) == 46)
stopifnot(sum(table_s_sampling$N) == 359)
stopifnot(!any(is.na(table_s_sampling$`Sampling area`)))

cat("\n[OUTPUT] Supplementary populations:", nrow(table_s_sampling), "\n")
cat("[OUTPUT] Total individuals:", sum(table_s_sampling$N), "\n")

print(table_s_sampling)

# 6. Export

write_csv(
  table_s_sampling,
  file.path(
    dir_out,
    "Table_S_Sampling_Information.csv"
  )
)

cat(
  "\nDone:",
  file.path(
    dir_out,
    "Table_S_Sampling_Information.csv"
  ),
  "\n"
)
