# 科学目的：从已完成的 sNMF project 绘制 K=3 的 Region_L2 × Habitat_new 地图。
# 输入：最终 Genlight、K=1:8 的 LEA project；输出：map data CSV 与最终 PDF/PNG。
# 流程：最低 cross-entropy run → Q matrix → ancestry/坐标汇总 → 手工版面偏移 → 绘图。

library(adegenet)
library(LEA)
library(ggplot2)
library(dplyr)
library(scatterpie)
library(sf)
library(rnaturalearth)
library(ggspatial)
library(ggsci)

dir_out <- "04_Results_Archive/04_structure"

# Load data
gl_main <- readRDS("00_rds_modules/02_final_genlight.rds")
stopifnot(nInd(gl_main)==359, nLoc(gl_main)==8625,
          identical(as.character(gl_main@other$meta$SampleID), as.character(indNames(gl_main))))
project_ld <- load.snmfProject(file.path(dir_out, "populations_final_359_8625.snmfProject"))
ce_k3 <- as.numeric(cross.entropy(project_ld, K=3))
stopifnot(length(ce_k3)==1000)
best_run_k3 <- which.min(ce_k3)
q_mat_3 <- Q(project_ld, K=3, run=best_run_k3)
stopifnot(nrow(q_mat_3)==359)

# Prepare data: pie offsets are display-only; Lon/Lat remain true sampling centroids.
df_q <- as.data.frame(q_mat_3)
names(df_q) <- paste0("Cluster", 1:3)
df_q$SampleID <- indNames(gl_main)
meta <- gl_main@other$meta
habitat_count <- meta %>% distinct(Region_L2, Habitat_new) %>%
  count(Region_L2, name="n_habitat")
map_pie_summary <- left_join(df_q, meta, by="SampleID") %>%
  group_by(Region_L2, Habitat_new) %>%
  summarise(across(starts_with("Cluster"), ~mean(.x, na.rm=TRUE)),
            Lat=mean(as.numeric(Lat), na.rm=TRUE),
            Lon=mean(as.numeric(Lon), na.rm=TRUE), .groups="drop") %>%
  left_join(habitat_count, by="Region_L2") %>%
  mutate(MapLabel=if_else(n_habitat>1, paste(Region_L2, Habitat_new, sep="_"),
                          Region_L2)) %>%
  filter(!is.na(Lat), !is.na(Lon))

pie_radius <- 0.20
label_gap <- 0.05
offset_table <- tibble::tribble(
  ~MapLabel, ~dx, ~dy, ~label_side, ~label_dx, ~label_dy,
  "TC", -0.45, 0.15, "right", 0.00, 0.00,
  "IBR_Inland", -0.65, -0.10, "right", 0.00, 0.00,
  "IBR_Seashore", 0.45, 0.30, "right", 0.00, 0.00,
  "FUJ", -0.45, 0.25, "left", 0.00, 0.00,
  "HKN", -0.45, 0.15, "left", 0.00, 0.00,
  "MUR", -0.25, 0.30, "left", 0.00, 0.00,
  "BOS_Inland", 0.10, 0.40, "right", 0.00, 0.00,
  "BOS_Seashore", 0.55, 0.15, "right", 0.00, 0.00,
  "ISU_Inland", -0.55, -0.20, "left", 0.00, 0.00,
  "ISU_Seashore", -1.05, -0.35, "left", 0.00, 0.00,
  "OSH", 0.45, 0.20, "right", 0.00, 0.05,
  "KOZ", -0.70, -0.25, "left", 0.00, -0.02,
  "NI", 0.65, -0.10, "right", 0.00, -0.05,
  "SH", -0.50, 0.15, "left", 0.00, 0.05,
  "MIY", 0.50, -0.30, "right", 0.00, 0.00,
  "HJ", -0.50, -0.10, "right", 0.00, 0.00,
  "AOG", -0.50, -0.20, "right", 0.00, 0.00)
map_pie_core <- map_pie_summary %>% left_join(offset_table, by="MapLabel") %>%
  mutate(dx=coalesce(dx, 0), dy=coalesce(dy, 0),
         label_dx=coalesce(label_dx, 0), label_dy=coalesce(label_dy, 0),
         PieRadius=pie_radius, plot_lon=Lon+dx, plot_lat=Lat+dy,
         label_lon=case_when(label_side=="left" ~ plot_lon-PieRadius-label_gap,
                             label_side=="right" ~ plot_lon+PieRadius+label_gap,
                             TRUE ~ plot_lon)+label_dx,
         label_lat=plot_lat+label_dy,
         label_hjust=case_when(label_side=="left" ~ 1,
                               label_side=="right" ~ 0, TRUE ~ 0.5))

# Analysis: display original locations with connecting segments to shifted pies.
japan_sf <- if (requireNamespace("rnaturalearthhires", quietly=TRUE))
  ne_countries(scale="large", country="japan", returnclass="sf") else
  ne_countries(scale="medium", country="japan", returnclass="sf")
p_map_k3 <- ggplot() +
  geom_sf(data=japan_sf, fill="#e0e4e8", color="#bdc3c7", linewidth=0.2) +
  geom_segment(data=map_pie_core, aes(x=Lon, y=Lat, xend=plot_lon, yend=plot_lat),
               linewidth=0.5, colour="grey30") +
  geom_point(data=map_pie_core, aes(Lon, Lat), size=1.5, colour="grey20") +
  geom_scatterpie(data=map_pie_core,
                  aes(x=plot_lon, y=plot_lat, group=MapLabel, r=PieRadius),
                  cols=paste0("Cluster", 1:3), color="white", linewidth=0.3, alpha=0.95) +
  geom_label(data=map_pie_core,
             aes(x=label_lon, y=label_lat, label=MapLabel, hjust=label_hjust),
             size=3.2, fontface="bold", linewidth=0.2,
             label.padding=unit(0.15, "lines")) +
  scale_fill_igv(name="Ancestry (K = 3)") +
  annotation_north_arrow(location="tr", which_north="true", height=unit(0.7,"cm"),
                          width=unit(0.7,"cm"), pad_x=unit(0.25,"cm"),
                          pad_y=unit(0.25,"cm"),
                          style=north_arrow_orienteering(text_size=8, line_width=0.5)) +
  coord_sf(xlim=c(136.0,142.5), ylim=c(31.8,37.3), expand=FALSE) +
  theme_minimal(base_size=14) +
  theme(panel.background=element_rect(fill="#e3f0f7", color=NA),
        panel.grid.major=element_line(color="white", linewidth=0.3),
        axis.title=element_blank(), legend.position="right")

# Save results
write.csv(map_pie_core, file.path(dir_out, "Table_Structure_Map_K3_Data.csv"), row.names=FALSE)
ggsave(file.path(dir_out, "Structure_Map_K3_Core.pdf"), p_map_k3, width=10, height=8)
ggsave(file.path(dir_out, "Structure_Map_K3_Core.png"), p_map_k3, width=10, height=8, dpi=600)
