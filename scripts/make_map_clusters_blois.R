# ==============================================================================
# Carte des clusters sur tout Blois — k-means sur la grille forêt complète
# avec relabel par LAI croissant (P1 sparse → P4 dense) et palette harmonisée.
#
# Outputs : outputs/figs_MEB2026_final/fig_map_clusters_blois.{png,pdf}
#           + raster GeoTIFF (out_files/clusters/blois_clusters_P1P4.tif)
# ==============================================================================

suppressMessages({
  library(data.table); library(here); library(cli); library(terra)
  library(tidyverse); library(sf)
})

source(here::here("R/config.R"))
source(here::here("R/io.R"))
source(here::here("R/forest.R"))
source(here::here("R/cluster_relabel.R"))

OUT <- here::here("outputs/figs_MEB2026_final")
dir.create(OUT, recursive = TRUE, showWarnings = FALSE)

cli_h1("Carte clusters Blois")

# 1. Load LiDAR rasters with V2 correction (PAD → LAI ×2)
cli_alert("Loading LiDAR rasters (V2, agg 20 m)...")
rasters <- load_lidar_rasters(CFG$in_dir, CFG$agg_factor)
forest  <- build_forest_dataframe(rasters$stack)
df_forest_all <- forest$df
cli_alert("Forest pixels : {nrow(df_forest_all)}")

# 2. K-means k=4 on LAI / Hmax / fCover with fixed seed (same as cLHS pipeline)
set.seed(42)
df_in <- df_forest_all %>% dplyr::select(LAI, Hmax, fCover)
ok    <- complete.cases(df_in)
df_scaled <- scale(df_in[ok, ])
km <- kmeans(df_scaled, centers = 4, iter.max = 100, nstart = 25)
df_forest_all$Cluster <- NA_integer_
df_forest_all$Cluster[ok] <- km$cluster
cli_alert("k-means done : {sum(ok)} pixels clustered into 4 groups")

# 3. Order clusters by ASCENDING mean LAI → P1 = sparsest, P4 = densest
cl_stats <- df_forest_all %>%
  filter(!is.na(Cluster)) %>%
  group_by(Cluster) %>%
  summarise(LAI_mean = mean(LAI, na.rm = TRUE),
            Hmax_mean = mean(Hmax, na.rm = TRUE),
            fCover_mean = mean(fCover, na.rm = TRUE),
            n = n(), .groups = "drop") %>%
  arrange(LAI_mean) %>%
  mutate(Profile = paste0("P", row_number()))
cli_alert("Cluster → Profile mapping (LAI ascending) :")
print(cl_stats)

map_old_to_P <- setNames(cl_stats$Profile, as.character(cl_stats$Cluster))
df_forest_all$Profile <- factor(map_old_to_P[as.character(df_forest_all$Cluster)],
                                 levels = paste0("P", 1:4))

# 4. Save GeoTIFF (Profile encoded as integer 1..4 in LAI ascending order)
out_dir_rast <- here::here("out_files/clusters")
dir.create(out_dir_rast, recursive = TRUE, showWarnings = FALSE)

profile_idx <- as.integer(df_forest_all$Profile)
template <- rasters$stack[[1]]
r_clu <- terra::rast(template); values(r_clu) <- NA_integer_
cells <- terra::cellFromXY(r_clu, as.matrix(df_forest_all[, c("x","y")]))
values(r_clu)[cells] <- profile_idx
names(r_clu) <- "Profile"
terra::writeRaster(r_clu, file.path(out_dir_rast, "blois_clusters_P1P4.tif"),
                    overwrite = TRUE, datatype = "INT1U", NAflag = 255)
cli_alert_success("Saved out_files/clusters/blois_clusters_P1P4.tif")

# 5. HOBO sensor locations (overlay)
hobo_pts <- sf::st_read(CFG$hobo_geojson, quiet = TRUE) %>%
  filter(!id_plot %in% CFG$ids_to_remove)
hobo_xy <- as.data.frame(sf::st_coordinates(hobo_pts))
names(hobo_xy) <- c("x", "y")
hobo_xy$layer <- "HOBO sensors (n=53)"

# 5b. cLHS sample locations (overlay)
clhs_xy <- as.data.frame(
  readRDS(here::here("out_files/Sensitivity_Analysis/clhs_sample_floor05_v2.rds"))
)[, c("x", "y")]
clhs_xy$layer <- "cLHS plots (n=400)"
cli_alert("HOBO points : {nrow(hobo_xy)}, cLHS plots : {nrow(clhs_xy)}")

points_all <- rbind(clhs_xy, hobo_xy)
points_all$layer <- factor(points_all$layer,
                            levels = c("cLHS plots (n=400)", "HOBO sensors (n=53)"))

# 6. Build map ggplot
cli_h2("Building map figure")

# Pixel counts and proportion per profile
prop_tab <- df_forest_all %>%
  filter(!is.na(Profile)) %>%
  count(Profile) %>%
  mutate(prop = 100 * n / sum(n))
print(prop_tab)
fwrite(prop_tab, file.path(OUT, "tab_blois_clusters_proportions.csv"))

p_map <- ggplot(filter(df_forest_all, !is.na(Profile)),
                 aes(x = x, y = y, fill = Profile)) +
  geom_raster() +
  # cLHS plots : small black dots (n=400)
  geom_point(data = points_all[points_all$layer == "cLHS plots (n=400)", ],
             aes(x = x, y = y, shape = layer, size = layer,
                  colour = layer, stroke = layer),
             inherit.aes = FALSE) +
  # HOBO sensors : larger white circles with black outline (n=53)
  geom_point(data = points_all[points_all$layer == "HOBO sensors (n=53)", ],
             aes(x = x, y = y, shape = layer, size = layer,
                  colour = layer, stroke = layer),
             inherit.aes = FALSE, fill = "white") +
  scale_fill_manual(values = PAL_CLUSTER, drop = FALSE,
                     name = "Profile",
                     guide = guide_legend(order = 1, override.aes = list(size = 8))) +
  scale_shape_manual(values = c("cLHS plots (n=400)" = 16,
                                  "HOBO sensors (n=53)" = 21),
                       name = "Plots",
                       guide = guide_legend(order = 2,
                                              override.aes = list(
                                                size = c(5, 5),
                                                stroke = c(0, 1.2),
                                                colour = c("black", "black"),
                                                fill = c("black", "white")))) +
  scale_size_manual(values = c("cLHS plots (n=400)" = 3.5,
                                 "HOBO sensors (n=53)" = 4),
                      guide = "none") +
  scale_colour_manual(values = c("cLHS plots (n=400)" = "black",
                                   "HOBO sensors (n=53)" = "black"),
                        guide = "none") +
  scale_discrete_manual("stroke",
                          values = c("cLHS plots (n=400)" = 0,
                                      "HOBO sensors (n=53)" = 1.2),
                          guide = "none") +
  coord_sf(crs = sf::st_crs(hobo_pts), datum = sf::st_crs(hobo_pts),
            expand = FALSE) +
  labs(x = "Easting (m, UTM 31N)", y = "Northing (m, UTM 31N)") +
  theme_bw(base_size = 30) +
  theme(panel.grid       = element_blank(),
        axis.text.x  = element_text(size = 20, angle = 35, hjust = 1),
        axis.text.y  = element_text(size = 20),
        axis.title.x = element_text(size = 26, face = "bold",
                                      margin = margin(t = 12)),
        axis.title.y = element_text(size = 26, face = "bold",
                                      margin = margin(r = 12)),
        axis.ticks = element_line(linewidth = 0.5),
        legend.text  = element_text(size = 26),
        legend.title = element_text(size = 30, face = "bold"),
        legend.key.size = unit(1.8, "cm"),
        legend.spacing.y = unit(0.6, "cm"),
        legend.position = "right",
        legend.box = "vertical")

# Slide 16:9 — large dpi
ggsave(file.path(OUT, "fig_map_clusters_blois.png"),
        p_map, width = 16, height = 9, dpi = 600, bg = "white")
ggsave(file.path(OUT, "fig_map_clusters_blois.pdf"),
        p_map, width = 16, height = 9, device = cairo_pdf)
cli_alert_success("Saved fig_map_clusters_blois.{{png,pdf}} (16:9 @ 600 dpi)")

cli_h1("Done.")
