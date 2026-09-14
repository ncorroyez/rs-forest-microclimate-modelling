# ==============================================================================
# A9 — "Illusion optique" Sentinel-2 (pont vers le Chapitre 2).
# H3 : le LAI Sentinel-2 sature en couvert dense (perd la structure intra-canopée)
# et son lien avec le buffering n'est qu'un PROXY de la structure LiDAR — il
# disparaît une fois le LAI LiDAR contrôlé. Donnée réelle, 53 placettes.
#
#   (A) S2 LAI vs LiDAR LAI  -> saturation optique (S2 plafonne quand LiDAR monte)
#   (B) buffering ΔTmax_obs vs S2 LAI : corrélation brute MAIS partielle | LiDAR ≈ 0
#
# Sortie : outputs/figures_pipeline/annex/fig_annex_s2_optical_illusion.{png,pdf}
# ==============================================================================
suppressMessages({ library(here); library(data.table); library(terra); library(tidyverse); library(patchwork) })
source(here::here("scripts/_article_style.R"))
OUT <- here::here("outputs/figures_pipeline/annex"); dir.create(OUT, recursive = TRUE, showWarnings = FALSE)
OAK <- path.expand("~/Documents/NC_Full/03_RESULTS/Blois/Metrics/Oak_Only")

# buffering observé par placette
V <- readRDS(here::here("outputs/figures_pipeline_z05/data/ref_validation.rds"))
obs <- as.data.table(V$obs_daily)[, .(dTmax_obs = mean(Delta_obs, na.rm = TRUE)), by = id_plot]
lad <- fread(here::here("in_files/lad_z05/Blois_lad_z05_r25.csv"))[, .(id_plot, x, y)]
D <- merge(obs, lad, by = "id_plot")

# extraction S2 LAI (été) et LiDAR LAI au point
ras <- c(s2_lai = file.path(OAK, "s2lai_summer_res_10_m.tif"),
         lidar_lai = file.path(OAK, "lidarlai_res_10_m.tif"))
pts <- vect(as.data.frame(D[, .(x, y, id_plot)]), geom = c("x","y"), crs = "EPSG:32631")
for (nm in names(ras)) { r <- rast(ras[[nm]]); p <- if (crs(r)!=crs(pts)) project(pts, crs(r)) else pts
  D[[nm]] <- terra::extract(r, p, ID = FALSE)[[1]] }
D <- D[is.finite(s2_lai) & is.finite(lidar_lai) & is.finite(dTmax_obs)]

# stats : saturation (compression + découplage) ; puis pouvoir prédictif LiDAR vs S2
r_sat   <- cor(D$lidar_lai, D$s2_lai)
p_sat   <- cor.test(D$lidar_lai, D$s2_lai)$p.value   # ~0.4 → non significatif
sd_s2   <- sd(D$s2_lai);  sd_li <- sd(D$lidar_lai)
r_li_b  <- cor(D$lidar_lai, D$dTmax_obs)
r_s2_b  <- cor(D$s2_lai,    D$dTmax_obs)
cat(sprintf("\nN=%d | S2~LiDAR LAI r=%.2f | SD S2=%.2f vs LiDAR=%.2f | buffering~LiDAR r=%.2f | buffering~S2 r=%.2f\n",
            nrow(D), r_sat, sd_s2, sd_li, r_li_b, r_s2_b))

# (A) saturation : S2 comprimé/découplé de la structure réelle (LiDAR)
lim <- range(c(D$lidar_lai, D$s2_lai), na.rm = TRUE)
pA <- ggplot(D, aes(lidar_lai, s2_lai)) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", colour = "grey50") +
  geom_point(size = 2.4, colour = "grey25") +
  coord_equal(xlim = lim, ylim = lim) +
  labs(x = "LiDAR LAI (one-sided)", y = "Sentinel-2 LAI (summer)",
       subtitle = sprintf("Optical saturation: S2 compressed (SD %.2f vs LiDAR %.2f), uncorrelated with structure (r=%.2f, ns p=%.2f)", sd_s2, sd_li, r_sat, p_sat)) +
  theme_article()

# (B) pouvoir prédictif du buffering : LiDAR LAI marche, S2 LAI non
B <- rbind(D[, .(lai = lidar_lai, dTmax_obs, source = sprintf("LiDAR LAI  (r=%.2f)", r_li_b))],
           D[, .(lai = s2_lai,    dTmax_obs, source = sprintf("Sentinel-2 LAI  (r=%.2f)", r_s2_b))])
pB <- ggplot(B, aes(lai, dTmax_obs)) +
  geom_hline(yintercept = 0, linetype = "dashed", colour = "grey60") +
  geom_smooth(method = "lm", se = TRUE, colour = "grey30", fill = "grey85", linewidth = 0.7) +
  geom_point(size = 2.2, colour = "#1A9850") +
  facet_wrap(~ source, nrow = 1, scales = "free_x") +
  labs(x = "LAI", y = expression(bar(Delta*T[max])~obs~"HOBO  ("*degree*"C)"),
       subtitle = "Buffering vs LAI: LiDAR carries the structural signal; saturated S2 does not → Chapter 2") +
  theme_article()

p <- patchwork::wrap_plots(pA, pB, ncol = 2, widths = c(1, 1.5))
ggsave_article(file.path(OUT, "fig_annex_s2_optical_illusion"), p, 13.5, 5.2)
cli::cli_alert_success("Saved fig_annex_s2_optical_illusion")
