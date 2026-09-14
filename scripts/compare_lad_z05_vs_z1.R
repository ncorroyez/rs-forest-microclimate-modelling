# ==============================================================================
# Compare new z0.5 LAD/LAI (raw 2-las_utm + on-the-fly normalise, captures
# understory) vs the legacy z1 product (pre-normalised cloud, 10m raster) at the
# Blois field plots. Effect on LAI, Hmax, understory, and cluster assignment.
# Output: outputs/figures_pipeline/annex/fig_lad_z05_vs_z1_blois.{png,pdf}
#         tables/tab_lad_z05_vs_z1_blois.csv
# ==============================================================================

suppressMessages({ library(data.table); library(here); library(sf); library(terra)
  library(tidyverse); library(patchwork) })
source(here::here("R/config.R")); source(here::here("R/io.R")); source(here::here("R/validation.R"))
OUT <- here::here("outputs/figures_pipeline/annex"); dir.create(OUT, recursive=TRUE, showWarnings=FALSE)
TAB <- here::here("outputs/figures_pipeline/tables")

# OLD (z1) per-plot from current pipeline inputs
rasters <- load_lidar_rasters(CFG$in_dir, CFG$agg_factor)
old <- as.data.table(build_hobo_inputs(CFG$hobo_geojson, rasters$stack, CFG$ids_to_remove))
hp <- sf::st_read(CFG$hobo_geojson, quiet=TRUE) %>% filter(!id_plot %in% CFG$ids_to_remove)
old$id_plot <- hp$id_plot
old <- old[, .(id_plot, LAI_z1 = LAI, Hmax_z1 = Hmax)]

# NEW (z0.5)
new <- fread(here::here("in_files/lad_z05/Blois_lad_z05_r25.csv"))
ladc <- grep("LAD_Layer_", names(new), value=TRUE); zc <- as.numeric(sub("LAD_Layer_","",ladc))
new[, LAI_under2m := rowSums(.SD[, ..ladc][, .SD, .SDcols = which(zc<2)]) * 0.5]
new2 <- new[, .(id_plot, LAI_z05 = LAI, Hmax_z05 = Hmax, LAI_under2m)]

D <- merge(old, new2, by="id_plot")
D[, dLAI := LAI_z05 - LAI_z1][, dHmax := Hmax_z05 - Hmax_z1]
fwrite(D, file.path(TAB, "tab_lad_z05_vs_z1_blois.csv"))

cli::cli_h2("z0.5 vs z1 (Blois, n={nrow(D)})")
cat(sprintf("LAI : r=%.3f | mean Δ=%+.2f | median z1=%.2f z05=%.2f\n",
            cor(D$LAI_z1,D$LAI_z05), mean(D$dLAI), median(D$LAI_z1), median(D$LAI_z05)))
cat(sprintf("Hmax: r=%.3f | mean Δ=%+.2f m\n", cor(D$Hmax_z1,D$Hmax_z05), mean(D$dHmax)))
cat(sprintf("Understory LAI <2m : mean=%.2f (%.0f%% of z05 LAI) — absent in z1\n",
            mean(D$LAI_under2m), 100*mean(D$LAI_under2m/D$LAI_z05)))

lim <- range(c(D$LAI_z1,D$LAI_z05))
pL <- ggplot(D, aes(LAI_z1, LAI_z05)) + geom_abline(linetype="dashed",colour="grey55") +
  geom_point(size=2.6, colour="#1A9850", alpha=0.8) + coord_equal(xlim=lim,ylim=lim) +
  labs(x="LAI (legacy z1)", y="LAI (new z0.5)",
       subtitle=sprintf("r=%.2f, mean Δ=%+.2f", cor(D$LAI_z1,D$LAI_z05), mean(D$dLAI))) + theme_bw(base_size=14)
limh <- range(c(D$Hmax_z1,D$Hmax_z05))
pH <- ggplot(D, aes(Hmax_z1, Hmax_z05)) + geom_abline(linetype="dashed",colour="grey55") +
  geom_point(size=2.6, colour="#2166AC", alpha=0.8) + coord_equal(xlim=limh,ylim=limh) +
  labs(x="Hmax (legacy z1)", y="Hmax (new z0.5)",
       subtitle=sprintf("r=%.2f, mean Δ=%+.2f m", cor(D$Hmax_z1,D$Hmax_z05), mean(D$dHmax))) + theme_bw(base_size=14)
pU <- ggplot(D, aes(LAI_z05, LAI_under2m)) + geom_point(size=2.6, colour="#D7191C", alpha=0.8) +
  labs(x="LAI (z0.5)", y="Understory LAI (<2 m)", subtitle="newly captured by z0.5") + theme_bw(base_size=14)

p <- (pL | pH | pU) + plot_annotation(title="Blois: new z0.5 (understory-aware) vs legacy z1 LAD/LAI",
       caption="z0.5 from raw 2-las_utm + on-the-fly normalise (25 m footprint) vs legacy z1 (pre-normalised 10 m raster).")
ggsave(file.path(OUT,"fig_lad_z05_vs_z1_blois.png"), p, width=15, height=5, dpi=300, bg="white")
ggsave(file.path(OUT,"fig_lad_z05_vs_z1_blois.pdf"), p, width=15, height=5, device=cairo_pdf)
cli::cli_alert_success("Saved comparison figure + table")
