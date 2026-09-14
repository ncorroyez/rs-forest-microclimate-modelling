# ==============================================================================
# CONTROL of the scan-angle correction against the lidR-computed profiles.
# For each field plot: lidR-standard LAD profile (k=0.5, vertical beams, l=dz) vs the
# scan-angle-corrected profile, with a PER-LAYER <sec theta> (Eq. 2: l = dz*<sec theta>,
# <sec theta> averaged over the beams reaching each layer, i.e. returns at/below layer top).
# Expected (director's control): a slight LAI decrease along the profiles, or equality
# where all beams are near-vertical. Output: figure + table.
# ==============================================================================
suppressMessages({ library(lidR); library(sf); library(data.table); library(ggplot2); library(patchwork) })
source("scripts/_article_style.R"); options(lidR.progress = FALSE)
CTG_DIR <- "/media/corroyez/MyPassport/01_DATA/Blois/LiDAR/leaf_on/2-las_utm"
R_FOOT <- 25; Z0 <- 0.5; DZ <- 0.5; KEXT <- 0.5; HMAX_CAP <- 40

ctg <- readLAScatalog(CTG_DIR); opt_select(ctg) <- "*"; opt_progress(ctg) <- FALSE
g <- st_read("in_files/data_Blois_utm31n.geojson", quiet = TRUE); g <- st_transform(g, 32631)
xy <- st_coordinates(g); ids <- if ("id_plot" %in% names(g)) g$id_plot else seq_len(nrow(g))

one <- function(x, y) {
  las <- tryCatch(suppressMessages(clip_circle(ctg, x, y, R_FOOT)), error = function(e) NULL)
  if (is.null(las) || is.empty(las)) return(NULL)
  las <- tryCatch(filter_duplicates(las), error = function(e) las)
  las <- tryCatch(suppressMessages(normalize_height(las, knnidw())), error = function(e) NULL)
  if (is.null(las)) return(NULL)
  fin <- is.finite(las@data$Z) & las@data$Z >= 0
  Z <- las@data$Z[fin]; sa <- las@data$ScanAngle[fin]
  if (length(Z) < 50 || is.null(sa)) return(NULL)
  d <- tryCatch(lidR::LAD(Z, dz = DZ, k = KEXT, z0 = Z0), error = function(e) NULL)
  if (is.null(d) || nrow(d) == 0) return(NULL)
  # per-layer <sec theta>: beams reaching layer top (returns with Z <= z + dz/2)
  secL <- vapply(d$z, function(zz) { m <- Z <= (zz + DZ/2); if (sum(m) < 10) NA_real_ else mean(1/cos(sa[m] * pi/180)) }, numeric(1))
  secL[!is.finite(secL)] <- mean(1/cos(sa * pi/180))
  data.table(z = d$z, lad_lidr = d$lad, lad_corr = d$lad / secL, sec = secL)
}

message("Control on ", nrow(xy), " field plots ...")
res <- rbindlist(lapply(seq_len(nrow(xy)), function(i) { r <- one(xy[i,1], xy[i,2]); if (!is.null(r)) r$id <- ids[i]; r }), fill = TRUE)
fwrite(res, "out_files/Chapter1/tables/tab_scanangle_control.csv")
lai <- res[, .(LAI_lidr = sum(lad_lidr, na.rm=T)*DZ, LAI_corr = sum(lad_corr, na.rm=T)*DZ), by = id]
cat(sprintf("plots: %d ; LAI lidR mean %.2f -> corrected %.2f (%+.1f%%)\n",
            uniqueN(res$id), mean(lai$LAI_lidr), mean(lai$LAI_corr), 100*(mean(lai$LAI_corr)/mean(lai$LAI_lidr)-1)))
agg <- res[, .(lidr = mean(lad_lidr, na.rm=T), corr = mean(lad_corr, na.rm=T),
               sec = mean(sec, na.rm=T), sec_lo = quantile(sec,.1,na.rm=T), sec_hi = quantile(sec,.9,na.rm=T)), by = z][order(z)]
cat(sprintf("per-layer <sec theta>: range [%.3f, %.3f] (near-uniform if flat)\n", min(agg$sec,na.rm=T), max(agg$sec,na.rm=T)))

# --- figure ------------------------------------------------------------------
pa <- ggplot(agg, aes(y = z)) +
  geom_path(aes(x = lidr, colour = "lidR (vertical beams, l=dz)"), linewidth = 1) +
  geom_path(aes(x = corr, colour = "scan-angle corrected (l=dz·<sec θ>)"), linewidth = 1, linetype = "22") +
  scale_colour_manual(values = c("lidR (vertical beams, l=dz)" = "grey30", "scan-angle corrected (l=dz·<sec θ>)" = "#D7191C"), name = NULL) +
  labs(x = expression("Mean leaf-area density  (m"^2*" m"^{-3}*")"), y = "Height (m)") +
  theme_article(12) + theme(legend.position = "bottom")
pb <- ggplot(agg, aes(y = z)) +
  geom_ribbon(aes(xmin = sec_lo, xmax = sec_hi), fill = "grey85") +
  geom_path(aes(x = sec), colour = "#D7191C", linewidth = 1) +
  geom_vline(xintercept = 1, linetype = "dashed", colour = "grey50") +
  labs(x = expression("<sec θ>  per layer"), y = NULL) + theme_article(12)
p <- pa + pb + plot_layout(widths = c(2,1)) + plot_annotation(tag_levels = "a")
OUT <- "outputs/figures_pipeline/annex"; dir.create(OUT, recursive=TRUE, showWarnings=FALSE)
# DISTINCT OUTPUT NAME 2026-07-31. This LAS-based control and the CSV-based
# scripts/fig_scanangle_control_csv.R both used to write fig_scanangle_control.png,
# so whichever ran last silently decided which figure the manuscript shipped. The
# CSV script is the one the ledger names and the one Fig. F1 must come from; this
# diagnostic now writes its own file.
ggsave(file.path(OUT, "fig_scanangle_control_LASdiag.png"), p, width = 9, height = 5.5, dpi = 300, bg = "white")
ggsave(file.path(OUT, "fig_scanangle_control_LASdiag.pdf"), p, width = 9, height = 5.5, device = cairo_pdf)
cat("Saved fig_scanangle_control (png+pdf) + tab_scanangle_control.csv\n")
