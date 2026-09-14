# ==============================================================================
# DECISIVE VCI TEST — is the low VCI (raster mean 0.53) an absolute-elevation bug?
# For each plot footprint (r=25 m), compute lidR::VCI two ways:
#   (1) vci_raw  = VCI on RAW Z (absolute elevation) — the suspected method
#                  (main_Blois_detailed.R:123 runs VCI BEFORE normalize_height)
#   (2) vci_norm = VCI after normalize_height (height above ground) — the FIX
# Hypothesis: vci_raw ≈ 0.53 (matches raster), vci_norm ≈ 0.75–0.8 (Eva Gril).
#   Rscript c1_vci_normalize_check.R
# Out: out_files/Chapter1/tables/tab_vci_normalize_check.csv
# ==============================================================================
suppressPackageStartupMessages({ library(lidR); library(sf); library(data.table) })
CTG <- "/media/corroyez/MyPassport/01_DATA/Blois/LiDAR/leaf_on/2-las_utm"
ctg <- readLAScatalog(CTG); opt_progress(ctg) <- FALSE
g <- st_read("in_files/data_Blois_utm31n.geojson", quiet = TRUE)
co <- data.table(id = g$id_plot, x = g$coord_x_utm31n, y = g$coord_y_utm31n)
co <- co[is.finite(x) & is.finite(y)]
cat(sprintf("plots: %d ; clipping r=25 m ...\n", nrow(co)))

clips <- clip_circle(ctg, co$x, co$y, 25)              # one catalog pass → list of LAS
if (inherits(clips, "LAS")) clips <- list(clips)

res <- rbindlist(lapply(seq_along(clips), function(i) {
  pt <- clips[[i]]
  if (is.null(pt) || npoints(pt) < 200) return(NULL)
  zr <- pt$Z; zr <- zr[is.finite(zr)]
  vci_raw <- tryCatch(VCI(zr, zmax = max(zr)), error = function(e) NA_real_)   # absolute elevation
  ptn <- tryCatch(normalize_height(pt, tin()), error = function(e) NULL)
  vci_norm <- NA_real_
  if (!is.null(ptn)) {
    zn <- ptn$Z; zn <- zn[is.finite(zn) & zn >= 0]
    if (length(zn) > 200) vci_norm <- tryCatch(VCI(zn, zmax = max(zn)), error = function(e) NA_real_)
  }
  data.table(id = co$id[i], n = npoints(pt), elev = round(median(zr), 1),
             hmax = round(max(if (!is.null(ptn)) ptn$Z else NA), 1),
             vci_raw = round(vci_raw, 3), vci_norm = round(vci_norm, 3))
}), fill = TRUE)

fwrite(res, "out_files/Chapter1/tables/tab_vci_normalize_check.csv")
cat("\n=== per-plot (first 10) ===\n"); print(head(res, 10))
cat(sprintf("\n=== MEDIANS over %d plots ===\n", sum(is.finite(res$vci_norm))))
cat(sprintf("  vci_raw  (absolute elevation, suspected bug) : %.3f   [raster mean = 0.53]\n",
            median(res$vci_raw,  na.rm = TRUE)))
cat(sprintf("  vci_norm (height above ground, FIX)          : %.3f   [Eva Gril ~ 0.75-0.8]\n",
            median(res$vci_norm, na.rm = TRUE)))
cat(sprintf("  ratio raw/norm = %.3f   [predicted ln(H)/ln(H+elev) ~ 0.68]\n",
            median(res$vci_raw, na.rm = TRUE) / median(res$vci_norm, na.rm = TRUE)))
