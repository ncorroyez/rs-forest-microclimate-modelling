# ==============================================================================
# Where the summer Sentinel-2 LAI decline lives, band by band (Blois 2021).
#
# DIAGNOSTIC_s2_summer_decline.md establishes that the retrieved LAI falls ~39 %
# from mid-June to late August, coherently on 59 of 60 plots, before real
# senescence, and that the ATBD/opt ratio drifts (1.98 -> 2.53). It names the
# test that separates the candidate causes but does not run it: look at the
# reflectance itself, band by band.
#
# The three signatures are distinguishable:
#   - every band moves by a similar relative amount -> a scene-level factor
#     (illumination geometry, atmospheric correction), not the canopy;
#   - visible and red edge move while NIR (B08/B8A) stays flat -> pigmentation
#     (chlorophyll decline, brown pigments), i.e. a retrieval degeneracy;
#   - red rises while NIR falls -> real loss of leaf area.
# SWIR (B11/B12) is included because it responds to canopy water content, which
# turns the "summer 2021 was wet, so drought is implausible" argument into a
# measurement.
#
# Relative change is used throughout because it is invariant to a multiplicative
# scene factor; the fraction of plots agreeing in sign is the robustness check,
# mirroring the 59/60 already established on the LAI product.
#
# Read-only on 03_RESULTS. Out: chapter3_S2_LAI/tables/TableB_s2_decline_bands.csv
#   Rscript c3_s2_decline_bands.R
# ==============================================================================
suppressPackageStartupMessages({ library(terra); library(data.table) })
source("Chapter3_config.R")

TAB  <- "/home/corroyez/Documents/NC_Full/manuscripts/ch3/tables"
ROOT <- "/home/corroyez/Documents/NC_Full/03_RESULTS/Blois"
BND  <- c("B02","B03","B04","B05","B06","B07","B08","B8A","B11","B12")
INV  <- c("B03","B04","B08")           # the bands the inversion receives

SCN <- data.table(
  date  = as.Date(c("2021-06-14","2021-07-19","2021-08-26","2021-09-22")),
  scene = c("L2A_T31TCN_A031222_20210614T105443",
            "L2A_T31TCN_A022814_20210719T104622",
            "L2A_T31TCN_A032266_20210826T110450",
            "L2A_T31TCN_A032652_20210922T105348"))

df  <- as.data.table(readRDS(file.path(CFG_C3$out_dir, "lai_prep",
                                       "df_plots_real53.rds")))
pts <- vect(as.data.frame(df[, .(x, y)]), geom = c("x", "y"), crs = "EPSG:32631")

# ---- extract the ten bands at the 53 plots, on each summer scene -------------
ref <- rbindlist(lapply(seq_len(nrow(SCN)), function(i) {
  p <- file.path(ROOT, SCN$scene[i], "Reflectance", "res_10_m",
                 paste0(SCN$scene[i], "_Refl"))
  stopifnot(file.exists(p))
  r <- rast(p); stopifnot(nlyr(r) == length(BND)); names(r) <- BND
  v <- as.data.table(terra::extract(r, pts))[, -1]
  cat(sprintf("%s: %d plots, incomplete rows %d, median B08 = %.0f\n",
              SCN$date[i], nrow(v), sum(!complete.cases(v)),
              median(v$B08, na.rm = TRUE)))
  cbind(data.table(date = SCN$date[i], id_plot = df$id_plot,
                   lai_als = df$LAI_ALS), v)
}))
stopifnot(nrow(ref) == 53 * nrow(SCN), !anyNA(ref[, ..BND]))

# ---- per-plot relative change, for each interval and each band --------------
L <- melt(ref, id.vars = c("date","id_plot","lai_als"), measure.vars = BND,
          variable.name = "band", value.name = "rho")
W <- dcast(L, id_plot + lai_als + band ~ date, value.var = "rho")
dates <- as.character(SCN$date)

chg <- function(t0, t1, label) {
  d <- W[, .(id_plot, lai_als, band,
             rel = (get(t1) - get(t0)) / get(t0))]
  d[, .(interval = label,
        med_rel_pct = 100 * median(rel),
        q25 = 100 * quantile(rel, .25), q75 = 100 * quantile(rel, .75),
        frac_same_sign = max(mean(rel > 0), mean(rel < 0)),
        n = .N), by = band]
}
res <- rbind(chg(dates[1], dates[3], "Jun 14 -> Aug 26"),
             chg(dates[1], dates[2], "Jun 14 -> Jul 19"),
             chg(dates[2], dates[3], "Jul 19 -> Aug 26"),
             chg(dates[1], dates[4], "Jun 14 -> Sep 22"))
res[, inversion_band := band %in% INV]
setcolorder(res, c("interval","band","inversion_band","med_rel_pct","q25","q75",
                   "frac_same_sign","n"))
fwrite(res, file.path(TAB, "TableB_s2_decline_bands.csv"))

for (iv in unique(res$interval)) {
  cat(sprintf("\n=== %s ===\n", iv))
  cat(sprintf("%-5s %7s %18s %9s\n", "band", "med %", "IQR %", "sign agr"))
  for (k in seq_len(nrow(res[interval == iv]))) {
    x <- res[interval == iv][k]
    cat(sprintf("%-5s %+7.1f  [%+6.1f,%+6.1f]  %6.0f%%%s\n", x$band,
                x$med_rel_pct, x$q25, x$q75, 100 * x$frac_same_sign,
                if (x$inversion_band) "  <- inversion" else ""))
  }
}

# ---- the discriminating contrast: visible/red-edge against NIR --------------
cat("\n=== contrast (Jun 14 -> Aug 26) ===\n")
a <- res[interval == "Jun 14 -> Aug 26"]
vis <- a[band %in% c("B02","B03","B04"), median(med_rel_pct)]
re  <- a[band %in% c("B05","B06","B07"), median(med_rel_pct)]
nir <- a[band %in% c("B08","B8A"),      median(med_rel_pct)]
swir<- a[band %in% c("B11","B12"),      median(med_rel_pct)]
cat(sprintf("visible %+.1f %%   red edge %+.1f %%   NIR %+.1f %%   SWIR %+.1f %%\n",
            vis, re, nir, swir))
spread <- diff(range(a$med_rel_pct))
cat(sprintf("spread across the ten bands: %.1f points\n", spread))
cat("(raw change only; see the invariant-target stage below before reading it)\n")
cat(sprintf("\nwrote %s\n", file.path(TAB, "TableB_s2_decline_bands.csv")))

# ---- second stage: is the change canopy-dependent, or scene-wide? -----------
# A scene-level factor (aerosol residual in the atmospheric correction, sun
# geometry) moves every plot by the same relative amount whatever its canopy.
# A canopy change scales with how much canopy there is. So regress each plot's
# relative change on its LiDAR leaf area, and report the spread across plots.
cat("\n=== canopy dependence of the Jun 14 -> Aug 26 change ===\n")
W2 <- W[, .(id_plot, lai_als, band,
            rel = 100 * (get(dates[3]) - get(dates[1])) / get(dates[1]))]
dep <- W2[, {
  ft <- lm(rel ~ lai_als)
  .(slope_pct_per_lai = unname(coef(ft)[2]),
    p = summary(ft)$coefficients[2, 4],
    r = cor(rel, lai_als),
    cv_abs = sd(rel) / abs(mean(rel)))
}, by = band]
print(dep[, .(band, slope = round(slope_pct_per_lai, 2), r = round(r, 2),
              p = signif(p, 2), cv = round(cv_abs, 2))])

# ---- NDVI, the quantity the retrieval actually leans on ---------------------
cat("\n=== NDVI per date (median over the 53 plots) ===\n")
nd <- ref[, .(NDVI = median((B08 - B04) / (B08 + B04)),
              NDVI_q25 = quantile((B08 - B04) / (B08 + B04), .25),
              NDVI_q75 = quantile((B08 - B04) / (B08 + B04), .75)), by = date]
print(nd[, lapply(.SD, function(x) if (is.numeric(x)) round(x, 3) else x)])
nn <- dcast(ref[, .(date, id_plot, ndvi = (B08 - B04) / (B08 + B04))],
            id_plot ~ date, value.var = "ndvi")
cat(sprintf("NDVI fell on %d of 53 plots between Jun 14 and Aug 26\n",
            sum(nn[[dates[3]]] < nn[[dates[1]]])))
fwrite(dep, file.path(TAB, "TableB_s2_decline_canopy_dependence.csv"))
fwrite(nd,  file.path(TAB, "TableB_s2_decline_ndvi.csv"))
cat("wrote TableB_s2_decline_{canopy_dependence,ndvi}.csv\n")

# ---- third stage: remove the scene factor with an invariant target ----------
# The raw between-date change mixes the canopy with everything that differs
# between two overpasses: sun geometry, view geometry, and the atmospheric
# correction. Surfaces that cannot change over ten weeks measure that common
# factor directly. We take pixels that are non-vegetated on BOTH dates (NDVI
# < 0.30) and bright enough for a ratio to mean something (B04 > 500, which
# drops water), i.e. roofs, roads and bare ground. Crops are excluded by
# construction: they are green in June. Whatever those pixels do between the
# two dates is the scene factor, and dividing it out leaves the canopy.
cat("\n=== scene factor, measured on surfaces that cannot change ===\n")
sc <- function(s) {
  r <- rast(file.path(ROOT, s, "Reflectance", "res_10_m", paste0(s, "_Refl")))
  names(r) <- BND; r
}
ra <- sc(SCN$scene[1]); rb <- sc(SCN$scene[3])
set.seed(1L); ii <- sample(ncell(ra), 120000L)
A <- as.data.table(ra[ii]); B <- as.data.table(rb[ii])
keep <- complete.cases(A) & complete.cases(B); A <- A[keep]; B <- B[keep]
ndvi <- function(D) (D$B08 - D$B04) / (D$B08 + D$B04)
inv  <- ndvi(A) < 0.30 & ndvi(B) < 0.30 & A$B04 > 500 & B$B04 > 500
cat(sprintf("invariant pixels: n = %d of %d sampled\n", sum(inv), nrow(A)))
stopifnot(sum(inv) > 300)

fac <- rbindlist(lapply(BND, function(bd) data.table(band = bd,
  scene_pct = 100 * median((B[[bd]][inv] - A[[bd]][inv]) / A[[bd]][inv]))))
raw <- res[interval == "Jun 14 -> Aug 26", .(band, forest_pct = med_rel_pct)]
cmp <- merge(raw, fac, by = "band", sort = FALSE)
# ratio correction: what the canopy did once the common factor is divided out
cmp[, canopy_pct := 100 * ((1 + forest_pct / 100) / (1 + scene_pct / 100) - 1)]
cmp[, inversion_band := band %in% INV]
cat(sprintf("%-5s %10s %10s %10s\n", "band", "forest %", "scene %", "canopy %"))
for (k in seq_len(nrow(cmp))) with(cmp[k], cat(sprintf(
  "%-5s %+10.1f %+10.1f %+10.1f%s\n", band, forest_pct, scene_pct, canopy_pct,
  if (inversion_band) "  <- inversion" else "")))
fwrite(cmp, file.path(TAB, "TableB_s2_decline_scenecorrected.csv"))

# the scene factor must be near-uniform across bands, otherwise it is not one
sp_scene <- diff(range(cmp$scene_pct))
cat(sprintf("\nscene factor spread across the ten bands: %.1f points (median %+.1f %%)\n",
            sp_scene, median(cmp$scene_pct)))
cv <- cmp[band %in% c("B03","B04"),   median(canopy_pct)]
cn <- cmp[band %in% c("B08","B8A"),   median(canopy_pct)]
cat(sprintf("canopy-only change: visible %+.1f %%   NIR %+.1f %%\n", cv, cn))
cat(if (abs(cn) < 10 && cv > 15)
      "-> visible up, NIR flat once the scene factor is out: pigmentation, not leaf area.\n"
    else if (cv > 0 && cn < -15) "-> red up, NIR down: real leaf-area loss.\n"
    else "-> mixed pattern, read the table.\n")
