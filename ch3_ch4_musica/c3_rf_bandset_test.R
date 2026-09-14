# ==============================================================================
# Band-set control for the raw-bands random forest (Chapter 3, Appendix B).
#
# The chapter compares a random forest trained on the TEN raw Sentinel-2 bands
# against the PROSAIL-retrieved LAI, and reads the gain as "part of the deficit
# arises in the retrieval". The two models do not receive the same input: the
# inversion runs on THREE bands only (B03 green, B04 red, B08 near infrared,
# Main_02_produceLAI.R:26). The gain therefore mixes two causes, the inversion
# itself and its restricted band set.
#
# This script isolates the band set by re-running the SAME recipe on B03/B04/B08
# alone (LOO over 53 plots, ntree = 500, seed 42) and bootstrapping the PAIRED
# difference in dense-stratum R2 on shared plot resamples. Only the 10-band
# minus 3-band contrast isolates the band set: a 3-band forest is still trained
# directly on the observed offset, whereas the retrieval goes bands -> physical
# LAI, so "3 bands vs retrieved LAI" would still conflate the two.
#
# Primary stratum: LAI_ALS >= 3.86 (n = 25), the split the published sentence
# rests on. Secondary: archetype P4 (n = 20), the framing the chapter is moving
# to. mtry is swept (1, 2, 3) because at three features mtry = 2 samples a much
# larger share of the feature space than it does at ten.
#
# Read-only on 03_RESULTS. Out: chapter3_S2_LAI/tables/TableB_rf_bandset.csv
#   Rscript c3_rf_bandset_test.R
# ==============================================================================
suppressPackageStartupMessages({
  library(terra); library(ncdf4); library(lubridate)
  library(data.table); library(randomForest)
})
source("Chapter3_config.R")

TAB   <- "/home/corroyez/Documents/NC_Full/manuscripts/ch3/tables"
REFL  <- paste0("/home/corroyez/Documents/NC_Full/03_RESULTS/Blois/",
                "L2A_T31TCN_A031222_20210614T105443/Reflectance/res_10_m/",
                "L2A_T31TCN_A031222_20210614T105443_Refl")
FORC  <- "in_files/musica_in_Blois_pblh.nc"
CLUST <- "out_files/Chapter1/tables/tab_hobo_perplot_cluster.csv"
CROSS <- 3.86
ds    <- seq(as.Date("2021-06-01"), as.Date("2021-09-30"), by = "day")
BND   <- c("B02","B03","B04","B05","B06","B07","B08","B8A","B11","B12")
BND3  <- c("B03","B04","B08")          # exactly what the inversion receives

# ---- plots + raw-band extraction (identical to c3_rf_rawbands.R) -------------
df  <- as.data.table(readRDS(file.path(CFG_C3$out_dir, "lai_prep",
                                       "df_plots_real53.rds")))
r   <- rast(REFL); names(r) <- BND
pts <- vect(as.data.frame(df[, .(x, y)]), geom = c("x", "y"), crs = "EPSG:32631")
ex  <- as.data.table(terra::extract(r, pts))[, -1]
df  <- cbind(df, ex)
cat(sprintf("extracted %d bands at %d plots; incomplete rows: %d\n",
            ncol(ex), nrow(df), sum(!complete.cases(ex))))

# ---- target `do` = mean(Tmax_obs - Tmax_macro) over JJAS ---------------------
ncf <- nc_open(FORC)
tu  <- ncatt_get(ncf, "time", "units")$value
th  <- ncvar_get(ncf, "time")
t0  <- as.POSIXct(sub(".*since ", "", tu), tz = "UTC")
MA  <- data.table(time = floor_date(t0 + th * 3600, "hour"),
                  Tm = as.numeric(ncvar_get(ncf, "Tair")) - 273.15)
nc_close(ncf)
MA   <- MA[as.Date(time) %in% ds]
mday <- MA[, .(Tmx = max(Tm, na.rm = TRUE)), by = .(date = as.Date(time))]
hb   <- as.data.table(read.csv(CFG_C3$hobo_temp_csv))
hb[, time := as.POSIXct(datetime, format = "%Y-%m-%d %H:%M:%S", tz = "UTC")]
hb <- hb[position_sensor == "a" & !(id_plot %in% CFG_C3$ids_to_remove) &
           as.Date(time) %in% ds,
         .(id_plot = as.character(id_plot), Tobs = t_hobo, date = as.Date(time))]
OBS <- merge(hb[, .(Tmx_o = max(Tobs, na.rm = TRUE)), by = .(id_plot, date)],
             mday, by = "date")[, .(do = mean(Tmx_o - Tmx, na.rm = TRUE)),
                                by = id_plot]

M <- merge(OBS, df, by = "id_plot")
cl <- fread(CLUST)[, .(id_plot = as.character(id_plot), P)]
M  <- merge(M, cl, by = "id_plot")
stopifnot(nrow(M) == 53)
dense <- M$LAI_ALS >= CROSS
p4    <- M$P == "P4"
cat(sprintf("n = %d plots; dense (LAI_ALS >= %.2f) n = %d; P4 n = %d\n",
            nrow(M), CROSS, sum(dense), sum(p4)))

# ---- LOO random forest, verbatim recipe of c3_rf_rawbands.R ------------------
r2 <- function(a, b) {
  ok <- is.finite(a) & is.finite(b)
  if (sum(ok) < 4) return(NA_real_)
  cor(a[ok], b[ok])^2
}
loo_rf <- function(feat, mtry = 2L, seed = 42L) {
  X <- as.data.frame(M[, ..feat]); y <- M$do
  set.seed(seed); p <- rep(NA_real_, nrow(X))
  for (i in seq_len(nrow(X))) {
    rf <- randomForest(x = X[-i, , drop = FALSE], y = y[-i],
                       ntree = 500, mtry = mtry)
    p[i] <- as.numeric(predict(rf, X[i, , drop = FALSE]))
  }
  p
}
# paired bootstrap: the SAME resampled plots score both feature sets, so the
# CI is on the difference and not on two independent quantities.
boot_pair <- function(obs, pa, pb, sel, B = 2000L, seed = 7L) {
  o <- obs[sel]; a <- pa[sel]; b <- pb[sel]
  set.seed(seed)
  d <- replicate(B, {
    i <- sample(length(o), replace = TRUE)
    r2(o[i], a[i]) - r2(o[i], b[i])
  })
  as.numeric(quantile(d, c(.025, .975), na.rm = TRUE))
}
boot_r2 <- function(obs, p, sel, B = 2000L, seed = 7L) {
  o <- obs[sel]; q <- p[sel]; set.seed(seed)
  as.numeric(quantile(replicate(B, { i <- sample(length(o), replace = TRUE)
                                     r2(o[i], q[i]) }), c(.025, .975), na.rm = TRUE))
}

PRED <- list()
rows <- rbindlist(lapply(list(
    list(lab = "bands10", feat = BND,  mtry = 2L),
    list(lab = "bands3",  feat = BND3, mtry = 2L),
    list(lab = "bands3",  feat = BND3, mtry = 1L),
    list(lab = "bands3",  feat = BND3, mtry = 3L),
    list(lab = "bands10", feat = BND,  mtry = 3L)
  ), function(v) {
  key <- paste0(v$lab, "_mtry", v$mtry)
  p   <- loo_rf(v$feat, mtry = v$mtry)
  PRED[[key]] <<- p
  ci_d <- boot_r2(M$do, p, dense); ci_4 <- boot_r2(M$do, p, p4)
  out <- data.table(features = key, n_feat = length(v$feat), mtry = v$mtry,
    pooled_r2 = r2(M$do, p),
    dense_r2 = r2(M$do[dense], p[dense]),
    dense_lo = ci_d[1], dense_hi = ci_d[2],
    p4_r2 = r2(M$do[p4], p[p4]), p4_lo = ci_4[1], p4_hi = ci_4[2])
  cat(sprintf("%-16s pooled=%.3f  dense=%.3f [%.2f,%.2f]  P4=%.3f [%.2f,%.2f]\n",
      key, out$pooled_r2, out$dense_r2, ci_d[1], ci_d[2], out$p4_r2,
      ci_4[1], ci_4[2]))
  out
}))
fwrite(rows, file.path(TAB, "TableB_rf_bandset.csv"))

# ---- the number that matters: paired Delta R2 (10 bands - 3 bands) -----------
cat("\n=== paired Delta R2 (10 bands - 3 bands), same resampled plots ===\n")
cmp <- rbindlist(lapply(list(c("bands10_mtry2", "bands3_mtry2"),
                             c("bands10_mtry2", "bands3_mtry1"),
                             c("bands10_mtry2", "bands3_mtry3"),
                             c("bands10_mtry3", "bands3_mtry3")),
  function(k) {
    dd <- boot_pair(M$do, PRED[[k[1]]], PRED[[k[2]]], dense)
    d4 <- boot_pair(M$do, PRED[[k[1]]], PRED[[k[2]]], p4)
    o <- data.table(a = k[1], b = k[2],
      d_dense = r2(M$do[dense], PRED[[k[1]]][dense]) -
                r2(M$do[dense], PRED[[k[2]]][dense]),
      dense_lo = dd[1], dense_hi = dd[2],
      d_p4 = r2(M$do[p4], PRED[[k[1]]][p4]) - r2(M$do[p4], PRED[[k[2]]][p4]),
      p4_lo = d4[1], p4_hi = d4[2])
    cat(sprintf("%-14s - %-14s  dense dR2=%+.3f [%+.2f,%+.2f]  P4 dR2=%+.3f [%+.2f,%+.2f]\n",
        k[1], k[2], o$d_dense, dd[1], dd[2], o$d_p4, d4[1], d4[2]))
    o
  }))
fwrite(cmp, file.path(TAB, "TableB_rf_bandset_delta.csv"))
fwrite(cbind(data.table(id_plot = M$id_plot, obs = M$do, lai_als = M$LAI_ALS,
                        P = M$P), as.data.table(PRED)),
       file.path(TAB, "TableB_rf_bandset_preds.csv"))
cat("\nreference (chapter): retrieved S2 LAI dense R2 = 0.01, LiDAR = 0.56\n")
cat("wrote TableB_rf_bandset{,_delta,_preds}.csv\n")
