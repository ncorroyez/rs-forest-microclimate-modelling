# ==============================================================================
# Chapter 3 — RF option B (2B): predict the S2 LAI (target) from LiDAR STRUCTURE
# features, trained OUT-OF-SAMPLE on 5000 pixels with a UNIFORM LAI distribution
# in [2, p98] (Chapter-2 sampling). Two targets: LAI_S2_ATBD and LAI_S2_opt.
# The RF gives a per-plot magnitude; the scenario is DYNAMIC = (RF magnitude) x
# (S2 ATBD temporal shape). Validated on the 53 HOBO (summer).
# Run from z_Example root:  Rscript c3_summer_rf_b.R
# ==============================================================================
suppressPackageStartupMessages({
  library(terra); library(data.table); library(lubridate); library(dplyr); library(ncdf4)
  library(stringr); library(purrr); library(randomForest); library(rmusica); library(musica.tools)
})
src <- list.files("R", pattern="\\.R$", full.names=TRUE); src <- src[!grepl("/(h1_|lovb_)", src)]
invisible(lapply(src, source)); source("Chapter3_config.R")
prep <- load_lai_prep(CFG_C3); df <- as.data.table(prep$df_plots); tb <- prep$ts_by_plot
df[, pid := sprintf("X%d_Y%d", round(x), round(y))]
df[is.na(LAI_ALS_DOPT), LAI_ALS_DOPT := LAI_ALS]
ropt0 <- mean(df$LAI_S2_DOPT,na.rm=TRUE)/mean(df$LAI_S2_ATBD,na.rm=TRUE); df[is.na(LAI_S2_DOPT), LAI_S2_DOPT := LAI_S2_ATBD*ropt0]
feats <- c("LAI_ALS","Hmax","fCover","VCI","LCV"); lai_min <- 2

# ---- pixel feature/target stack (Blois) ----
ind <- CFG_C3$in_dir
R <- list(LAI_ALS=rast(file.path(ind,"lai_z1_res_10_m.tif")),
          Hmax   =rast(file.path(ind,"max_res_10_m.tif")),
          fCover =rast(file.path(ind,"fCover_res_10_m.tif")),
          VCI    =rast(file.path(ind,"vci_res_10_m.tif")),
          LCV    =rast(file.path(ind,"lcv_res_10_m.tif")),
          LAI_S2_ATBD=rast(file.path(ind,"s2lai_summer_atbd_res_10_m.tif")),
          LAI_S2_DOPT=rast(make_lai_paths(CFG_C3)$lai_s2_dopt))
base <- R$LAI_ALS
for (nm in names(R)) if (!compareGeom(base, R[[nm]], stopOnError=FALSE)) R[[nm]] <- resample(R[[nm]], base)
px <- as.data.table(setNames(lapply(R, function(r) as.numeric(values(r))), names(R)))
px <- px[complete.cases(px)]
cat(sprintf("valid pixels: %d\n", nrow(px)))

# ---- Chapter-2 uniform-LAI sampling: LAI in [2, p98], bins of 1, equal per bin ----
samp_uniform <- function(d, n) {
  hi <- as.numeric(quantile(d$LAI_ALS, 0.98, names=FALSE))
  br <- seq(lai_min, hi, by=1.0); if (tail(br,1) < hi) br <- c(br, hi)
  d <- d[LAI_ALS >= lai_min & LAI_ALS <= hi]; d[, bin := cut(LAI_ALS, br, include.lowest=TRUE)]
  per <- ceiling(n/(length(br)-1))
  set.seed(42); d[, .SD[sample(.N, min(.N, per))], by=bin][, bin:=NULL][]
}
train <- samp_uniform(px, 5000)
cat(sprintf("training pixels: %d | LAI_ALS range [%.1f, %.1f]\n", nrow(train), min(train$LAI_ALS), max(train$LAI_ALS)))

# ---- 2 RF (target = S2), trained on the uniform pixel sample ----
fit_predict <- function(target){
  set.seed(42)
  rf <- randomForest(x=as.data.frame(train[, ..feats]), y=train[[target]], ntree=500, mtry=2)
  oob <- 1 - rf$mse[rf$ntree]/var(train[[target]])
  pred <- as.numeric(predict(rf, as.data.frame(df[, ..feats])))
  cat(sprintf("  RF(target=%s) on 5000 uniform px | OOB R²=%.3f | pred 53 HOBO mean=%.2f\n", target, oob, mean(pred)))
  pred
}
df[, RF_S2_ATBD := fit_predict("LAI_S2_ATBD")]
df[, RF_S2_OPT  := fit_predict("LAI_S2_DOPT")]

# ---- DYNAMIC series: RF magnitude x S2 ATBD temporal shape ----
get_atbd <- function(p) tb$atbd[[p]]
shape_ts <- function(p, mag){
  a <- get_atbd(p)
  if (is.null(a) || !all(c("doy","lai") %in% names(a)) || nrow(a)==0)
    return(data.frame(doy=1:365, lai=rep(mag,365)))
  mu <- mean(a$lai, na.rm=TRUE); if (!is.finite(mu) || mu<=0) mu <- 1
  data.frame(doy=a$doy, lai=a$lai/mu*mag)            # ATBD shape, mean rescaled to RF magnitude
}
ts_rf_atbd <- setNames(lapply(seq_len(nrow(df)), function(i) shape_ts(df$pid[i], df$RF_S2_ATBD[i])), df$pid)
ts_rf_opt  <- setNames(lapply(seq_len(nrow(df)), function(i) shape_ts(df$pid[i], df$RF_S2_OPT[i])),  df$pid)
.fn <- function(col) function(pr) as.numeric(pr[[col]]); mk <- function(ts) make_phenology_fn_factory(ts, CFG_C3$list_year)
sc <- list(
  RF_S2_ATBD = list(name="RF_S2_ATBD", lai_fn=.fn("RF_S2_ATBD"), hmax_fn=.fn("Hmax"), fcover_fn=.fn("fCover"), lad_fn=make_lad_real, phenology_fn=mk(ts_rf_atbd)),
  RF_S2_OPT  = list(name="RF_S2_OPT",  lai_fn=.fn("RF_S2_OPT"),  hmax_fn=.fn("Hmax"), fcover_fn=.fn("fCover"), lad_fn=make_lad_real, phenology_fn=mk(ts_rf_opt)))
ds <- seq(as.Date("2021-06-01"), as.Date("2021-09-30"), by="day")
dm <- extract_macro_daily(CFG_C3$forcing_file, ds); hd <- read_hobo_daily(CFG_C3$hobo_temp_csv, ds, dm, CFG_C3$ids_to_remove)
validate_scenarios_at_hobos(as.data.frame(df), hd, sc, file.path(CFG_C3$out_dir,"nc"), dm, ds, CFG_C3$forcing_file, CFG_C3$musica_cmd, FALSE)
fwrite(df[, .(id_plot, RF_S2_ATBD, RF_S2_OPT, LAI_S2_ATBD, LAI_S2_DOPT, LAI_ALS)],
       "/home/corroyez/Documents/NC_Full/manuscripts/ch3/tables/Table4h_rfB_predictions.csv")
cat("\nDONE RF-B MuSICA runs\n")
