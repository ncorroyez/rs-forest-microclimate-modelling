# ==============================================================================
# gen_gedi_matrix.R — Chapter 4: the 3 missing sensor-combination scenarios to
# complete the summer matrix {ALS, S2, GEDI} x {static, temporal}:
#
#  DYN_ALS_GEDIONLY : ALS + GEDI, no S2. Per-plot LAI_ALS magnitude x the
#      SITE-LEVEL GEDI seasonal shape (Blois adjusted acquisitions, Table14,
#      weekly-interpolated + whit-smoothed, converted to leaf fraction with the
#      April pre-leaf-out floor). Real LAD, ALS Hmax.
#  DYN_S2GEDI       : S2 + GEDI, ALS-FREE. Per-plot GEDI-RF magnitude (Table7)
#      x the plot's own S2 GEDIMAX-corrected shape (dips filled + GEDI fall),
#      normalized by its plateau. NAIVE recipe (FORMS-H height, uniform LAD).
#  DYN_GEDIONLY     : GEDI only, ALS-FREE. GEDI-RF magnitude x the site GEDI
#      shape. NAIVE recipe.
#
# Everything else in the matrix already exists in nc_genuine53_windcorr
# (CONST/STATIC_ALS, STATIC/DYN_S2_ATBD, NAIVE_S2_FORMSH, STATIC_GEDI_*,
# DYN_S2_ANNUAL_FIX, DYN_ALS_GEDIFALL/GEDIMAX, FUSION_H).
#   DRY=1 Rscript gen_gedi_matrix.R ; Rscript gen_gedi_matrix.R
# ==============================================================================
suppressPackageStartupMessages({
  library(ncdf4); library(lubridate); library(dplyr); library(tidyr)
  library(stringr); library(purrr); library(data.table); library(parallel)
  library(rmusica); library(musica.tools)
})
src <- list.files("R", pattern="\\.R$", full.names=TRUE); src <- src[!grepl("/(h1_|lovb_)", src)]
invisible(lapply(src, source)); source("Chapter3_config.R"); source("R/wind_correction.R")
BIN  <- normalizePath("in_files/model-3.2.3/musica", mustWork = TRUE)
FORC_BASE <- "in_files/musica_in_Blois_pblh.nc"; stopifnot(file.exists(FORC_BASE))
ABL  <- list("abl_flag" = '"iter"'); LY <- 2020:2022; NCORES <- 10L
DRY  <- nzchar(Sys.getenv("DRY"))
NCROOT <- file.path(CFG_C3$out_dir, "nc_genuine53_windcorr")
WFDIR  <- file.path(CFG_C3$out_dir, "windcorr_forc_pblh")
NCFULL <- "/home/corroyez/Documents/NC_Full"

windcorr_pblh <- function(hmax){
  f <- wind_factor(pmax(hmax, 2)); ff <- file.path(WFDIR, sprintf("pblh_f%.2f.nc", round(f, 2)))
  if(!file.exists(ff) || file.size(ff) < 1e5){ file.copy(FORC_BASE, ff, overwrite=TRUE)
    nc <- nc_open(ff, write=TRUE); for(v in c("Wind_E","Wind_N")) ncvar_put(nc, v, ncvar_get(nc, v) * f); nc_close(nc) }
  ff
}
inject_seed <- function(ph){ d<-ph[ph$year==2020 & ph$Julian_day==365,,drop=FALSE]
  if(nrow(d)!=1L) stop("missing 2020/365"); d$Julian_day<-366L; rbind(ph,d) }
wrap_iter <- function(sc){ orig<-sc$phenology_fn; lai_fn<-sc$lai_fn
  sc$phenology_fn <- function(pr){ ph<-if(!is.null(orig))orig(pr) else NULL
    if(is.null(ph)) ph<-calc_phenology(list.year=LY,nleafage=1,budburst_date=115,leaf_age_max_in=0.56,
        relative_age_firstmax=0.10,relative_age_lastmax=0.75,LAI_max_per_cohort=as.numeric(lai_fn(pr)))
    inject_seed(ph) }; sc }

# ---- site-level GEDI leaf-fraction shape (Blois, Table14 adjusted) -----------
t14 <- fread(file.path(NCFULL, "chapter4_GEDI/tables/Table14_c4_blois_dates.csv"))
t14 <- t14[pool == "power"][as.Date(date) != as.Date("2021-05-29")]
t14[, doy := as.integer(format(as.Date(date), "%j"))]
setorder(t14, doy)
wk  <- seq(min(t14$doy), max(t14$doy), by = 7)
sm  <- as.numeric(phenofit::whit2(approx(t14$doy, t14$pai_adj, wk)$y, lambda = 10))
flo <- t14[doy == min(doy), pai_adj]                    # April pre-leaf-out floor
frc <- pmin(pmax((sm - flo) / (max(sm) - flo), 0), 1)
f_gedi_shape <- approxfun(
  c(1, min(wk) - 1, wk, 350, 365),
  c(0.02, 0.02, frc, 0.05, 0.05), rule = 2)
cat("GEDI site shape: floor", round(flo, 2), "| max", round(max(sm), 2),
    "| frac @ doy 200/280/318:", round(f_gedi_shape(c(200, 280, 318)), 2), "\n")

# ---- per-plot data ------------------------------------------------------------
prep <- load_lai_prep(CFG_C3); df <- prep$df_plots; tb <- prep$ts_by_plot
df$pid <- sprintf("X%d_Y%d", round(df$x), round(df$y))
t7 <- fread(file.path(NCFULL, "chapter4_GEDI/tables/Table7_c4_perplot_gedi_lai.csv"))
df <- merge(df, t7[, .(id_plot = as.character(id_plot), LAI_GEDI_RF)],
            by = "id_plot", sort = FALSE)
stopifnot(nrow(df) == 53, !anyNA(df$LAI_GEDI_RF))

# GEDIMAX per-plot S2 shape (dips filled + GEDI fall), normalized by plateau
f_gedi_fall <- approxfun(c(250, 282, 318, 350), c(1.00, 0.951, 0.321, 0.08), rule = 2)
shape_s2max <- setNames(lapply(df$pid, function(p) {
  a <- tb$annual[[p]]
  if (is.null(a) || !nrow(a)) return(NULL)
  mag <- mean(a$lai[a$doy %in% 180:240], na.rm = TRUE)
  lai <- a$lai
  core <- a$doy >= 180 & a$doy < 250
  lai[core] <- pmax(lai[core], mag)
  fall <- a$doy >= 250 & a$doy <= 350
  lai[fall] <- mag * f_gedi_fall(a$doy[fall])
  data.frame(doy = a$doy, frac = lai / mag)
}), df$pid)

mk_ts <- function(mag_col, shape) setNames(lapply(seq_len(nrow(df)), function(i) {
  p <- df$pid[i]; mag <- df[[mag_col]][i]
  if (identical(shape, "gedi_site"))
    return(data.frame(doy = 1:365, lai = mag * f_gedi_shape(1:365)))
  s <- shape_s2max[[p]]; if (is.null(s)) return(NULL)
  data.frame(doy = s$doy, lai = mag * s$frac)
}), df$pid)
drop_null <- function(ts) ts[!vapply(ts, is.null, TRUE)]

.fn <- function(col) function(pr) as.numeric(pr[[col]])
mk  <- function(ts) make_phenology_fn_factory(drop_null(ts), LY)
SCEN <- list(
  DYN_ALS_GEDIONLY = wrap_iter(list(name="DYN_ALS_GEDIONLY",
    lai_fn=.fn("LAI_ALS"), hmax_fn=.fn("Hmax"), fcover_fn=.fn("fCover"),
    lad_fn=make_lad_real, phenology_fn=mk(mk_ts("LAI_ALS", "gedi_site")))),
  DYN_S2GEDI = wrap_iter(list(name="DYN_S2GEDI",
    lai_fn=.fn("LAI_GEDI_RF"), hmax_fn=.fn("FORMS_H"), fcover_fn=.fn("fCover"),
    lad_fn=make_lad_uniform, phenology_fn=mk(mk_ts("LAI_GEDI_RF", "s2max")))),
  DYN_GEDIONLY = wrap_iter(list(name="DYN_GEDIONLY",
    lai_fn=.fn("LAI_GEDI_RF"), hmax_fn=.fn("FORMS_H"), fcover_fn=.fn("fCover"),
    lad_fn=make_lad_uniform, phenology_fn=mk(mk_ts("LAI_GEDI_RF", "gedi_site")))))
for (s in names(SCEN)) dir.create(file.path(NCROOT, s), recursive=TRUE, showWarnings=FALSE)

invisible(lapply(unique(df$Hmax), windcorr_pblh))
rows <- if (DRY) seq_len(4) else seq_len(nrow(df))
jobs <- CJ(i = rows, sc = names(SCEN))
run_one <- function(j){
  i <- jobs$i[j]; scn <- jobs$sc[j]; prow <- as.data.frame(df[i, ])
  nc <- file.path(NCROOT, scn, sprintf("musica_out_HOBO_%s.nc", prow$id_plot))
  if (file.exists(nc) && file.size(nc) > 1000) return(sprintf("skip %s %s", scn, prow$id_plot))
  if (file.exists(nc)) file.remove(nc)
  forc <- windcorr_pblh(prow$Hmax)
  ok <- tryCatch({ run_musica_one(prow, SCEN[[scn]], nc, forc, BIN, extra_setup=ABL)
                   file.exists(nc) && file.size(nc) > 1000 }, error=function(e) FALSE)
  sprintf("%s %s %s", if (ok) "OK" else "FAIL", scn, prow$id_plot)
}
cat(sprintf("%s: %d sims on %d cores\n", if (DRY) "DRY" else "FULL", nrow(jobs), NCORES))
t0 <- Sys.time(); res <- unlist(mclapply(seq_len(nrow(jobs)), run_one,
                                         mc.cores=NCORES, mc.preschedule=FALSE))
dt <- as.numeric(difftime(Sys.time(), t0, units="mins")); tab <- table(sub(" .*","",res))
cat(sprintf("DONE %.1f min | %s\n", dt, paste(names(tab), tab, sep="=", collapse=" ")))
fails <- grep("^FAIL", res, value=TRUE)
if (length(fails)) { cat("FAILS:\n"); cat(fails, sep="\n") }
