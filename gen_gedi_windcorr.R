# ==============================================================================
# gen_gedi_windcorr.R — Chapter 4 Part B: force MuSICA with GEDI-anchored LAI
# at the 53 HOBO plots (Blois), genuine-53 machinery (v3.2.3-iter, per-plot wind
# correction, 2xLAI inside run_musica_one).
#
# Two new scenarios, both strictly ALS-free (NAIVE_S2_FORMSH recipe: FORMS-H
# height + UNIFORM LAD + parametric static phenology), only the LAI magnitude
# source changes:
#   STATIC_GEDI_RF    : per-plot RF(GEDI PAI ~ S2 + FORMS-H) trained on Blois
#                       power-beam footprints (NC_Full c4 matchup, zero ALS)
#   STATIC_GEDI_RATIO : LAI_S2_ATBD x [mean(GEDI PAI)/mean(S2)] (Blois anchor)
# FUSION_GEDI (dense<-GEDI, open<-S2) is assembled POST-HOC at extraction time
# from these nc + existing STATIC_S2_ATBD nc (no extra sims), like FigK.
#
# Fair comparators on the leaderboard: NAIVE_S2_FORMSH (same recipe, raw S2
# magnitude), STATIC_S2_ATBD, STATIC_ALS (upper ref).
#   DRY=1 Rscript gen_gedi_windcorr.R   # 6 plots x 2 scenarios
#   Rscript gen_gedi_windcorr.R         # full 53 x 2
# ==============================================================================
suppressPackageStartupMessages({
  library(terra); library(sf); library(ncdf4); library(lubridate); library(dplyr); library(tidyr)
  library(stringr); library(purrr); library(data.table); library(parallel)
  library(randomForest); library(rmusica); library(musica.tools)
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

# ---- GEDI anchor: train on Blois power-beam footprints (zero ALS) -----------
fp <- readRDS(file.path(NCFULL, "output/intermediate/c4/gedi_matchup_3site.rds"))
fp <- fp[site == "Blois" & power == "full" & !is.na(lai_s2)]
formsh <- rast(file.path(NCFULL, "01_DATA/FORMS-H_Height_10m_cm.tif"))
pts <- vect(fp[, .(x, y)], geom=c("x","y"), crs="EPSG:32631")
fp[, FORMS_H := terra::extract(formsh, project(pts, "EPSG:2154"))[, 2] / 100]
fp <- fp[is.finite(FORMS_H)]
set.seed(42)
rf_gedi <- randomForest(x = fp[, .(LAI_S2 = lai_s2, FORMS_H)], y = fp$pai_gedi, ntree = 400)
ratio_gedi <- fp[, mean(pai_gedi) / mean(lai_s2)]
cat(sprintf("GEDI anchor: n=%d footprints, RF OOB R2=%.3f, ratio=%.3f\n",
            nrow(fp), tail(rf_gedi$rsq, 1), ratio_gedi))

prep <- load_lai_prep(CFG_C3); df <- prep$df_plots
df$pid <- sprintf("X%d_Y%d", round(df$x), round(df$y))
df$LAI_GEDI_RF    <- as.numeric(predict(rf_gedi,
                      data.frame(LAI_S2 = df$LAI_S2_ATBD, FORMS_H = df$FORMS_H)))
df$LAI_GEDI_RATIO <- df$LAI_S2_ATBD * ratio_gedi
cat("per-plot LAI fed (mean): ALS", round(mean(df$LAI_ALS),2),
    "| S2", round(mean(df$LAI_S2_ATBD),2),
    "| GEDI_RF", round(mean(df$LAI_GEDI_RF),2),
    "| GEDI_RATIO", round(mean(df$LAI_GEDI_RATIO),2), "\n")
fwrite(as.data.table(df)[, .(id_plot, pid, LAI_ALS, LAI_S2_ATBD, FORMS_H, Hmax,
                             LAI_GEDI_RF, LAI_GEDI_RATIO)],
       file.path(NCFULL, "chapter4_GEDI/tables/Table7_c4_perplot_gedi_lai.csv"))

# ---- scenarios (NAIVE recipe: FORMS-H height + uniform LAD, static pheno) ---
.fn <- function(col) function(pr) as.numeric(pr[[col]])
SCEN <- list(
  STATIC_GEDI_RF    = wrap_iter(list(name="STATIC_GEDI_RF",
    lai_fn=.fn("LAI_GEDI_RF"),    hmax_fn=.fn("FORMS_H"), fcover_fn=.fn("fCover"),
    lad_fn=make_lad_uniform, phenology_fn=NULL)),
  STATIC_GEDI_RATIO = wrap_iter(list(name="STATIC_GEDI_RATIO",
    lai_fn=.fn("LAI_GEDI_RATIO"), hmax_fn=.fn("FORMS_H"), fcover_fn=.fn("fCover"),
    lad_fn=make_lad_uniform, phenology_fn=NULL)))
for (s in names(SCEN)) dir.create(file.path(NCROOT, s), recursive=TRUE, showWarnings=FALSE)

invisible(lapply(unique(df$Hmax), windcorr_pblh))
rows <- if (DRY) seq_len(6) else seq_len(nrow(df))
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
cat(sprintf("%s: %d sims on %d cores\n", if (DRY) "DRY RUN" else "FULL RUN",
            nrow(jobs), NCORES))
t0 <- Sys.time(); res <- unlist(mclapply(seq_len(nrow(jobs)), run_one,
                                         mc.cores=NCORES, mc.preschedule=FALSE))
dt <- as.numeric(difftime(Sys.time(), t0, units="mins")); tab <- table(sub(" .*","",res))
cat(sprintf("DONE %.1f min | %s\n", dt, paste(names(tab), tab, sep="=", collapse=" ")))
fails <- grep("^FAIL", res, value=TRUE)
if (length(fails)) { cat("FAILS:\n"); cat(fails, sep="\n"); cat("\n") }
