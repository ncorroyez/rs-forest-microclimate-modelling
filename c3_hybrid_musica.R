# ==============================================================================
# Chapter 3 — LAYERED-FUSION dynamic LAI as a MuSICA scenario, validated across
# seasonal windows (Blois). The hybrid gives each layer its own phenology:
#   top    (S2-visible, d_opt) -> the OBSERVED S2 phenology (frac_top)
#   bottom (S2-invisible)      -> a WIDER understory phenology (leaf-out earlier,
#                                 senescence later) — stretch factor s.
#   hybrid_LAI(pid,doy) = LAI_dopt(pid)*frac_top(pid,doy) + LAI_below(pid)*frac_bot(doy)
# If frac_bot == frac_top it reduces EXACTLY to DYN_S2_ANNUAL (full x S2 shape);
# a wider bottom adds leaf area at the SHOULDERS — the seasonal test.
# Compares DYN_HYBRID (s=1.0 and s=1.2) vs DYN_S2_ANNUAL / CONST_ALS / STATIC_S2_ATBD.
# Run from z_Example root:  Rscript c3_hybrid_musica.R
# ==============================================================================
suppressPackageStartupMessages({
  library(terra); library(sf); library(ncdf4); library(lubridate)
  library(dplyr); library(tidyr); library(stringr); library(purrr)
  library(rmusica); library(musica.tools)
})
src <- list.files("R", pattern="\\.R$", full.names=TRUE); src <- src[!grepl("/(h1_|lovb_)", src)]
invisible(lapply(src, source)); source("Chapter3_config.R")
prep <- load_lai_prep(CFG_C3); df <- prep$df_plots; tb <- prep$ts_by_plot

pid_of <- function(r) sprintf("X%d_Y%d", round(r$x), round(r$y))
key <- df %>% mutate(pid = sprintf("X%d_Y%d", round(x), round(y))) %>%
  transmute(pid, LAI_ALS, dopt = LAI_ALS_DOPT, below = pmax(LAI_ALS - LAI_ALS_DOPT, 0))

# stretch the phenology shape horizontally around DOY 200 (s>1 = wider: earlier
# leaf-out + later senescence). frac in [0,1].
stretch_frac <- function(doy, frac, s) {
  ds <- 200 + (doy - 200) / s
  approx(doy, frac, xout = ds, rule = 2)$y
}
build_hybrid <- function(s_bottom) {
  pids <- intersect(names(tb$annual), key$pid)
  setNames(lapply(pids, function(p) {
    a <- tb$annual[[p]]                       # doy, lai = full_LiDAR * frac_top
    k <- key[key$pid == p, ]
    frac_top <- pmin(pmax(a$lai / k$LAI_ALS, 0), 1)
    frac_bot <- pmin(pmax(stretch_frac(a$doy, frac_top, s_bottom), 0), 1)
    data.frame(doy = a$doy, lai = k$dopt * frac_top + k$below * frac_bot)
  }), pids)
}
ts_hyb_10 <- build_hybrid(1.0)     # control: should reproduce DYN_S2_ANNUAL
ts_hyb_12 <- build_hybrid(1.2)     # wider understory phenology

mk <- function(ts_list, nm) make_phenology_fn_factory(ts_list, CFG_C3$list_year)
all_sc <- make_all_scenarios_c3(tb, CFG_C3$list_year, CFG_C3$d_opt_m, CFG_C3$scenarios_mode)
.fn_col <- function(col) function(pr) as.numeric(pr[[col]])
lad_full <- make_lad_real
sc <- list(
  DYN_S2_ANNUAL  = all_sc[["DYN_S2_ANNUAL"]],
  CONST_ALS      = all_sc[["CONST_ALS"]],
  STATIC_S2_ATBD = all_sc[["STATIC_S2_ATBD"]],
  DYN_HYBRID_s10 = list(name="DYN_HYBRID_s10", lai_fn=.fn_col("LAI_ALS"), hmax_fn=.fn_col("Hmax"),
                        fcover_fn=.fn_col("fCover"), lad_fn=lad_full, phenology_fn=mk(ts_hyb_10)),
  DYN_HYBRID_s12 = list(name="DYN_HYBRID_s12", lai_fn=.fn_col("LAI_ALS"), hmax_fn=.fn_col("Hmax"),
                        fcover_fn=.fn_col("fCover"), lad_fn=lad_full, phenology_fn=mk(ts_hyb_12)))
sc <- Filter(Negate(is.null), sc)
nc_parent <- file.path(CFG_C3$out_dir, "nc")

# sanity: how different is the hybrid from DYN_S2_ANNUAL? (shoulder LAI)
p1 <- names(ts_hyb_12)[1]
cat(sprintf("plot %s: DOY110  annual=%.2f hyb_s10=%.2f hyb_s12=%.2f | DOY300 annual=%.2f hyb_s12=%.2f\n",
  p1, tb$annual[[p1]]$lai[110], ts_hyb_10[[p1]]$lai[110], ts_hyb_12[[p1]]$lai[110],
  tb$annual[[p1]]$lai[300], ts_hyb_12[[p1]]$lai[300]))

mkseq <- function(a,b) seq(as.Date(a), as.Date(b), by="day")
WIN <- list(leafout=mkseq("2021-04-01","2021-05-31"), summer=mkseq("2021-06-01","2021-09-30"),
            autumn=mkseq("2021-10-01","2021-10-31"),
            winter=c(mkseq("2021-01-01","2021-03-31"),mkseq("2021-11-01","2021-12-31")),
            fullyear=mkseq("2021-01-01","2021-12-31"))
res <- list()
for (wn in names(WIN)) { ds <- WIN[[wn]]
  dm <- extract_macro_daily(CFG_C3$forcing_file, ds)
  hd <- read_hobo_daily(CFG_C3$hobo_temp_csv, ds, dm, CFG_C3$ids_to_remove)
  val <- validate_scenarios_at_hobos(df, hd, sc, nc_parent, dm, ds,
                                     CFG_C3$forcing_file, CFG_C3$musica_cmd, FALSE)
  m <- val$metrics; m$window <- wn; res[[wn]] <- m
  cat(sprintf("\n--- %s ---\n", wn)); print(as.data.frame(m[,c("scenario","n","r2","rmse","bias")]), row.names=FALSE, digits=3)
}
out <- bind_rows(res)
write.csv(out, file.path(CFG_C3$out_dir,"tables","c3_hybrid_musica_windows.csv"), row.names=FALSE)
write.csv(out, "/home/corroyez/Documents/NC_Full/manuscripts/ch3/tables/Table4c_hybrid_musica_windows.csv", row.names=FALSE)
cat("\nDONE\n")
