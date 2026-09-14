# Does d_opt explain MuSICA's per-plot residual error (ΔTmax)? Tests d_opt at
# the microclimate level (not via LAI). Best scenario = DYN_RF. No new sim.
suppressPackageStartupMessages({
  library(terra); library(sf); library(ncdf4); library(lubridate)
  library(dplyr); library(tidyr); library(stringr); library(purrr)
  library(rmusica); library(musica.tools)
})
src <- list.files("R", pattern="\\.R$", full.names=TRUE); src <- src[!grepl("/(h1_|lovb_)", src)]
invisible(lapply(src, source)); source("Chapter3_config.R")

prep <- load_lai_prep(CFG_C3)
sc <- make_all_scenarios_c3(prep$ts_by_plot, list_year=CFG_C3$list_year,
                            d_opt=CFG_C3$d_opt_m, mode=CFG_C3$scenarios_mode)["DYN_RF"]
ds <- CFG_C3$date_seq
dm <- extract_macro_daily(CFG_C3$forcing_file, ds)
hd <- read_hobo_daily(CFG_C3$hobo_temp_csv, ds, dm, CFG_C3$ids_to_remove)
val <- validate_scenarios_at_hobos(prep$df_plots, hd, sc, file.path(CFG_C3$out_dir,"nc"),
                                   dm, ds, CFG_C3$forcing_file, CFG_C3$musica_cmd, FALSE)

pp <- val$daily %>% group_by(id_plot) %>%
  summarise(err=mean(Delta_sim-Delta_obs), mae=mean(abs(Delta_sim-Delta_obs)), .groups="drop")
hb <- prep$df_plots
idcol <- intersect(c("id_plot","plot_id"), names(hb))[1]
hb$id_plot <- hb[[idcol]]
m <- merge(pp, hb[,c("id_plot","LAI_ALS","LAI_ALS_DOPT","fCover","Hmax","VCI")], by="id_plot")
m$compact <- m$LAI_ALS_DOPT / m$LAI_ALS

cat(sprintf("plots with residual + d_opt: %d\n", nrow(m)))
cat("\n=== correlation of DYN_RF per-plot residual with structural quantities ===\n")
cat("(compact = LAI_ALS_dopt/LAI_ALS = how top-heavy the canopy is)\n")
for (v in c("LAI_ALS_DOPT","compact","LAI_ALS","fCover","Hmax","VCI")) {
  cat(sprintf("  %-14s  cor(bias)=%+.2f  cor(|err|)=%+.2f\n", v,
      cor(m[[v]], m$err, use="complete.obs"), cor(m[[v]], m$mae, use="complete.obs")))
}
