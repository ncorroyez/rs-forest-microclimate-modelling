# Hottest-days regime: does d_opt matter on heat extremes (where the upper
# canopy intercepts most radiation)? Judged by match-to-HOBO (RMSE), not "true"
# buffering. Re-analyses existing daily ΔTmax; no new sim.
suppressPackageStartupMessages({
  library(terra); library(sf); library(ncdf4); library(lubridate)
  library(dplyr); library(tidyr); library(stringr); library(purrr)
  library(rmusica); library(musica.tools)
})
src <- list.files("R", pattern="\\.R$", full.names=TRUE); src <- src[!grepl("/(h1_|lovb_)", src)]
invisible(lapply(src, source)); source("Chapter3_config.R")

prep <- load_lai_prep(CFG_C3); ds <- CFG_C3$date_seq
want <- c("DYN_RF","STATIC_ALS","STATIC_S2_ATBD","STATIC_ALS_DOPT",
          "DYN_S2_RESCALED_DOPT","DYN_RF_CLHS","DYN_RF_CLHS_DOPT")
sc_all <- make_all_scenarios_c3(prep$ts_by_plot, list_year=CFG_C3$list_year,
                                d_opt=CFG_C3$d_opt_m, mode=CFG_C3$scenarios_mode)
sc <- sc_all[intersect(want, names(sc_all))]
dm <- extract_macro_daily(CFG_C3$forcing_file, ds)
hd <- read_hobo_daily(CFG_C3$hobo_temp_csv, ds, dm, CFG_C3$ids_to_remove)
val <- validate_scenarios_at_hobos(prep$df_plots, hd, sc, file.path(CFG_C3$out_dir,"nc"),
                                   dm, ds, CFG_C3$forcing_file, CFG_C3$musica_cmd, FALSE)
d <- as.data.frame(val$daily)
dm <- as.data.frame(dm)
dm$date <- as.Date(dm$date); d$date <- as.Date(d$date)

# day subsets by macro Tmax (defined from dm, applied via date membership)
q <- quantile(dm$Tmax_macro, c(.75,.90), na.rm=TRUE)
date_all <- dm$date
date25   <- dm$date[dm$Tmax_macro >= q[1]]
date10   <- dm$date[dm$Tmax_macro >= q[2]]
subsets  <- list(all=date_all, hot_top25=date25, hot_top10=date10)
cat(sprintf("hot-day thresholds: top25%% Tmax_macro>=%.1f (%d d)  top10%%>=%.1f (%d d)\n",
            q[1], length(date25), q[2], length(date10)))

cat("\n=== (a) match-to-HOBO RMSE per scenario, by day regime ===\n")
out <- do.call(rbind, lapply(names(subsets), function(sn){
  dd <- d[d$date %in% subsets[[sn]],]
  dd %>% group_by(scenario) %>%
    summarise(rmse=sqrt(mean((Delta_sim-Delta_obs)^2,na.rm=TRUE)), .groups="drop") %>%
    mutate(regime=sn)
}))
tab <- tidyr::pivot_wider(out, names_from=regime, values_from=rmse)
tab <- tab[order(tab$all),]
print(as.data.frame(tab), row.names=FALSE, digits=3)

# (b) does d_opt predict observed buffering better on hot days?
hb <- prep$df_plots; idcol <- intersect(c("id_plot","plot_id"),names(hb))[1]; hb$id_plot <- hb[[idcol]]
cat("\n=== (b) predictor of OBSERVED ΔTmax buffering, by regime (per-plot) ===\n")
for(sn in names(subsets)){
  dd <- d[d$date %in% subsets[[sn]] & d$scenario==want[1],]   # Delta_obs identical across scenarios
  pp <- dd %>% group_by(id_plot) %>% summarise(dtmax_obs=mean(Delta_obs), .groups="drop")
  m <- merge(pp, hb[,c("id_plot","LAI_ALS","LAI_ALS_DOPT")], by="id_plot")
  cat(sprintf("  %-10s cor(ΔTmax_obs, LAI_full)=%+.2f  cor(.., LAI_dopt)=%+.2f\n",
              sn, cor(m$dtmax_obs,m$LAI_ALS), cor(m$dtmax_obs,m$LAI_ALS_DOPT)))
}
