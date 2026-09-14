# Buffering SLOPE = hourly Tmicro ~ Tmacro (NOT Tmax). Tests which scenario
# reproduces the observed hourly buffering slope, and whether d_opt explains the
# slope residual. Sim slope via extract_hourly_slope_one; obs slope = HOBO hourly
# t_hobo ~ ERA5 hourly. No new sim.
suppressPackageStartupMessages({
  library(terra); library(sf); library(ncdf4); library(lubridate)
  library(dplyr); library(tidyr); library(stringr); library(purrr)
  library(rmusica); library(musica.tools)
})
src <- list.files("R", pattern="\\.R$", full.names=TRUE); src <- src[!grepl("/(h1_|lovb_)", src)]
invisible(lapply(src, source)); source("Chapter3_config.R")

prep <- load_lai_prep(CFG_C3)
ds   <- CFG_C3$date_seq
era5 <- build_era5_hourly(CFG_C3$forcing_file, ds)

# ---- observed hourly buffering slope per plot: t_hobo ~ Tair_era5 ----
hobo <- read.csv(CFG_C3$hobo_temp_csv) %>%
  mutate(time = floor_date(as.POSIXct(datetime, format="%Y-%m-%d %H:%M:%S", tz="UTC"), "hour"),
         date = as.Date(time)) %>%
  filter(position_sensor=="a", date %in% ds, !(id_plot %in% CFG_C3$ids_to_remove)) %>%
  group_by(id_plot, time) %>% summarise(t_hobo=mean(t_hobo, na.rm=TRUE), .groups="drop")
hm <- merge(hobo, era5, by="time")
slp <- function(y,x){ v<-var(x,na.rm=TRUE); if(is.na(v)||v==0) NA else cov(y,x,use="complete.obs")/v }
obs <- hm %>% group_by(id_plot) %>%
  summarise(slope_obs=slp(t_hobo, Tair_era5), n=n(), .groups="drop") %>% filter(n>=24)
cat(sprintf("obs hourly slope: %d plots, mean=%.2f range=[%.2f,%.2f]  (<1 = buffering)\n\n",
            nrow(obs), mean(obs$slope_obs), min(obs$slope_obs), max(obs$slope_obs)))

# ---- simulated hourly slope per (scenario, plot) ----
want <- c("DYN_RF","STATIC_ALS","STATIC_S2_ATBD","STATIC_ALS_DOPT","DYN_RF_CLHS","DYN_RF_CLHS_DOPT")
ncdir <- file.path(CFG_C3$out_dir, "nc")
sim <- do.call(rbind, lapply(want, function(scn) {
  d <- file.path(ncdir, scn); fs <- list.files(d, pattern="\\.nc$", full.names=TRUE)
  do.call(rbind, lapply(fs, function(f){
    id <- str_extract(basename(f), "(?<=HOBO_).*(?=\\.nc)")
    r <- extract_hourly_slope_one(f, era5, ds); if (is.null(r)) return(NULL)
    data.frame(scenario=scn, id_plot=id, slope_sim=r$slope)
  }))
}))
sl <- merge(sim, obs[,c("id_plot","slope_obs")], by="id_plot")

cat("=== per-scenario reproduction of observed hourly buffering slope ===\n")
res <- sl %>% group_by(scenario) %>%
  summarise(n=sum(!is.na(slope_sim)),
            slope_sim_mean=mean(slope_sim,na.rm=TRUE),
            rmse=sqrt(mean((slope_sim-slope_obs)^2,na.rm=TRUE)),
            bias=mean(slope_sim-slope_obs,na.rm=TRUE),
            r=cor(slope_sim,slope_obs,use="complete.obs"), .groups="drop") %>%
  arrange(rmse)
print(as.data.frame(res), row.names=FALSE, digits=3)

best <- res$scenario[1]
hb <- prep$df_plots; idcol <- intersect(c("id_plot","plot_id"),names(hb))[1]; hb$id_plot <- hb[[idcol]]
b <- merge(sl[sl$scenario==best,], hb[,c("id_plot","LAI_ALS","LAI_ALS_DOPT","fCover","Hmax","VCI")], by="id_plot")
b$resid <- b$slope_sim - b$slope_obs; b$compact <- b$LAI_ALS_DOPT/b$LAI_ALS
cat(sprintf("\n=== slope residual vs structure (best: %s) ===\n", best))
for(v in c("LAI_ALS_DOPT","compact","LAI_ALS","fCover","Hmax","VCI"))
  cat(sprintf("  %-14s cor(resid)=%+.2f\n", v, cor(b[[v]], b$resid, use="complete.obs")))
cat(sprintf("\n=== obs slope vs structure (does buffering track LAI/d_opt at all?) ===\n"))
ob <- merge(obs, hb[,c("id_plot","LAI_ALS","LAI_ALS_DOPT","fCover","Hmax","VCI")], by="id_plot")
ob$compact <- ob$LAI_ALS_DOPT/ob$LAI_ALS
for(v in c("LAI_ALS_DOPT","compact","LAI_ALS","fCover","Hmax","VCI"))
  cat(sprintf("  %-14s cor(slope_obs)=%+.2f\n", v, cor(ob[[v]], ob$slope_obs, use="complete.obs")))
