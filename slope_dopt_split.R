# Does the d_opt two-layer split of LAI (top vs below d_opt) predict the
# observed buffering slope better than bulk LAI? Physical d_opt use: top canopy
# intercepts more than shaded lower canopy. Compares nested regressions on the
# OBSERVED hourly slope (n HOBO plots). No MuSICA.
suppressPackageStartupMessages({
  library(terra); library(ncdf4); library(lubridate); library(dplyr); library(stringr)
})
src <- list.files("R", pattern="\\.R$", full.names=TRUE); src <- src[!grepl("/(h1_|lovb_)", src)]
invisible(lapply(src, source)); source("Chapter3_config.R")
prep <- load_lai_prep(CFG_C3); ds <- CFG_C3$date_seq
era5 <- build_era5_hourly(CFG_C3$forcing_file, ds)

hobo <- read.csv(CFG_C3$hobo_temp_csv) %>%
  mutate(time=floor_date(as.POSIXct(datetime, format="%Y-%m-%d %H:%M:%S", tz="UTC"),"hour"),
         date=as.Date(time)) %>%
  filter(position_sensor=="a", date %in% ds, !(id_plot %in% CFG_C3$ids_to_remove)) %>%
  group_by(id_plot,time) %>% summarise(t_hobo=mean(t_hobo,na.rm=TRUE),.groups="drop")
hm <- merge(hobo, era5, by="time")
slp <- function(y,x){v<-var(x,na.rm=TRUE); if(is.na(v)||v==0) NA else cov(y,x,use="complete.obs")/v}
obs <- hm %>% group_by(id_plot) %>% summarise(slope=slp(t_hobo,Tair_era5),n=n(),.groups="drop") %>% filter(n>=24)

hb <- prep$df_plots; idcol <- intersect(c("id_plot","plot_id"),names(hb))[1]; hb$id_plot <- hb[[idcol]]
m <- merge(obs, hb[,c("id_plot","LAI_ALS","LAI_ALS_DOPT","VCI","fCover","Hmax")], by="id_plot")
m <- m[complete.cases(m[,c("slope","LAI_ALS","LAI_ALS_DOPT")]),]
m$LAI_top   <- m$LAI_ALS_DOPT
m$LAI_below <- m$LAI_ALS - m$LAI_ALS_DOPT
cat(sprintf("n=%d plots | LAI_top mean=%.2f sd=%.2f | LAI_below mean=%.2f sd=%.2f\n",
            nrow(m), mean(m$LAI_top), sd(m$LAI_top), mean(m$LAI_below), sd(m$LAI_below)))

r2 <- function(f) summary(lm(f, data=m))$r.squared
cat("\n=== nested models for OBSERVED buffering slope ===\n")
cat(sprintf("  M1 slope~LAI_full              R2=%.3f\n", r2(slope~LAI_ALS)))
cat(sprintf("  M2 slope~LAI_top+LAI_below     R2=%.3f   (d_opt split)\n", r2(slope~LAI_top+LAI_below)))
cat(sprintf("  M3 slope~LAI_full+VCI          R2=%.3f\n", r2(slope~LAI_ALS+VCI)))
cat(sprintf("  M4 slope~LAI_top+LAI_below+VCI R2=%.3f   (d_opt split + VCI)\n", r2(slope~LAI_top+LAI_below+VCI)))

cat("\n=== M2 coefficients (do top & below buffer differently?) ===\n")
print(round(summary(lm(slope~LAI_top+LAI_below, data=m))$coefficients,4))
cat("\nF-test: does the d_opt split (M2) improve on bulk LAI (M1)?\n")
print(anova(lm(slope~LAI_ALS,data=m), lm(slope~LAI_top+LAI_below,data=m)))
