# ==============================================================================
# Chapter 3 — Fig 5 (premise): observed hourly buffering slope ~ total LiDAR LAI.
# OBSERVED only (no MuSICA): slope_obs = lm(t_hobo ~ Tair_era5) per plot.
# English; writes into the NC_Full chapter figures folder.
# Run from z_Example root:  Rscript c3_fig5_premise_slope.R
# ==============================================================================
suppressPackageStartupMessages({ library(ncdf4); library(dplyr); library(lubridate); library(ggplot2); library(rmusica); library(musica.tools) })
src <- list.files("R", pattern="\\.R$", full.names=TRUE); src <- src[!grepl("/(h1_|lovb_)", src)]
invisible(lapply(src, source)); source("Chapter3_config.R")
prep <- load_lai_prep(CFG_C3)
ds   <- seq(as.Date("2021-06-01"), as.Date("2021-09-30"), by = "day")   # summer
era5 <- build_era5_hourly(CFG_C3$forcing_file, ds)

hobo <- read.csv(CFG_C3$hobo_temp_csv) %>%
  mutate(time = floor_date(as.POSIXct(datetime, format="%Y-%m-%d %H:%M:%S", tz="UTC"), "hour"),
         date = as.Date(time)) %>%
  filter(position_sensor=="a", date %in% ds, !(id_plot %in% CFG_C3$ids_to_remove)) %>%
  group_by(id_plot, time) %>% summarise(t_hobo=mean(t_hobo, na.rm=TRUE), .groups="drop")
hm  <- merge(hobo, era5, by="time")
slp <- function(y,x){ v<-var(x,na.rm=TRUE); if(is.na(v)||v==0) NA else cov(y,x,use="complete.obs")/v }
obs <- hm %>% group_by(id_plot) %>%
  summarise(slope_obs=slp(t_hobo, Tair_era5), n=n(), .groups="drop") %>% filter(n>=24)
b <- merge(obs, prep$df_plots[, c("id_plot","LAI_ALS")], by = "id_plot")
b <- b[is.finite(b$slope_obs) & is.finite(b$LAI_ALS), ]
r <- cor(b$slope_obs, b$LAI_ALS)
cat(sprintf("n=%d  r(slope_obs, LAI_ALS)=%.2f\n", nrow(b), r))

g <- ggplot(b, aes(LAI_ALS, slope_obs)) +
  geom_smooth(method = "lm", colour = "black", fill = "grey80") +
  geom_point(colour = "#4C78A8", size = 2.5, alpha = 0.85) +
  geom_hline(yintercept = 1, linetype = "dotted") +
  labs(title = "Premise: buffering is governed by total LiDAR LAI",
       subtitle = sprintf("Observed slope vs LiDAR LAI | r = %.2f (slope < 1 = buffering)", r),
       x = expression(LiDAR~LAI~(m^2/m^2)),
       y = expression(observed~slope~(HOBO %~% ERA5))) +
  theme_minimal(base_size = 12)
for (d in c("/home/corroyez/Documents/NC_Full/manuscripts/ch3/figures",
            "/home/corroyez/Documents/NC_Full/output/figures"))
  ggsave(file.path(d, "Fig5_premise_slope_LAI.png"), g, width = 8, height = 5.5, dpi = 150)
cat("DONE\n")
