# Diagnostic (supervisor request #1, 2026-09-10): do the simulated vertical
# air-temperature profiles differ by LAI, and do they converge at 1 m (the
# sensor readout height)? Fix the Gaussian peak (mu) and spread (sigma) at the
# design centre, vary one-sided LAI 2..6 (cover 0.87, Hmax 25 m, renormalised).
# Pure post-processing of on-disk H2 gaussian-grid NetCDFs (CHS41-Rmerge/no-wind).
# In : out_files/H2_gaussian_grid_chs41/LAI{2..6}_mu050_s024.nc
# Out: out_files/Chapter1/figures/diag_2026-09-10/diag_vertical_Tprofile_by_LAI.png
suppressPackageStartupMessages({library(ncdf4);library(lubridate);library(data.table);library(ggplot2)})
ds <- seq(as.Date("2021-06-01"), as.Date("2021-09-30"), by = "day")
prof <- function(f){
  nc <- tryCatch(nc_open(f), error = function(e) NULL); if (is.null(nc)) return(NULL)
  on.exit(nc_close(nc))
  tu <- ncatt_get(nc,"time","units")$value; t0 <- as.POSIXct(sub(".*since ","",tu), tz="UTC")
  tt <- t0 + dhours(ncvar_get(nc,"time"))
  Tk <- ncvar_get(nc,"Tair_z"); rh <- ncvar_get(nc,"relative_height")
  vh <- median(ncvar_get(nc,"veget_height_top"), na.rm = TRUE)
  keep <- as.Date(tt) %in% ds & hour(tt) >= 10 & hour(tt) <= 16
  data.table(z_m = rh*vh, Tair = rowMeans(Tk[,keep,drop=FALSE], na.rm = TRUE) - 273.15)
}
MU <- "050"; SG <- "024"   # design-centre peak height and spread
D <- rbindlist(lapply(2:6, function(l){
  f <- sprintf("out_files/H2_gaussian_grid_chs41/LAI%d_mu%s_s%s.nc", l, MU, SG)
  d <- prof(f); if (is.null(d)) return(NULL); d[, LAI := l]
}))
if (!nrow(D)) stop("no H2 grid NetCDFs found for mu", MU, " s", SG)
D[, LAI := factor(LAI)]
# value at the 1 m readout per LAI, and the spread across LAI at 1 m vs at the
# TRUE canopy top (z ~ Hmax = 25 m; note the air column extends to ~43 m of free air).
HTOP <- 25
at1 <- D[, .(T1m = Tair[which.min(abs(z_m - 1))]), by = LAI]
attop<- D[, .(Ttop = Tair[which.min(abs(z_m - HTOP))]), by = LAI]
cat("=== air temperature at the 1 m readout, by LAI (mu=0.5, sigma=0.24) ===\n"); print(at1)
cat(sprintf("spread across LAI AT 1 m            : %.3f °C (max-min)\n", diff(range(at1$T1m))))
cat(sprintf("spread across LAI AT canopy top (%d m): %.3f °C (max-min)\n", HTOP, diff(range(attop$Ttop))))
g <- ggplot(D, aes(Tair, z_m, colour = LAI, group = LAI)) +
  geom_hline(yintercept = 1, linetype = "dotted", colour = "grey55") +
  annotate("text", x = -Inf, y = 1.4, label = "1 m readout", hjust = -0.05, size = 3, colour = "grey40") +
  geom_path(linewidth = .7) +
  scale_colour_viridis_d(option = "C", end = .9, name = "one-sided LAI") +
  labs(x = "Air temperature (°C), daytime JJAS mean", y = "Height (m)",
       subtitle = "Gaussian LAD, peak mu=0.5 Hmax, spread sigma=0.24, cover 0.87, Hmax 25 m") +
  theme_bw(base_size = 12) + theme(panel.grid.minor = element_blank())
ggsave("out_files/Chapter1/figures/diag_2026-09-10/diag_vertical_Tprofile_by_LAI.png",
       g, width = 6.5, height = 5.5, dpi = 300, bg = "white")
cat("DONE -> diag_vertical_Tprofile_by_LAI.png\n")
