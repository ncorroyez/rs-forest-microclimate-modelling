# Diagnostic (supervisor request #3, 2026-09-10): go beyond the 1 m air readout.
# Report the daytime JJAS-mean vertical structure per archetype: canopy-top air,
# 1 m air, and soil-surface temperature, real canopy. Foliage (leaf) temperature
# is NOT emitted by the v3.2.0 output (no Tleaf variable), so it is a stated
# limitation, not an analysis. Pure post-processing (CHS41-Rmerge/no-wind).
# In : out_files/Chapter1/nc_archetype_chs41/P{1..4}_real.nc
# Out: out_files/Chapter1/figures/diag_2026-09-10/diag_full_profile_air_soil.png
suppressPackageStartupMessages({library(ncdf4);library(lubridate);library(data.table);library(ggplot2);library(patchwork)})
source("R/cluster_relabel.R")
ds <- seq(as.Date("2021-06-01"), as.Date("2021-09-30"), by = "day")
grab <- function(p){
  f <- sprintf("out_files/Chapter1/nc_archetype_chs41/%s_real.nc", p)
  nc <- tryCatch(nc_open(f), error = function(e) NULL); if (is.null(nc)) return(NULL)
  on.exit(nc_close(nc))
  tu <- ncatt_get(nc,"time","units")$value; t0 <- as.POSIXct(sub(".*since ","",tu), tz="UTC")
  tt <- t0 + dhours(ncvar_get(nc,"time")); keep <- as.Date(tt) %in% ds & hour(tt) >= 10 & hour(tt) <= 16
  Tk <- ncvar_get(nc,"Tair_z"); rh <- ncvar_get(nc,"relative_height")
  vh <- median(ncvar_get(nc,"veget_height_top"), na.rm = TRUE)
  Ts <- ncvar_get(nc,"T_soil")                      # (time, nsoil); surface = layer 1
  air <- data.table(P = p, z_m = rh*vh, Tair = rowMeans(Tk[,keep,drop=FALSE], na.rm=TRUE) - 273.15)
  soil0 <- mean(Ts[1, keep], na.rm = TRUE) - 273.15
  list(air = air, soil = soil0, vh = vh)
}
res <- lapply(paste0("P",1:4), grab)
air <- rbindlist(lapply(res, `[[`, "air")); air[, P := factor(P, levels = paste0("P",1:4))]
tab <- rbindlist(lapply(seq_along(res), function(i){
  a <- res[[i]]$air; vh <- res[[i]]$vh
  data.table(P = paste0("P",i),
             canopy_top = a$Tair[which.min(abs(a$z_m - vh))],   # z ~ Hmax (true canopy top)
             free_air = a$Tair[which.max(a$z_m)],               # top of the air column (aloft)
             air_1m = a$Tair[which.min(abs(a$z_m - 1))],
             soil_surface = res[[i]]$soil)
}))
cat("=== daytime JJAS-mean temperatures per archetype (°C) ===\n"); print(tab)
cat("NOTE: leaf/foliage temperature is not in the v3.2.0 output (no Tleaf variable);\n")
cat("      reporting it would require recompiling MuSICA. Stated as a limitation.\n")
g <- ggplot(air, aes(Tair, z_m, colour = P)) +
  geom_hline(yintercept = 1, linetype = "dotted", colour = "grey55") +
  geom_path(linewidth = .7) +
  geom_point(data = tab, aes(x = soil_surface, y = -0.5, colour = P), shape = 15, size = 2.4, inherit.aes = FALSE) +
  annotate("text", x = Inf, y = -0.5, label = "soil surface", hjust = 1.05, size = 2.8, colour = "grey40") +
  scale_colour_manual(values = PAL_CLUSTER, name = NULL) +
  facet_wrap(~P, scales = "free_x", nrow = 1) +
  labs(x = "Temperature (°C), daytime JJAS mean", y = "Height (m); soil at −0.5",
       subtitle = "Air profile (line) with soil-surface temperature (square). Leaf temperature not emitted by the model.") +
  theme_bw(base_size = 11) + theme(panel.grid.minor = element_blank(), legend.position = "none",
                                   strip.text = element_text(face = "bold"))
ggsave("out_files/Chapter1/figures/diag_2026-09-10/diag_full_profile_air_soil.png",
       g, width = 12, height = 4, dpi = 300, bg = "white")
cat("DONE -> diag_full_profile_air_soil.png\n")
