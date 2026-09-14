# ==============================================================================
# The log-profile wind correction, drawn: U(h+2m)/U(10m) across canopy height, for four
# assumed ERA5 roughness lengths. This is the correction that is in the chapter's baseline.
#
# Reads : nothing on disk. The curves are analytic (d = 0.7h, z0 = 0.1h); the 0.44 m value
#         is the CDS forecast_surface_roughness for the Blois cell, recorded inline below.
# Writes: outputs/figures_chap1/G1_wind_profile_correction.{png,pdf}
#   Rscript c1_wind_profile_correction.R
# ==============================================================================
# Appendix G — theoretical log-wind-profile correction factor U(h+2m)/U(10m) that
# would map ERA5's 10 m open-field wind to the height just above the canopy, versus
# canopy height h. Three assumed ERA5 roughness lengths z0,ERA = 0.01 / 0.1 / 0.40 m;
# 0.40 m is the value ACTUALLY used at Blois (CDS forecast_surface_roughness, summer
# 2021, range 0.30-0.58 m shown as a band). NOT applied in the simulations.
#   → outputs/figures_chap1/G1_wind_profile_correction.{png,pdf}
suppressPackageStartupMessages({library(data.table);library(ggplot2);library(ggtext);source("scripts/_article_style.R")})
h <- seq(1, 40, by=0.25); d <- 0.7*h; z0 <- 0.1*h                # canopy displacement + roughness
#' Log-profile wind ratio U(h+2m)/U(10m) over the canopy-height grid `h`
#' @param z0ERA assumed ERA5 roughness length (m) for the open-field 10 m wind
#' @return numeric vector, one ratio per height in `h`
ratio <- function(z0ERA) (log(2 + h - d) - log(z0)) / (log(10) - log(z0ERA))
LV  <- c("0.01 (grassland)","0.1 (cropland)","0.44 (ERA5, Blois)","1 (closed forest)")
dat <- rbind(data.table(h=h, r=ratio(0.01), z=LV[1]),
             data.table(h=h, r=ratio(0.10), z=LV[2]),
             data.table(h=h, r=ratio(0.44), z=LV[3]),   # ERA5 fsr Blois cell, Jun-Sep 2021 (0.438-0.443, ~constant)
             data.table(h=h, r=ratio(1.00), z=LV[4]))
dat[, z := factor(z, levels=LV)]
PAL <- setNames(c("#c6dbef","#6baed6","#2171b5","#08306b"), LV)
p <- ggplot() +
  geom_hline(yintercept=1, linetype="dotted", colour="grey40") +
  geom_line(data=dat, aes(h, r, colour=z, linewidth=z)) +
  scale_colour_manual(values=PAL, name=expression(z["0,ERA"]~"(m)")) +
  scale_linewidth_manual(values=c(0.7,0.7,1.4,0.7), guide="none") +
  scale_x_log10(breaks=c(1,2,5,10,20,40)) +
  annotation_logticks(sides="b", colour="grey45", size=0.4,
                      short=unit(3,"pt"), mid=unit(5,"pt"), long=unit(7,"pt")) +
  coord_cartesian(ylim=c(0, 2.4)) +
  labs(x="Canopy height h (m)", y=expression(U["h+2m"]~"/"~U["10m"])) +
  theme_article(12) + legend_corner(0.99, 0.99)
ggsave_article("outputs/figures_chap1/G1_wind_profile_correction", p, 6.6, 4.4)
cat("factor over 15-35 m:\n")
for(ze in c(0.01,0.1,0.44,1)) cat(sprintf("  z0,ERA=%.2f : h15=%.2f  h35=%.2f\n", ze,
  ratio(ze)[which.min(abs(h-15))], ratio(ze)[which.min(abs(h-35))]))
cat("DONE\n")
