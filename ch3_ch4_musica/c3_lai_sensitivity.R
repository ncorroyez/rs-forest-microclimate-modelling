# ==============================================================================
# Chapter 3 — JUSTIFICATION figure: how strongly does sub-canopy microclimate
# respond to LAI? Synthetic canopy (fixed Hmax, fCover, UNIFORM LAD), LAI swept
# 1 → 9 by 0.5 (MuSICA receives 2 × LAI = 2 → 18, the model convention).
# y: ΔTmax (Tmax_micro − Tmax_macro, daily) and the hourly buffering slope
# (Tmicro ~ Tmacro). Shows where buffering responds to LAI and where it plateaus.
# Run from z_Example root:  Rscript c3_lai_sensitivity.R
# ==============================================================================
suppressPackageStartupMessages({
  library(terra); library(sf); library(ncdf4); library(lubridate)
  library(dplyr); library(tidyr); library(ggplot2); library(patchwork)
  library(rmusica); library(musica.tools)
})
src <- list.files("R", pattern="\\.R$", full.names=TRUE); src <- src[!grepl("/(h1_|lovb_)", src)]
invisible(lapply(src, source)); source("Chapter3_config.R")
prep <- load_lai_prep(CFG_C3)

HMAX <- 25; FCOVER <- 0.95
tmpl <- prep$df_plots[1, ]; tmpl$Hmax <- HMAX; tmpl$fCover <- FCOVER   # synthetic canopy
ds  <- seq(as.Date("2021-06-01"), as.Date("2021-09-30"), by="day")
dm  <- extract_macro_daily(CFG_C3$forcing_file, ds)
era5 <- build_era5_hourly(CFG_C3$forcing_file, ds)
outdir <- file.path(CFG_C3$out_dir, "nc_laisens"); dir.create(outdir, showWarnings=FALSE, recursive=TRUE)

lais <- seq(1, 9, 0.5)
res <- lapply(lais, function(L) {
  sc <- list(name=sprintf("LAI_%.1f", L),
             lai_fn=function(pr) L, hmax_fn=function(pr) HMAX,
             fcover_fn=function(pr) FCOVER, lad_fn=make_lad_uniform, phenology_fn=NULL)
  nc <- file.path(outdir, sprintf("musica_LAI_%04.1f.nc", L))
  run_musica_one(tmpl, sc, nc, CFG_C3$forcing_file, CFG_C3$musica_cmd)
  if (!file.exists(nc)) return(NULL)
  dT <- extract_deltatmax_one(nc, dm, ds)            # date, Tmax_micro, Delta_Tmax
  sl <- extract_hourly_slope_one(nc, era5, ds)       # $slope
  data.table::data.table(LAI=L, lai_musica=2*L,
    dTmax = mean(dT$Delta_Tmax, na.rm=TRUE),          # Tmax_micro - Tmax_macro (<0 = buffering)
    Tmax_micro = mean(dT$Tmax_micro, na.rm=TRUE),
    slope = if (!is.null(sl)) sl$slope else NA_real_)
})
R <- data.table::rbindlist(Filter(Negate(is.null), res))
cat("=== LAI sensitivity (synthetic canopy, Hmax=25, fCover=0.95, uniform LAD) ===\n")
print(as.data.frame(R), row.names=FALSE, digits=3)
write.csv(R, "/home/corroyez/Documents/NC_Full/manuscripts/ch3/tables/Table0c_lai_sensitivity.csv", row.names=FALSE)

# ---- figure: ΔTmax and buffering slope vs LAI ----
SAT <- 4   # S2 saturation ceiling
dspan <- round(R$dTmax[R$LAI==SAT] - R$dTmax[R$LAI==max(R$LAI)], 1)
g1 <- ggplot(R, aes(LAI, dTmax)) +
  annotate("rect", xmin=SAT, xmax=9.3, ymin=-Inf, ymax=Inf, fill="#F4C7C3", alpha=.45) +
  annotate("text", x=6.6, y=max(R$dTmax)*0.85,
           label=sprintf("Sentinel-2 saturates (caps ~%d)\nyet ΔTmax still spans ~%.1f°C here\n→ only LiDAR resolves it", SAT, dspan),
           size=2.7, colour="#7B241C") +
  geom_hline(yintercept=0, linetype="dotted", colour="grey60") +
  geom_line(colour="#C0392B", linewidth=1) + geom_point(colour="#C0392B", size=2) +
  scale_x_continuous(breaks=1:9) +
  labs(title="Daytime buffering", x=expression(LAI~(m^2/m^2)),
       y=expression(Delta*T[max]~"="~T[max]^micro-T[max]^macro~(degree*C))) +
  theme_minimal(base_size=11)
g2 <- ggplot(R, aes(LAI, slope)) +
  annotate("rect", xmin=SAT, xmax=9.3, ymin=-Inf, ymax=Inf, fill="#C6DBEF", alpha=.45) +
  annotate("text", x=6.6, y=1.03, label="S2-saturation range", size=2.7, colour="#1B3A5B") +
  geom_hline(yintercept=1, linetype="dotted", colour="grey60") +
  geom_line(colour="#2C7FB8", linewidth=1) + geom_point(colour="#2C7FB8", size=2) +
  scale_x_continuous(breaks=1:9) +
  labs(title="Coupling to macroclimate", x=expression(LAI~(m^2/m^2)),
       y=expression(slope~(T[micro]%~%T[macro])~", <1 = decoupled")) +
  theme_minimal(base_size=11)
g <- (g1 | g2) +
  plot_annotation(title="Sub-canopy microclimate responds to LAI across the whole range — including where Sentinel-2 saturates",
    subtitle="Synthetic canopy: Hmax = 25 m, fCover = 0.95, uniform LAD; MuSICA forced with 2×LAI (model convention).",
    theme=theme(plot.subtitle=element_text(size=8.5)))
ggsave("/home/corroyez/Documents/NC_Full/manuscripts/ch3/figures/Fig0_lai_sensitivity.png", g, width=10, height=4.2, dpi=150)
cat("\nDONE\n")
