# ==============================================================================
# Chapter 3 — SUMMER MuSICA comparison of LAI-magnitude products / corrections,
# validated vs HOBO ΔTmax (2×LAI convention). Static canopy (parametric summer
# phenology), real LAD shape rescaled to each scenario's magnitude.
# Scenarios: S2_ATBD, S2_opt(=LAI_S2_DOPT), ALS, ALS_dopt, S2_rescaled (prorata),
#   ratioA (S2×global ratio), hybrid (h<d_opt→S2 else LiDAR-bottom+calibrated-S2-top), RF.
# Run from z_Example root:  Rscript c3_summer_lai_scenarios.R
# ==============================================================================
suppressPackageStartupMessages({
  library(terra); library(sf); library(ncdf4); library(lubridate)
  library(dplyr); library(tidyr); library(purrr); library(stringr); library(data.table); library(ggplot2)
  library(rmusica); library(musica.tools)
})
src <- list.files("R", pattern="\\.R$", full.names=TRUE); src <- src[!grepl("/(h1_|lovb_)", src)]
invisible(lapply(src, source)); source("Chapter3_config.R")
prep <- load_lai_prep(CFG_C3); df <- as.data.frame(prep$df_plots)
dopt <- CFG_C3$d_opt_m   # 7

# Fill the 6 open/sparse HOBO plots (near-zero LiDAR LAI, canopy 3-5 m) where the
# d_opt products are NA, so ALL 53 plots run in every scenario:
#  - LAI_ALS_DOPT (top-d_opt LiDAR) = full LAI_ALS  (a short/sparse canopy fits within d_opt)
#  - LAI_S2_DOPT  (opt LUT) imputed from ATBD via the global opt/ATBD ratio
na_dopt <- is.na(df$LAI_ALS_DOPT)
df$LAI_ALS_DOPT[na_dopt] <- df$LAI_ALS[na_dopt]
.ropt <- mean(df$LAI_S2_DOPT, na.rm=TRUE) / mean(df$LAI_S2_ATBD, na.rm=TRUE)
na_s2d <- is.na(df$LAI_S2_DOPT)
df$LAI_S2_DOPT[na_s2d] <- df$LAI_S2_ATBD[na_s2d] * .ropt
cat(sprintf("filled %d plots LAI_ALS_DOPT (=full) and %d plots LAI_S2_DOPT (=ATBD*%.3f)\n",
            sum(na_dopt), sum(na_s2d), .ropt))

# precomputed magnitude columns.
# Prorata correction = rescale each raw S2 LUT to the LiDAR full-canopy magnitude by
# its OWN global ratio (mean LAI_ALS / mean LAI_S2_lut). Applied to BOTH LUTs so the
# comparison never feeds raw (magnitude-biased) S2 into MuSICA.
rA   <- mean(df$LAI_ALS, na.rm=TRUE) / mean(df$LAI_S2_ATBD, na.rm=TRUE)   # ATBD -> LiDAR magnitude
rOpt <- mean(df$LAI_ALS, na.rm=TRUE) / mean(df$LAI_S2_DOPT, na.rm=TRUE)   # opt  -> LiDAR magnitude
bT   <- mean(df$LAI_ALS_DOPT, na.rm=TRUE) / mean(df$LAI_S2_ATBD, na.rm=TRUE)  # S2 -> top calib
# SIMPLE summer scenario set (53 HOBO, ΔTmax). Two corrections of S2 ATBD by prorata,
# and one height-switch fusion (S2 opt where the canopy is short enough that S2 sees it
# all, full LiDAR where it is tall and S2 saturates).
df$LAI_ATBD_rfull <- df$LAI_S2_ATBD * rA      # S2 ATBD rescaled to FULL-LiDAR magnitude (ratio-A)
df$LAI_ATBD_rdopt <- df$LAI_S2_ATBD * bT      # S2 ATBD rescaled to TOP-d_opt magnitude (ratio-B)
df$LAI_fusion_h   <- ifelse(df$Hmax < dopt, df$LAI_S2_DOPT, df$LAI_ALS)  # short->S2 opt, tall->full LiDAR
cat(sprintf("ratio->full=%.3f  ratio->dopt=%.3f  | plots Hmax<d_opt (use S2 opt): %d/%d\n",
            rA, bT, sum(df$Hmax<dopt), nrow(df)))

.fn <- function(col) function(pr) as.numeric(pr[[col]])
mk <- function(nm, col) list(name=nm, lai_fn=.fn(col), hmax_fn=.fn("Hmax"),
                             fcover_fn=.fn("fCover"), lad_fn=make_lad_real, phenology_fn=NULL)
sc <- list(
  SUM_S2_ATBD    = mk("SUM_S2_ATBD",   "LAI_S2_ATBD"),       # 1. raw S2 ATBD
  SUM_S2_opt     = mk("SUM_S2_opt",    "LAI_S2_DOPT"),       # 2. raw S2 opt
  SUM_ALS        = mk("SUM_ALS",       "LAI_ALS"),           # 3. raw LiDAR full
  SUM_ALS_dopt   = mk("SUM_ALS_dopt",  "LAI_ALS_DOPT"),      # 4. raw LiDAR d_opt
  SUM_ATBD_rfull = mk("SUM_ATBD_rfull","LAI_ATBD_rfull"),    # 5. S2 ATBD x ratio -> full
  SUM_ATBD_rdopt = mk("SUM_ATBD_rdopt","LAI_ATBD_rdopt"),    # 6. S2 ATBD x ratio -> d_opt
  SUM_fusion_h   = mk("SUM_fusion_h",  "LAI_fusion_h"))      # 7. height-switch fusion

ds <- seq(as.Date("2021-06-01"), as.Date("2021-09-30"), by="day")
dm <- extract_macro_daily(CFG_C3$forcing_file, ds)
hd <- read_hobo_daily(CFG_C3$hobo_temp_csv, ds, dm, CFG_C3$ids_to_remove)
ncp <- file.path(CFG_C3$out_dir, "nc")
val <- validate_scenarios_at_hobos(df, hd, sc, ncp, dm, ds, CFG_C3$forcing_file, CFG_C3$musica_cmd, FALSE)
m <- as.data.table(val$metrics)
# mean fed LAI per scenario (for context)
laimean <- sapply(sc, function(s) mean(sapply(seq_len(nrow(df)), function(i) s$lai_fn(df[i,])), na.rm=TRUE))
m[, mean_LAI := round(laimean[scenario], 2)]
m <- m[order(rmse)]
cat("\n=== SUMMER ΔTmax vs HOBO (2×LAI), by LAI scenario ===\n")
print(as.data.frame(m[, .(scenario, mean_LAI, n, r2=round(r2,3), rmse=round(rmse,3), bias=round(bias,3))]), row.names=FALSE)
fwrite(m, "/home/corroyez/Documents/NC_Full/manuscripts/ch3/tables/Table4e_summer_simple.csv")

lab <- c(SUM_S2_ATBD="S2 ATBD (raw)", SUM_S2_opt="S2 opt (raw)",
         SUM_ALS="LiDAR full", SUM_ALS_dopt="LiDAR d_opt",
         SUM_ATBD_rfull="S2 ATBD ×ratio→full", SUM_ATBD_rdopt="S2 ATBD ×ratio→d_opt",
         SUM_fusion_h="fusion (S2/LiDAR by height)")
m[, lbl := factor(lab[scenario], levels=lab[order(m$rmse)])]
g <- ggplot(m, aes(lbl, rmse, fill=r2)) + geom_col(width=0.7) +
  geom_text(aes(label=sprintf("R²=%.2f\nbias=%+.2f", r2, bias)), vjust=-0.2, size=2.6) +
  scale_fill_viridis_c(option="D", name="R²") +
  labs(title="Summer sub-canopy ΔTmax error by LAI magnitude scenario (Blois, 2×LAI)",
       subtitle="Static canopy, real LAD; MuSICA vs HOBO. Lower RMSE = better.",
       x=NULL, y="RMSE ΔTmax (°C)") +
  theme_minimal(base_size=11) + theme(axis.text.x=element_text(angle=25, hjust=1))
ggsave("/home/corroyez/Documents/NC_Full/manuscripts/ch3/figures/Fig4e_summer_simple.png", g, width=9, height=4.5, dpi=150)
cat("\nDONE\n")
