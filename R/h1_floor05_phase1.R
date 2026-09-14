# ==============================================================================
# Phase 1 — Prep floored inputs.
# Floor real fCover to max(real, 0.5) in df_sample and df_hobo_inputs.
# Recompute C1 archetype centroid with floored values.
# Identify affected (plot, coalition) tuples for selective re-run.
# ==============================================================================

suppressMessages({
  library(data.table); library(here); library(cli)
})

cli_h1("Phase 1 — Floor inputs + identify affected sims")

# --- 1. Floor cLHS sample -----------------------------------------------------
df_orig <- as.data.table(readRDS(here::here("out_files/Sensitivity_Analysis/clhs_sample.rds")))
if ("Archetype" %in% names(df_orig) && !"Cluster" %in% names(df_orig))
  df_orig[, Cluster := Archetype]

df_floor <- copy(df_orig)
df_floor[fCover < 0.5, fCover := 0.5]
cli_alert("cLHS : {sum(df_orig$fCover < 0.5)} plots floored to 0.5")
cat("Per cluster summary :\n")
print(df_floor[, .(n_floored = sum(df_orig$fCover < 0.5 & Cluster == .BY$Cluster),
                     fCover_mean_new = mean(fCover),
                     fCover_mean_old = mean(df_orig$fCover[df_orig$Cluster == .BY$Cluster])),
                by = Cluster][order(Cluster)])

out_dir <- here::here("out_files/Sensitivity_Analysis")
saveRDS(df_floor, file.path(out_dir, "clhs_sample_floor05.rds"))
cli_alert_success("Saved clhs_sample_floor05.rds")

# --- 2. Identify affected cLHS sims --------------------------------------------
# fCover bit=1 (real value used) -> these 8 coalitions need re-run for affected plots
COAL_FCOVER_REAL <- c("0010","0011","0110","0111","1010","1011","1110","1111")
# Filename → scenario suffix mapping (from R/lovb_01_load.R conventions)
bit_to_suffix <- function(b) {
  parts <- strsplit(b, "")[[1]]
  paste(ifelse(parts[1]=="1","r","m"),
         ifelse(parts[2]=="1","r","m"),
         ifelse(parts[3]=="1","r","a"),
         ifelse(parts[4]=="1","r","u"), sep = "_")
}
scenarios_affected <- vapply(COAL_FCOVER_REAL,
                               function(b) paste0("H1F_", bit_to_suffix(b)),
                               character(1))
affected_plots <- df_orig[fCover < 0.5, .(x, y)]
cli_alert("cLHS re-run scope : {nrow(affected_plots)} plots x {length(scenarios_affected)} coalitions = {nrow(affected_plots) * length(scenarios_affected)} sims")
saveRDS(list(plots = affected_plots, scenarios = scenarios_affected,
              coalitions = COAL_FCOVER_REAL),
         file.path(out_dir, "floor05_affected_clhs.rds"))

# --- 3. Floor HOBO sensors -----------------------------------------------------
DT_h <- as.data.table(readRDS(here::here("outputs/lovb/data/DT_HOBO_scalars.rds")))
n_hobo_floor <- sum(DT_h$fCover < 0.5)
cli_alert("HOBO : {n_hobo_floor}/{nrow(DT_h)} sensors with fCover < 0.5")
DT_h_floor <- copy(DT_h)
DT_h_floor[fCover < 0.5, fCover := 0.5]
saveRDS(DT_h_floor, here::here("outputs/lovb/data/DT_HOBO_scalars_floor05.rds"))
cli_alert_success("Saved DT_HOBO_scalars_floor05.rds")

# Affected HOBO sensors and coalitions
affected_hobo <- DT_h[fCover < 0.5, .(id_plot, fCover)]
# Coalitions where fCover is real (in HOBO 10-coalition set)
COAL_FCOVER_REAL_HOBO <- c("0010", "0111", "1011", "1110", "1111")
cli_alert("HOBO re-run scope : {nrow(affected_hobo)} sensors x {length(COAL_FCOVER_REAL_HOBO)} coalitions = {nrow(affected_hobo) * length(COAL_FCOVER_REAL_HOBO)} sims")
saveRDS(list(sensors = affected_hobo, coalitions = COAL_FCOVER_REAL_HOBO),
         here::here("outputs/lovb/data/floor05_affected_hobo.rds"))

# --- 4. C1 archetype centroid (recomputed) -------------------------------------
c1_mean <- df_floor[Cluster == 1, .(LAI = mean(LAI), Hmax = mean(Hmax),
                                       fCover = mean(fCover))]
c1_orig <- df_orig[Cluster == 1, .(LAI = mean(LAI), Hmax = mean(Hmax),
                                      fCover = mean(fCover))]
cli_alert("C1 archetype centroid : fCover {round(c1_orig$fCover, 3)} -> {round(c1_mean$fCover, 3)}")
cli_alert("C1 archetype re-run scope : 16 coalitions (centroid changed)")
saveRDS(c1_mean, here::here("outputs/lovb/data/floor05_C1_centroid.rds"))

cli_alert_success("Phase 1 done. Total sims to re-run : {nrow(affected_plots)*length(scenarios_affected) + 16 + nrow(affected_hobo)*length(COAL_FCOVER_REAL_HOBO)}")
