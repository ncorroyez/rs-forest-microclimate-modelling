# ==============================================================================
# Chapter 3 — (A) per-window paired-bootstrap significance to augment Table 4
#             (B) temporal LAI <-> buffering correlation (annex robustness):
#                 raw vs season-removed vs two-way-demeaned (pure temporal anomaly)
#
# NO re-simulation: reuses cached .nc (validate skips when nc exists).
# Run from z_Example repo root:  Rscript c3_signif_temporal.R
# ==============================================================================
suppressPackageStartupMessages({
  library(terra); library(sf); library(ncdf4); library(lubridate)
  library(dplyr); library(tidyr); library(stringr); library(purrr)
  library(rmusica); library(musica.tools)
})
src <- list.files("R", pattern = "\\.R$", full.names = TRUE)
src <- src[!grepl("/(h1_|lovb_)", src)]
invisible(lapply(src, source))
source("Chapter3_config.R")

prep    <- load_lai_prep(CFG_C3)
df_hobo <- prep$df_plots
ts_list <- prep$ts_by_plot

REF <- "DYN_S2_ANNUAL"
T4_NC  <- file.path(CFG_C3$out_dir, "tables", "Table4_scenarios_MuSICA_windows.csv")
t4_alt <- "/home/corroyez/Documents/NC_Full/output/tables/Table4_scenarios_MuSICA_windows.csv"
t4_ch3 <- "/home/corroyez/Documents/NC_Full/manuscripts/ch3/tables/Table4_scenarios_MuSICA_windows.csv"
# Source-of-truth Table4 (NC_Full); fall back to repo copy if needed.
t4_path <- if (file.exists(t4_alt)) t4_alt else T4_NC
t4 <- read.csv(t4_path, stringsAsFactors = FALSE)
keep <- unique(c(t4$scenario, REF))

all_sc <- make_all_scenarios_c3(ts_list, list_year = CFG_C3$list_year,
                                d_opt = CFG_C3$d_opt_m, mode = CFG_C3$scenarios_mode)
nc_parent <- file.path(CFG_C3$out_dir, "nc")
sc <- Filter(function(s) {
  s$name %in% keep &&
  dir.exists(file.path(nc_parent, s$name)) &&
  length(list.files(file.path(nc_parent, s$name), pattern = "\\.nc$")) > 0
}, all_sc)
cat(sprintf("Scenarios used (%d): %s\n", length(sc),
            paste(vapply(sc, function(s) s$name, ""), collapse = ", ")))

# Windows EXACTLY as in reextract_windows.R (so the new sig column aligns with the
# existing rmse/r2 rows in Table 4).
mkseq <- function(a, b) seq(as.Date(a), as.Date(b), by = "day")
WIN <- list(
  summer   = mkseq("2021-06-01", "2021-09-30"),
  autumn   = mkseq("2021-10-01", "2021-10-31"),
  winter   = c(mkseq("2021-01-01", "2021-03-31"), mkseq("2021-11-01", "2021-12-31")),
  fullyear = mkseq("2021-01-01", "2021-12-31")
)
WIN <- WIN[names(WIN) %in% unique(t4$window)]   # only windows present in Table4

set.seed(7); B <- 5000
boot_ci <- function(x) {
  x <- x[is.finite(x)]
  m <- replicate(B, mean(sample(x, replace = TRUE)))
  quantile(m, c(.025, .975))
}

sig_rows  <- list()
daily_ref <- NULL   # keep fullyear daily (REF) for the temporal annex
for (wn in names(WIN)) {
  ds <- WIN[[wn]]
  dm <- extract_macro_daily(CFG_C3$forcing_file, ds)
  hd <- read_hobo_daily(CFG_C3$hobo_temp_csv, ds, dm, CFG_C3$ids_to_remove)
  val <- validate_scenarios_at_hobos(df_hobo, hd, sc, nc_parent, dm, ds,
                                     CFG_C3$forcing_file, CFG_C3$musica_cmd, FALSE)
  d <- val$daily
  if (wn == "fullyear") daily_ref <- d
  pr <- d %>% group_by(scenario, id_plot) %>%
    summarise(rmse = sqrt(mean((Delta_sim - Delta_obs)^2, na.rm = TRUE)),
              .groups = "drop")
  W <- pr %>% pivot_wider(names_from = scenario, values_from = rmse)
  W <- W[stats::complete.cases(W), ]
  np <- nrow(W)
  for (s in unique(t4$scenario[t4$window == wn])) {
    if (!s %in% names(W)) {
      sig_rows[[length(sig_rows) + 1]] <- data.frame(
        window = wn, scenario = s, n_plots_paired = np,
        dRMSE_vs_DYN_ANNUAL = NA, ci95_lo = NA, ci95_hi = NA, sig = "NA")
      next
    }
    if (s == REF) {
      sig_rows[[length(sig_rows) + 1]] <- data.frame(
        window = wn, scenario = s, n_plots_paired = np,
        dRMSE_vs_DYN_ANNUAL = 0, ci95_lo = 0, ci95_hi = 0, sig = "ref")
      next
    }
    dlt <- W[[s]] - W[[REF]]
    ci  <- boot_ci(dlt)
    sig_rows[[length(sig_rows) + 1]] <- data.frame(
      window = wn, scenario = s, n_plots_paired = np,
      dRMSE_vs_DYN_ANNUAL = round(mean(dlt, na.rm = TRUE), 4),
      ci95_lo = round(ci[1], 4), ci95_hi = round(ci[2], 4),
      sig = if (ci[1] > 0 | ci[2] < 0) "SIG" else "ns")
  }
}
sig_df <- bind_rows(sig_rows)

# ---- Merge into Table 4 and write (both NC_Full locations) -------------------
t4_aug <- t4 %>% left_join(sig_df, by = c("window", "scenario")) %>%
  arrange(window, rmse)
write.csv(t4_aug, t4_alt, row.names = FALSE)
write.csv(t4_aug, t4_ch3, row.names = FALSE)
cat("\n=== Augmented Table 4 (fullyear rows) ===\n")
print(as.data.frame(t4_aug[t4_aug$window == "fullyear",
        c("scenario","rmse","r2","dRMSE_vs_DYN_ANNUAL","ci95_lo","ci95_hi","sig")]),
      row.names = FALSE, digits = 3)

# ============================================================================
# (B) Temporal LAI <-> buffering correlation (annex)
# LAI(t) per plot = cached annual S2 series ts_list$annual[[pid]] (doy, lai).
# Buffering(t)    = observed Delta_obs (ΔTmax) per (id_plot, date), REF daily.
# ============================================================================
ann <- ts_list$annual                                   # named by pid = X{x}_Y{y}
key <- df_hobo %>% mutate(pid = sprintf("X%d_Y%d", round(x), round(y))) %>%
  dplyr::select(id_plot, pid)

obs <- daily_ref %>% filter(scenario == REF) %>%
  dplyr::select(id_plot, date, Delta_obs) %>%
  mutate(doy = yday(as.Date(date))) %>%
  left_join(key, by = "id_plot") %>%
  filter(pid %in% names(ann))

lai_long <- bind_rows(lapply(unique(obs$pid), function(p) {
  s <- ann[[p]]; data.frame(pid = p, doy = s$doy, lai = s$lai)
}))
dat <- obs %>% inner_join(lai_long, by = c("pid", "doy")) %>%
  filter(is.finite(lai), is.finite(Delta_obs))

npl <- length(unique(dat$pid))
cat(sprintf("\n=== Temporal annex: %d plots with annual S2 series, %d plot-days ===\n",
            npl, nrow(dat)))

# 1) Pooled raw temporal correlation (dominated by the shared seasonal cycle)
r_raw <- cor(dat$lai, dat$Delta_obs, use = "complete.obs")

# 2) Per-plot temporal r, averaged (still seasonal within a plot)
per_plot <- dat %>% group_by(pid) %>%
  summarise(r = if (sd(lai) > 0 & sd(Delta_obs) > 0)
    cor(lai, Delta_obs) else NA_real_, nd = n(), .groups = "drop")
r_plot_mean <- mean(per_plot$r, na.rm = TRUE)

# 3a) Season-removed: subtract the cross-plot daily mean (the shared seasonal cycle)
day_mu <- dat %>% group_by(doy) %>%
  summarise(lai_d = mean(lai), del_d = mean(Delta_obs), .groups = "drop")
dat2 <- dat %>% left_join(day_mu, by = "doy") %>%
  mutate(lai_s = lai - lai_d, del_s = Delta_obs - del_d)
r_deseason <- cor(dat2$lai_s, dat2$del_s, use = "complete.obs")

# 3b) Two-way demeaned (remove plot mean = spatial AND day mean = season):
#     pure temporal-anomaly coupling.
gm_l <- mean(dat$lai); gm_d <- mean(dat$Delta_obs)
plot_mu <- dat %>% group_by(pid) %>%
  summarise(lai_p = mean(lai), del_p = mean(Delta_obs), .groups = "drop")
dat3 <- dat %>% left_join(day_mu, by = "doy") %>% left_join(plot_mu, by = "pid") %>%
  mutate(lai_tw = lai - lai_d - lai_p + gm_l,
         del_tw = Delta_obs - del_d - del_p + gm_d)
r_twoway <- cor(dat3$lai_tw, dat3$del_tw, use = "complete.obs")

# Variance accounting: how much LAI temporal variance survives season removal?
v_tot  <- var(dat$lai)
v_des  <- var(dat2$lai_s)
v_tw   <- var(dat3$lai_tw)
frac_season_lai <- 1 - v_des / v_tot     # share of LAI variance that is seasonal

annex <- data.frame(
  quantity = c("pooled_raw_r", "per_plot_mean_r",
               "season_removed_r", "two_way_demeaned_r",
               "LAI_var_share_seasonal", "n_plots", "n_plot_days"),
  value = c(round(r_raw, 3), round(r_plot_mean, 3),
            round(r_deseason, 3), round(r_twoway, 3),
            round(frac_season_lai, 3), npl, nrow(dat)))
annex_csv <- file.path(CFG_C3$out_dir, "tables", "c3_temporal_r_annex.csv")
write.csv(annex, annex_csv, row.names = FALSE)
write.csv(annex, "/home/corroyez/Documents/NC_Full/output/tables/TableA1_temporal_r_annex.csv",
          row.names = FALSE)
write.csv(annex, "/home/corroyez/Documents/NC_Full/manuscripts/ch3/tables/TableA1_temporal_r_annex.csv",
          row.names = FALSE)

cat("\n=== Temporal r annex ===\n")
print(annex, row.names = FALSE)
cat(sprintf(paste0(
  "\nReading: raw temporal cor(LAI, ΔTmax) = %.2f (seasonal co-variation);\n",
  "after removing the shared seasonal cycle it collapses to %.2f, and the pure\n",
  "two-way (plot+day demeaned) temporal-anomaly coupling is %.2f. %.0f%% of the\n",
  "LAI temporal variance is the shared seasonal shape — i.e. the temporal LAI\n",
  "signal is almost entirely season, leaving little independent day-to-day\n",
  "information to drive ΔTmax. Consistent with the scenario bootstrap (RMSE).\n"),
  r_raw, r_deseason, r_twoway, 100 * frac_season_lai))
cat("\nDONE\n")
