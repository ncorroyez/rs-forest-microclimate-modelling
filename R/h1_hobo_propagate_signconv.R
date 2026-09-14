# ==============================================================================
# Propagate new sign convention (x = Tmicro - Tmacro, °C) to :
#   - tab_HOBO_jonckheere_dTmax.csv + fig_HOBO_jonckheere_bins_dTmax.png
#   - tab_HOBO_spearman_daily_block.csv
#   - tab_HOBO_validation_nonparam.csv
# ==============================================================================

suppressMessages({
  library(data.table); library(here); library(cli); library(ggplot2)
})

# ---- JT helpers (same as before) --------------------------------------------
jt_test <- function(x, group, alternative = c("two.sided", "increasing", "decreasing")) {
  alternative <- match.arg(alternative)
  group <- as.factor(group); g_levels <- levels(group); k <- length(g_levels)
  J <- 0
  for (i in seq_len(k - 1L)) for (j in seq.int(i + 1L, k)) {
    x_i <- x[group == g_levels[i]]; x_j <- x[group == g_levels[j]]
    for (a in x_i) J <- J + sum(x_j > a) + 0.5 * sum(x_j == a)
  }
  n   <- length(x); n_g <- tabulate(group, nbins = k)
  mu  <- (n^2 - sum(n_g^2)) / 4
  sigma2 <- (n^2 * (2*n + 3) - sum(n_g^2 * (2*n_g + 3))) / 72
  z <- (J - mu) / sqrt(sigma2)
  p <- switch(alternative,
                "two.sided"  = 2 * (1 - pnorm(abs(z))),
                "increasing" = 1 - pnorm(z),
                "decreasing" = pnorm(z))
  list(statistic = J, z = z, p.value = p, alternative = alternative)
}
jt_boot_ci <- function(x, group, alternative, B = 1000L, seed = 42L) {
  set.seed(seed); n <- length(x); zb <- numeric(B)
  for (b in seq_len(B)) {
    idx <- sample.int(n, n, replace = TRUE)
    zb[b] <- tryCatch(jt_test(x[idx], group[idx], alternative)$z,
                       error = function(e) NA_real_)
  }
  quantile(zb, c(0.025, 0.975), na.rm = TRUE, type = 7)
}

# ---- JT on signed dT --------------------------------------------------------
cli_h1("Propagate sign convention to JT_dTmax")
DT_h <- as.data.table(readRDS(here::here("outputs/lovb/data/DT_HOBO_scalars.rds")))
cg <- readRDS(here::here("outputs/lovb/data/DT_cross_gam_targets.rds"))
hobo_obs <- as.data.table(cg$HOBO)[, .(id_plot = id, dTmax_obs = dTmax_mean)]
DT <- merge(DT_h, hobo_obs, by = "id_plot")
DT[, dT_signed := -dTmax_obs]   # Tmicro - Tmacro

vars  <- c("LAI", "Hmax", "fCover", "LAD_FPC1", "LAD_L2")
trait_cols <- c(LAI = "LAI", Hmax = "Hmax", fCover = "fCover",
                 LAD_FPC1 = "FPC1", LAD_L2 = "LAD_L2")
# Physical : trait ↑  ->  buffering ↑  ->  Tmicro-Tmacro ↓  =>  "decreasing"
alt_lookup <- c(LAI = "decreasing", Hmax = "decreasing",
                 fCover = "decreasing",
                 LAD_FPC1 = "two.sided", LAD_L2 = "two.sided")

rows <- list(); bins_long <- list()
for (v in vars) {
  trait <- DT[[trait_cols[v]]]; y <- DT$dT_signed
  keep <- !is.na(trait) & !is.na(y); trait <- trait[keep]; y <- y[keep]
  n <- length(trait); rnk <- rank(trait, ties.method = "first")
  bin <- cut(rnk, breaks = c(0, n/3, 2*n/3, n),
              include.lowest = TRUE, labels = c("low","mid","high"))
  bin <- factor(bin, levels = c("low","mid","high"))
  res <- jt_test(y, bin, alternative = alt_lookup[v])
  ci  <- jt_boot_ci(y, bin, alt_lookup[v], B = 1000L,
                      seed = 42L + which(vars == v))
  med <- tapply(y, bin, median); iqr <- tapply(y, bin, IQR)
  trend <- if (res$p.value < 0.05) {
    if (res$z > 0) "Monotonic increase" else "Monotonic decrease"
  } else "No monotonic trend"
  rows[[length(rows)+1L]] <- data.table(
    Variable = v, Metric = "dTmax_signed",
    Alternative = alt_lookup[v],
    med_low = round(med["low"], 3),  med_mid = round(med["mid"], 3),
    med_high = round(med["high"], 3),
    iqr_low = round(iqr["low"], 3),  iqr_mid = round(iqr["mid"], 3),
    iqr_high = round(iqr["high"], 3),
    JT_statistic = res$statistic, JT_z = round(res$z, 3),
    JT_pvalue = signif(res$p.value, 4),
    JT_z_ci_lo = round(ci[1], 3), JT_z_ci_hi = round(ci[2], 3),
    Interpretation = trend
  )
  bins_long[[v]] <- data.table(Variable = v, bin = bin, y = y,
                                 p_value = res$p.value, z = res$z,
                                 ci_lo = ci[1], ci_hi = ci[2])
}
tab_jt <- rbindlist(rows)
cli_alert("JT (dTmax_signed = Tmicro - Tmacro) :")
print(tab_jt[, .(Variable, JT_z, JT_z_ci_lo, JT_z_ci_hi,
                   JT_pvalue, Interpretation)])
fwrite(tab_jt, here::here("outputs/lovb/tables/tab_HOBO_jonckheere_dTmax.csv"))
cli_alert_success("Saved tab_HOBO_jonckheere_dTmax.csv")

# Refresh JT bin figure
long <- rbindlist(bins_long)
long[, Variable := factor(Variable, levels = vars)]
ann <- tab_jt[, .(Variable, label = sprintf("JT p = %s\nz = %+.2f [%+.2f ; %+.2f]",
                                                ifelse(JT_pvalue < 1e-3, "<0.001",
                                                        sprintf("%.3f", JT_pvalue)),
                                                JT_z, JT_z_ci_lo, JT_z_ci_hi))]
ann[, Variable := factor(Variable, levels = vars)]
p_jt <- ggplot(long, aes(x = bin, y = y)) +
  geom_hline(yintercept = 0, linetype = "dotted", colour = "grey60") +
  geom_boxplot(fill = "#2C5F2D", alpha = 0.55, colour = "grey20",
                outlier.size = 1.3) +
  geom_jitter(width = 0.12, alpha = 0.55, size = 1, colour = "grey30") +
  geom_text(data = ann, aes(x = -Inf, y = Inf, label = label),
             hjust = -0.05, vjust = 1.3, size = 2.9, fontface = "bold",
             inherit.aes = FALSE) +
  facet_wrap(~ Variable, nrow = 1L, scales = "free_y") +
  labs(
    title    = "Jonckheere-Terpstra trend test (Tmicro,max − Tmacro,max)",
    subtitle = "Tertile binning by trait  |  n = 53 HOBO  |  JT z bootstrap CI95 (1000 reps)",
    x = NULL,
    y = "Observed  T_micro,max − T_macro,max  (°C)",
    caption = paste0("Negative y = canopy buffers (sub-canopy cooler than macro). ",
                       "Alternative = decreasing for LAI/Hmax/fCover (trait ↑  ->  buffering ↑  ->  y ↓) ; ",
                       "two.sided for LAD (FPC1, L2 scalars).")
  ) +
  theme_bw(base_size = 11) +
  theme(strip.background = element_rect(fill = "grey92", colour = NA),
         strip.text       = element_text(face = "bold"),
         plot.title       = element_text(face = "bold"),
         plot.subtitle    = element_text(colour = "grey25"),
         plot.caption     = element_text(size = 8, colour = "grey35", hjust = 0),
         panel.grid.minor = element_blank())
ggsave(here::here("outputs/lovb/figures/hobo/fig_HOBO_LOVB_jonckheere_dTmax.png"),
        p_jt, width = 13, height = 4.3, dpi = 300)
cli_alert_success("Saved fig_HOBO_jonckheere_bins_dTmax.png")

# ---- Daily block-bootstrap with signed dT -----------------------------------
cli_h1("Daily block-bootstrap with signed dT_signed = Tmicro - Tmacro")
DT_obs_d <- readRDS(here::here("outputs/lovb/data/DT_HOBO_daily_obs.rds"))
sim <- as.data.table(readRDS(here::here("outputs/lovb/data/DT_daily_HOBO.rds")))
DT_obs_d[, dT_signed := -dTmax_obs]
bit_map <- c(LAI="0111", Hmax="1011", fCover="1101", LAD="1110", REF="1111")
sim_ref <- sim[bit_code == bit_map[["REF"]],
                 .(id_plot, date, dTmax_REF = Delta_Tmax)]
sim_w <- merge(DT_obs_d[, .(id_plot, date, dT_signed)],
                sim_ref, by = c("id_plot","date"))
vars_v <- c("LAI","Hmax","fCover","LAD")
for (v in vars_v) {
  sim_v <- sim[bit_code == bit_map[[v]],
                .(id_plot, date, dTmax_LOVB = Delta_Tmax)]
  setnames(sim_v, "dTmax_LOVB", paste0("dTmax_LOVB_", v))
  sim_w <- merge(sim_w, sim_v, by = c("id_plot","date"))
}
# Delta_v_LOVB (raw temp diff: sim_LOVB - sim_REF)
for (v in vars_v) {
  sim_w[, paste0("Delta_", v, "_LOVB") := dTmax_REF - get(paste0("dTmax_LOVB_", v))]
}
spear_pooled <- function(x, y) suppressWarnings(cor(x, y, method = "spearman"))
sensors <- unique(sim_w$id_plot)
set.seed(42)
res_daily <- list()
for (v in vars_v) {
  dcol <- paste0("Delta_", v, "_LOVB")
  rho_hat <- spear_pooled(sim_w$dT_signed, sim_w[[dcol]])
  rb <- numeric(1000L)
  for (b in seq_len(1000L)) {
    idx <- sample(sensors, length(sensors), replace = TRUE)
    sub <- sim_w[id_plot %in% idx]
    rb[b] <- spear_pooled(sub$dT_signed, sub[[dcol]])
  }
  ci <- quantile(rb, c(0.025, 0.975), na.rm = TRUE, type = 7)
  res_daily[[v]] <- data.table(Variable = v,
                                  rho_daily = rho_hat,
                                  ci_lo_block = unname(ci[1]),
                                  ci_hi_block = unname(ci[2]),
                                  n_pairs = nrow(sim_w),
                                  n_sensors = length(sensors),
                                  B_boot = 1000L)
}
tab_d <- rbindlist(res_daily)
# Add sensor-level comparison
sp_sensor <- fread(here::here("outputs/lovb/tables/tab_HOBO_spearman_dTmax.csv"))
tab_d <- merge(tab_d,
                 sp_sensor[, .(Variable, rho_sensor = rho_S,
                                 sensor_ci_lo = ci_lo, sensor_ci_hi = ci_hi)],
                 by = "Variable")
tab_d[, Variable := factor(Variable, levels = vars_v)]
setorder(tab_d, Variable)
cli_alert("Daily block-bootstrap (signed dT) :")
print(tab_d[, .(Variable, rho_sensor, rho_daily, ci_lo_block, ci_hi_block)])
fwrite(tab_d, here::here("outputs/lovb/tables/tab_HOBO_spearman_daily_block.csv"))
cli_alert_success("Saved tab_HOBO_spearman_daily_block.csv")

# ---- Synthesis table --------------------------------------------------------
cli_h1("Regenerate tab_HOBO_validation_nonparam.csv")
classify <- function(rho, jt_p, jt_z) {
  bulk_dir <- (jt_z < 0 && rho > 0) || (jt_z > 0 && rho < 0) ||
                (jt_z < 0 && rho < 0 && abs(rho) > 0.3)  # signed-dT consistency
  rho_abs <- abs(rho)
  # For new convention, positive rho and negative JT z = bulk effect
  if (rho > 0.3 && jt_p < 0.05) "BULK CONFIRMED"
  else if ((rho_abs >= 0.2 && rho_abs <= 0.3) ||
            (jt_p >= 0.05 && jt_p <= 0.10))  "BULK MARGINAL"
  else if (rho_abs < 0.2 && jt_p > 0.10)     "NULL (H2 OK)"
  else                                         "AMBIGUOUS"
}
sp_slope <- fread(here::here("outputs/lovb/tables/tab_HOBO_spearman.csv"))
sp_dTmax <- fread(here::here("outputs/lovb/tables/tab_HOBO_spearman_dTmax.csv"))
jt_slope <- fread(here::here("outputs/lovb/tables/tab_HOBO_jonckheere_slope.csv"))
jt_dTmax <- fread(here::here("outputs/lovb/tables/tab_HOBO_jonckheere_dTmax.csv"))
jt_slope_main <- jt_slope[Variable %in% c("LAI","Hmax","fCover","LAD_FPC1")]
jt_dTmax_main <- jt_dTmax[Variable %in% c("LAI","Hmax","fCover","LAD_FPC1")]
jt_slope_main[, Variable := gsub("_FPC1$", "", Variable)]
jt_dTmax_main[, Variable := gsub("_FPC1$", "", Variable)]
build_one <- function(sp, jt, metric) {
  DT_m <- merge(sp[, .(Variable, rho_S, p_value, ci_lo, ci_hi)],
                  jt[, .(Variable, JT_z, JT_pvalue, JT_z_ci_lo, JT_z_ci_hi)],
                  by = "Variable")
  DT_m[, x_metric := metric]
  DT_m[, Interpretation := mapply(classify, rho_S, JT_pvalue, JT_z)]
  setcolorder(DT_m, c("Variable", "x_metric"))
  DT_m
}
tab_synth <- rbind(build_one(sp_slope, jt_slope_main, "slope_obs"),
                     build_one(sp_dTmax, jt_dTmax_main, "Tmicro-Tmacro"))
tab_synth[, Variable := factor(Variable, levels = c("LAI","Hmax","fCover","LAD"))]
setorder(tab_synth, x_metric, Variable)
fwrite(tab_synth,
        here::here("outputs/lovb/tables/tab_HOBO_validation_nonparam.csv"))
cli_alert("Synthesis table :")
print(tab_synth[, .(Variable, x_metric, rho_S, JT_z, JT_pvalue, Interpretation)])
cli_alert_success("Saved tab_HOBO_validation_nonparam.csv")
