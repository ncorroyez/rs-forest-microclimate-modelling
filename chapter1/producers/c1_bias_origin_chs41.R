# ==============================================================================
# #1 — Bias-origin test: is the +0.9 °C model warm bias inherited from the ERA5
# forcing, or generated inside MuSICA's canopy scheme?
#
# Logic (advisor-validated). The chapter's ΔTmax uses Tmax_macro = ERA5 for BOTH
# model and obs, so the macro cancels in the bias:
#     bias = ΔTmax_model − ΔTmax_obs = max(T_micro_MuSICA) − max(T_HOBO) = +0.9 °C
# The +0.9 is "MuSICA's 1 m air runs warm". MuSICA is driven by ERA5 aloft, so:
#   * if ERA5 daily-max runs WARM vs an independent open-field station (CHS41
#     RENECOFOR), the warm forcing propagates down → forcing-origin CONFIRMED;
#   * if ERA5 ≈ CHS41, the forcing is exonerated → the +0.9 is canopy-scheme
#     under-buffering → §4.2 "most likely the forcing" must be RETRACTED.
# CHS41 = RENECOFOR sessile-oak met station, dept 41 (Blois), open-field 1.5 m air.
#
#   Rscript c1_bias_origin_chs41.R
# Out: tab_bias_origin_chs41.csv + FigSh_bias_origin_chs41.png
# ==============================================================================
suppressPackageStartupMessages({
  library(ncdf4); library(data.table); library(ggplot2); library(patchwork)
  library(lubridate); library(dplyr)
})
source("R/config.R")   # CFG (forcing_file, date_seq)
source("R/io.R")       # extract_macro_daily(), force_utc_nc()

WIN   <- CFG$date_seq                       # 2021-06-01 .. 2021-09-30
MODEL_BIAS <- 0.9                           # +0.9 °C reported model warm bias (MuSICA_micro − HOBO)

# ---- ERA5 (the forcing aloft) : daily max, exactly as the chapter computes it -
era5 <- as.data.table(extract_macro_daily(CFG$forcing_file, WIN))[, .(date, ERA5_max = Tmax_macro)]

# ---- CHS41 (independent open-field RENECOFOR station) : hourly → daily max -----
raw <- fread("MetHor2021.txt", sep = ";", header = TRUE, encoding = "Latin-1")
setnames(raw, 1, "code"); tcol <- grep("instantan", names(raw), value = TRUE)[1]
chs <- raw[code == "CHS 41", .(code, Date, Heure = `Heure (TU)`, Tair = as.numeric(get(tcol)))]
chs[, date := as.Date(Date, format = "%d/%m/%Y")]
chs <- chs[date %in% WIN & is.finite(Tair)]
chs_d <- chs[, .(CHS41_max = max(Tair, na.rm = TRUE),
                 n_hr      = .N), by = date][n_hr >= 20]   # require near-complete days

# ---- merge & decompose ---------------------------------------------------------
D <- merge(era5, chs_d, by = "date")[is.finite(ERA5_max) & is.finite(CHS41_max)]
D[, dmax := ERA5_max - CHS41_max]                          # forcing − truth at daily max

fb <- function(x) c(mean = mean(x), sd = sd(x), med = median(x),
                    lo = quantile(x, .025, names = FALSE), hi = quantile(x, .975, names = FALSE))
b_max <- fb(D$dmax)

# flat offset vs conditional? regress daily-max bias on the day's CHS41 max
fit <- lm(dmax ~ CHS41_max, data = D)
slope <- coef(fit)[2]; slope_p <- summary(fit)$coefficients[2, 4]

# share of the +0.9 model bias that an ERA5 warm offset could account for
share <- b_max["mean"] / MODEL_BIAS

res <- data.table(
  metric          = "daily-max (ΔTmax-relevant)",
  n_days          = nrow(D),
  ERA5_minus_CHS41_mean = round(b_max["mean"], 3),
  sd              = round(b_max["sd"], 3),
  ci_lo           = round(b_max["lo"], 3),
  ci_hi           = round(b_max["hi"], 3),
  model_bias_ref  = MODEL_BIAS,
  share_of_model_bias = round(share, 3),
  slope_vs_CHS41max = round(slope, 4),
  slope_p         = round(slope_p, 4))
fwrite(res, "out_files/Chapter1/tables/tab_bias_origin_chs41.csv")

cat("\n=== #1 BIAS-ORIGIN TEST (ERA5 vs CHS41, daily max, summer 2021) ===\n")
print(res)
cat(sprintf("\nModel warm bias to explain (MuSICA_micro − HOBO): +%.2f °C\n", MODEL_BIAS))
cat(sprintf("ERA5 − CHS41 at daily max: %+.2f °C  [95%% %.2f, %.2f], n=%d days\n",
            b_max["mean"], b_max["lo"], b_max["hi"], nrow(D)))
cat(sprintf("=> ERA5 warm offset accounts for ~%.0f%% of the +%.1f model bias.\n",
            100 * share, MODEL_BIAS))
verdict <- if (b_max["mean"] >= 0.6) "FORCING-ORIGIN SUPPORTED (ERA5 runs warm)  -> §4.2 stands" else
           if (b_max["mean"] <= 0.3) "FORCING EXONERATED (ERA5 ~ CHS41)  -> §4.2 must be RETRACTED (canopy-scheme under-buffering)" else
           "PARTIAL: ERA5 warm offset explains only part -> §4.2 must be SOFTENED (forcing + canopy both)"
cat("VERDICT:", verdict, "\n")

# ---- figure : time series + scatter -------------------------------------------
DL <- melt(D[, .(date, ERA5 = ERA5_max, CHS41 = CHS41_max)], id.vars = "date")
p1 <- ggplot(DL, aes(date, value, colour = variable)) +
  geom_line(linewidth = .5) +
  scale_colour_manual(values = c(ERA5 = "#D6604D", CHS41 = "#4393C3"), name = NULL) +
  labs(x = NULL, y = "Daily max air T (°C)",
       title = "Forcing aloft vs independent station (daily max, summer 2021)",
       subtitle = sprintf("ERA5 − CHS41 = %+.2f °C  [95%% %.2f, %.2f] ; model warm bias to explain = +%.1f °C",
                          b_max["mean"], b_max["lo"], b_max["hi"], MODEL_BIAS)) +
  theme_bw(base_size = 11) + theme(legend.position = "bottom",
       plot.subtitle = element_text(size = 8, colour = "grey35"))
rng <- range(c(D$ERA5_max, D$CHS41_max))
p2 <- ggplot(D, aes(CHS41_max, ERA5_max)) +
  geom_abline(slope = 1, intercept = 0, linetype = 2, colour = "grey50") +
  geom_point(alpha = .6, colour = "#762A83") +
  coord_equal(xlim = rng, ylim = rng) +
  labs(x = "CHS41 daily max (°C)", y = "ERA5 daily max (°C)",
       title = "ERA5 vs CHS41 (1:1 dashed)",
       subtitle = sprintf("slope on CHS41max = %.2f (p=%.3f): %s offset",
                          slope, slope_p, ifelse(slope_p < .05, "conditional", "flat"))) +
  theme_bw(base_size = 11) + theme(plot.subtitle = element_text(size = 8, colour = "grey35"),
       aspect.ratio = 1)
ggsave("out_files/Chapter1/figures/FigSh_bias_origin_chs41.png",
       p1 + p2 + plot_layout(widths = c(1.4, 1)), width = 11, height = 4.6, dpi = 200, bg = "white")
cat("DONE -> tab_bias_origin_chs41.csv + FigSh_bias_origin_chs41.png\n")
