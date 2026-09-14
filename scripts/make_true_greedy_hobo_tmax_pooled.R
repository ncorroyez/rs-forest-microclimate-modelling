# ==============================================================================
# TRUE greedy forward — DAILY POOLED variant
#   - Each panel uses ALL (plot, day) pairs (~6500 points instead of 53)
#   - r computed pooled across plots AND days
#   - Same greedy logic : at each step add the variable maximising r on HOBO
#
# Complements scripts/make_true_greedy_hobo_tmax.R (per-plot mean variant).
# ==============================================================================

suppressMessages({
  library(data.table); library(here); library(cli); library(ncdf4)
  library(tidyverse); library(lubridate); library(musica.tools)
})

source(here::here("R/config.R"))
source(here::here("R/io.R"))
source(here::here("R/musica.R"))
source(here::here("R/h1_shapley_archetypes.R"))

OUT <- here::here("outputs/figs_MEB2026_final")
df_macro <- extract_macro_daily(CFG$forcing_file, CFG$date_seq)

# ---- HOBO observed daily ΔTmax (one row per (plot, date)) ------------------
hobo_obs_daily <- as.data.table(read_hobo_daily(CFG$hobo_temp_csv, CFG$date_seq,
                                                 df_macro, CFG$ids_to_remove))[
  , .(id_plot, date, Delta_obs)]
cli_alert("HOBO obs daily rows : {nrow(hobo_obs_daily)} | {length(unique(hobo_obs_daily$id_plot))} plots")

# ---- Sim daily ΔTmax per (plot, date) for one coalition --------------------
collect_step_daily <- function(bit) {
  v10 <- file.path("out_files/musica_hobo_v10_fcovmean", bit)
  v9  <- file.path("out_files/musica_hobo_v9", bit)
  subdir <- if (dir.exists(v10) && length(list.files(v10, "\\.nc$")) >= 53) v10 else v9
  nc_files <- list.files(subdir, pattern = "\\.nc$", full.names = TRUE)
  rows <- list()
  for (f in nc_files) {
    id <- sub("musica_out_HOBO_(.+)\\.nc$", "\\1", basename(f))
    res <- tryCatch(extract_deltatmax_one(f, df_macro, CFG$date_seq,
                                          z_target = CFG$tair_target_height),
                    error = function(e) NULL)
    if (is.null(res) || nrow(res) == 0) next
    d <- as.data.table(res)[, .(id_plot = id, date, Delta_sim = Delta_Tmax)]
    rows[[id]] <- d
  }
  rbindlist(rows)
}

# ---- Fit per coalition : pool all (plot, day) pairs ------------------------
fit_bit_daily <- function(bit) {
  s <- collect_step_daily(bit)
  if (nrow(s) == 0) return(list(r = NA, RMSE = NA, MAE = NA, bias = NA, n = 0))
  m <- merge(s, hobo_obs_daily, by = c("id_plot", "date"))
  ok <- !is.na(m$Delta_sim) & !is.na(m$Delta_obs)
  x <- m$Delta_obs[ok]; y <- m$Delta_sim[ok]
  list(r = suppressWarnings(cor(x, y)),
       RMSE = sqrt(mean((y - x)^2)),
       MAE = mean(abs(y - x)),
       bias = mean(y - x),
       n = length(x),
       df_long = m)
}

cli_h1("Fitting all 16 coalitions (DAILY POOLED ΔTmax)")
bit_table <- data.table(bit = unname(.ARCH_COAL_MAP))
fits <- list()
for (b in bit_table$bit) fits[[b]] <- fit_bit_daily(b)
all_fits <- data.table(bit = names(fits),
                        r    = sapply(fits, `[[`, "r"),
                        RMSE = sapply(fits, `[[`, "RMSE"),
                        MAE  = sapply(fits, `[[`, "MAE"),
                        bias = sapply(fits, `[[`, "bias"),
                        n    = sapply(fits, `[[`, "n"))
print(all_fits)
fwrite(all_fits, file.path(OUT, "tab_v10_HOBO_all16_fits_Tmax_daily_pooled.csv"))

# ---- Greedy forward ---------------------------------------------------------
cli_h1("TRUE greedy forward selection (HOBO Tmax, daily pooled)")

VARS <- c("LAI", "Hmax", "fCover", "LAD")
bit_from_set <- function(set) paste(ifelse(VARS %in% set, "1", "0"), collapse = "")

active <- character(0)
path <- list()

# Step 1
b0 <- bit_from_set(active)
f0 <- fits[[b0]]
path[[1]] <- data.table(step = 1, added = "Baseline", active_set = "",
                        bit = b0, r = f0$r, RMSE = f0$RMSE, MAE = f0$MAE,
                        bias = f0$bias, n = f0$n)
for (k in 2:5) {
  remaining <- setdiff(VARS, active)
  candidates <- sapply(remaining, function(v) bit_from_set(c(active, v)))
  cand_r <- sapply(candidates, function(b) fits[[b]]$r)
  best_i <- which.max(cand_r)
  best_v <- remaining[best_i]
  active <- c(active, best_v)
  bit <- candidates[best_i]
  path[[k]] <- data.table(step = k, added = best_v,
                          active_set = paste(active, collapse = "+"),
                          bit = bit,
                          r = cand_r[best_i],
                          RMSE = fits[[bit]]$RMSE,
                          MAE  = fits[[bit]]$MAE,
                          bias = fits[[bit]]$bias,
                          n    = fits[[bit]]$n)
  if (length(active) == length(VARS)) break
}
greedy <- rbindlist(path)
print(greedy)
fwrite(greedy, file.path(OUT, "tab_v10_HOBO_true_greedy_Tmax_daily_pooled.csv"))

# ---- Figure -----------------------------------------------------------------
cli_h1("Figure greedy forward (daily pooled)")

labels_path <- sapply(seq_len(nrow(greedy)), function(k) {
  if (k == 1) "Baseline"
  else if (k == nrow(greedy)) sprintf("REF (+%s)", greedy$added[k])
  else paste("+", paste(strsplit(greedy$active_set[k], "\\+")[[1]],
                         collapse = " + "))
})
labels_path <- gsub("Hmax", "H_max", labels_path)

# Pre-extract sim data once per bit
long <- list()
for (k in seq_len(nrow(greedy))) {
  bit <- greedy$bit[k]
  s <- collect_step_daily(bit)
  if (nrow(s) == 0) next
  m <- merge(s, hobo_obs_daily, by = c("id_plot", "date"))
  m[, step := k]; m[, label := labels_path[k]]
  long[[k]] <- m
}
long <- rbindlist(long, fill = TRUE)
long$label <- factor(as.character(long$label), levels = labels_path)

ann <- greedy[, .(label = labels_path,
                   lab_r    = sprintf("italic(r)==%+.2f", r),
                   lab_rmse = sprintf("RMSE==%.2f * ' °C'", RMSE),
                   lab_mae  = sprintf("MAE==%.2f * ' °C'", MAE),
                   lab_n    = sprintf("n==%d", n))]
ann$label <- factor(ann$label, levels = labels_path)
xy_lim <- range(c(long$Delta_sim, long$Delta_obs), na.rm = TRUE)
xy_pad <- diff(xy_lim) * 0.05; xy_lim <- xy_lim + c(-xy_pad, xy_pad)

p <- ggplot(long, aes(x = Delta_obs, y = Delta_sim)) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", colour = "grey55") +
  geom_hline(yintercept = 0, linetype = "dotted", colour = "grey80") +
  geom_vline(xintercept = 0, linetype = "dotted", colour = "grey80") +
  geom_hex(bins = 60, alpha = 0.85) +
  scale_fill_viridis_c(option = "viridis", trans = "log10",
                        name = "n (log)", labels = scales::label_log()) +
  geom_smooth(method = "lm", formula = y ~ x, se = FALSE,
              colour = "#FFB400", linewidth = 0.9) +
  geom_text(data = ann, aes(x = -Inf, y = Inf, label = lab_r),
            parse = TRUE, hjust = -0.08, vjust = 1.4, size = 5,
            inherit.aes = FALSE, colour = "grey15") +
  geom_text(data = ann, aes(x = -Inf, y = Inf, label = lab_rmse),
            parse = TRUE, hjust = -0.08, vjust = 2.8, size = 5,
            inherit.aes = FALSE, colour = "grey15") +
  geom_text(data = ann, aes(x = -Inf, y = Inf, label = lab_mae),
            parse = TRUE, hjust = -0.08, vjust = 4.2, size = 5,
            inherit.aes = FALSE, colour = "grey15") +
  geom_text(data = ann, aes(x = -Inf, y = Inf, label = lab_n),
            parse = TRUE, hjust = -0.08, vjust = 5.6, size = 4,
            inherit.aes = FALSE, colour = "grey40") +
  facet_wrap(~ label, nrow = 1) +
  coord_fixed(xlim = xy_lim, ylim = xy_lim) +
  labs(title = "V10 — TRUE greedy forward selection on HOBO ΔTmax (DAILY POOLED)",
       x = expression("Field " ~ Delta * T[max] ~ "(°C)"),
       y = expression("Simulated " ~ Delta * T[max] ~ "(°C)")) +
  theme_bw(base_size = 14) +
  theme(strip.background = element_rect(fill = "grey92", colour = NA),
        strip.text = element_text(face = "bold", size = 13),
        plot.title = element_text(face = "bold", size = 13),
        panel.grid.minor = element_blank(),
        legend.position = "right")
ggsave(file.path(OUT, "fig_forward_HOBO_v10_true_greedy_Tmax_daily_pooled.png"),
        p, width = 20, height = 5.5, dpi = 300)
cli_alert_success("Saved fig_forward_HOBO_v10_true_greedy_Tmax_daily_pooled.png")
