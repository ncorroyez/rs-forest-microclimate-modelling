# ==============================================================================
# TRUE greedy forward feature selection — answer to JB Féret.
#
# At each step :
#   - Take the *current* active set S
#   - For each variable v ∉ S, fit Tmax_sim_per_plot ~ Tmax_obs_per_plot
#     using simulated NCs at coalition = S ∪ {v}
#   - Pick the v that maximises validation r (or equivalently min RMSE)
#   - Add v to S
# Repeat until all 4 are in.
#
# Uses the COMPLETE 16-coalition HOBO factorial (V10 baseline where fCover=0).
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

# Helper : per-plot Tmax mean from V10 (fallback V9)
collect_step <- function(bit) {
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
    rows[[id]] <- data.table(id_plot = id,
                              sim_metric = mean(res$Delta_Tmax, na.rm = TRUE))
  }
  rbindlist(rows)
}

# HOBO observed daily ΔTmax mean per plot
hobo_obs <- as.data.table(read_hobo_daily(CFG$hobo_temp_csv, CFG$date_seq,
                                           df_macro, CFG$ids_to_remove))[
  , .(obs_metric = mean(Delta_obs, na.rm = TRUE)), by = id_plot]

# Fit one coalition vs HOBO
fit_bit <- function(bit) {
  d <- collect_step(bit)
  if (nrow(d) == 0) return(list(r = NA, RMSE = NA, MAE = NA, bias = NA, n = 0))
  d <- merge(d, hobo_obs, by = "id_plot")
  ok <- !is.na(d$sim_metric) & !is.na(d$obs_metric)
  x <- d$obs_metric[ok]; y <- d$sim_metric[ok]
  list(r = suppressWarnings(cor(x, y)),
       RMSE = sqrt(mean((y - x)^2)),
       MAE = mean(abs(y - x)),
       bias = mean(y - x),
       n = length(x))
}

# Build all 16 fits once
cli_h1("Fitting all 16 HOBO coalitions vs Tmax_obs")
bit_table <- data.table(bit = unname(.ARCH_COAL_MAP))
bit_table[, c("LAI", "Hmax", "fCover", "LAD") := tstrsplit(bit, "", fixed = TRUE)]
bit_table[, c("LAI", "Hmax", "fCover", "LAD") := lapply(.SD, as.integer),
          .SDcols = c("LAI", "Hmax", "fCover", "LAD")]
fits <- list()
for (b in bit_table$bit) fits[[b]] <- fit_bit(b)
all_fits <- data.table(bit = names(fits),
                       r    = sapply(fits, `[[`, "r"),
                       RMSE = sapply(fits, `[[`, "RMSE"),
                       MAE  = sapply(fits, `[[`, "MAE"),
                       bias = sapply(fits, `[[`, "bias"),
                       n    = sapply(fits, `[[`, "n"))
all_fits <- merge(all_fits, bit_table, by = "bit")
all_fits[, k := LAI + Hmax + fCover + LAD]
setorder(all_fits, k, -r)
print(all_fits)
fwrite(all_fits, file.path(OUT, "tab_v10_HOBO_all16_fits_Tmax.csv"))

# ============================================================================
# Greedy forward : at each step, add variable maximising r vs HOBO
# ============================================================================
cli_h1("TRUE greedy forward selection (HOBO Tmax)")

VARS <- c("LAI", "Hmax", "fCover", "LAD")
bit_from_set <- function(set) {
  paste(ifelse(VARS %in% set, "1", "0"), collapse = "")
}

active <- character(0)
path <- list()

# Step 1 : Baseline (bit 0000)
b0 <- bit_from_set(active)
f0 <- fits[[b0]]
path[[1]] <- data.table(step = 1, added = "Baseline",
                        active_set = "",
                        bit = b0,
                        r = f0$r, RMSE = f0$RMSE, MAE = f0$MAE, bias = f0$bias)

# Greedy iterations
for (k in 2:5) {
  if (length(active) == length(VARS)) break
  remaining <- setdiff(VARS, active)
  candidates <- lapply(remaining, function(v) c(active, v))
  cand_bits <- sapply(candidates, bit_from_set)
  cand_r    <- sapply(cand_bits, function(b) fits[[b]]$r)
  cand_rmse <- sapply(cand_bits, function(b) fits[[b]]$RMSE)
  best_i <- which.max(cand_r)
  best_v <- remaining[best_i]
  active <- c(active, best_v)
  path[[k]] <- data.table(step = k,
                          added = best_v,
                          active_set = paste(active, collapse = "+"),
                          bit = cand_bits[best_i],
                          r = cand_r[best_i], RMSE = cand_rmse[best_i],
                          MAE = fits[[cand_bits[best_i]]]$MAE,
                          bias = fits[[cand_bits[best_i]]]$bias)
}
greedy <- rbindlist(path)
print(greedy)
fwrite(greedy, file.path(OUT, "tab_v10_HOBO_true_greedy_Tmax.csv"))

# ============================================================================
# Figure : 5 panels Baseline → +v1 → +v1+v2 → +v1+v2+v3 → REF
# ============================================================================
cli_h1("Figure greedy forward HOBO Tmax")

# Build the long data for figure
labels_path <- sapply(seq_len(nrow(greedy)), function(k) {
  if (k == 1) "Baseline"
  else if (k == nrow(greedy)) sprintf("REF (+%s)", greedy$added[k])
  else paste("+", paste(strsplit(greedy$active_set[k], "\\+")[[1]],
                         collapse = " + "))
})
# Make labels short by replacing Hmax → H_max
labels_path <- gsub("Hmax", "H_max", labels_path)

long <- list()
for (k in seq_len(nrow(greedy))) {
  bit <- greedy$bit[k]
  d <- collect_step(bit)
  if (nrow(d) == 0) next
  d <- merge(d, hobo_obs, by = "id_plot")
  d[, step := k]; d[, label := labels_path[k]]
  long[[k]] <- d
}
long <- rbindlist(long, fill = TRUE)
long$label <- factor(as.character(long$label), levels = labels_path)

ann <- greedy[, .(label = labels_path,
                   lab_r    = sprintf("italic(r)==%+.2f", r),
                   lab_rmse = sprintf("RMSE==%.2f * ' °C'", RMSE),
                   lab_mae  = sprintf("MAE==%.2f * ' °C'", MAE))]
ann$label <- factor(ann$label, levels = labels_path)
xy_lim <- range(c(long$sim_metric, long$obs_metric), na.rm = TRUE)
xy_pad <- diff(xy_lim) * 0.06; xy_lim <- xy_lim + c(-xy_pad, xy_pad)

p <- ggplot(long, aes(x = obs_metric, y = sim_metric)) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", colour = "grey55") +
  geom_hline(yintercept = 0, linetype = "dotted", colour = "grey80") +
  geom_vline(xintercept = 0, linetype = "dotted", colour = "grey80") +
  geom_point(colour = "#2C5F2D", alpha = 0.75, size = 2.2) +
  geom_smooth(method = "lm", formula = y ~ x, se = TRUE,
              colour = "#FFB400", fill = "#FFB400",
              linewidth = 0.9, alpha = 0.22) +
  geom_text(data = ann, aes(x = -Inf, y = Inf, label = lab_r),
            parse = TRUE, hjust = -0.08, vjust = 1.4, size = 5.5,
            inherit.aes = FALSE, colour = "grey15") +
  geom_text(data = ann, aes(x = -Inf, y = Inf, label = lab_rmse),
            parse = TRUE, hjust = -0.08, vjust = 2.8, size = 5.5,
            inherit.aes = FALSE, colour = "grey15") +
  geom_text(data = ann, aes(x = -Inf, y = Inf, label = lab_mae),
            parse = TRUE, hjust = -0.08, vjust = 4.2, size = 5.5,
            inherit.aes = FALSE, colour = "grey15") +
  facet_wrap(~ label, nrow = 1) +
  coord_cartesian(xlim = xy_lim, ylim = xy_lim) +
  labs(title = "V10 — TRUE greedy forward selection on HOBO ΔTmax (each step: best remaining variable)",
       x = expression("Field " ~ Delta * T[max] ~ "(°C)"),
       y = expression("Simulated " ~ Delta * T[max] ~ "(°C)")) +
  theme_bw(base_size = 14) +
  theme(strip.background = element_rect(fill = "grey92", colour = NA),
        strip.text = element_text(face = "bold", size = 13),
        plot.title = element_text(face = "bold", size = 13),
        panel.grid.minor = element_blank())
ggsave(file.path(OUT, "fig_forward_HOBO_v10_true_greedy_Tmax.png"),
        p, width = 18, height = 5, dpi = 300)
cli_alert_success("Saved fig_forward_HOBO_v10_true_greedy_Tmax.png")

cli_h1("Done.")
