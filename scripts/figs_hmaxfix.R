# ==============================================================================
# Regenerate the MEB figures with the Hmax-fix data :
#   - fig_heatmap_LOO_profiles_gradient.png  (LOO × 4 profiles, with Avg col)
#   - fig_forward_selection_HOBO.png         (5 panels, explicit labels)
# ==============================================================================

suppressMessages({
  library(data.table); library(here); library(cli); library(ncdf4); library(ggplot2)
  library(sf); library(terra); library(dplyr); library(tidyr); library(tibble)
  library(musica.tools); library(rmusica); library(lubridate)
})

source(here::here("R/config.R"))
source(here::here("R/io.R"))
source(here::here("R/musica.R"))
source(here::here("R/lad.R"))
source(here::here("R/validation.R"))
source(here::here("R/h1_shapley_archetypes.R"))

OUT <- here::here("outputs/figs_MEB2026_final")
VARS <- c("LAI","Hmax","fCover","LAD")

# ----------------------------------------------------------------------------
# 1) Heatmap LOO × 4 profiles
# ----------------------------------------------------------------------------
cli_h1("Heatmap LOO × profiles (Hmax fix)")

# Load archetype Tmax data from Phase 2 output
ah <- readRDS(here::here("outputs/hmaxfix_archetype_heatmap_data.rds"))
long <- copy(ah$delta)   # variable, Cluster, delta

# Convert to the figure-friendly format (with Profile P1-P4)
# Cluster = RAW k-means code; map via cluster_relabel (raw 1->P4, 2->P2, 3->P1, 4->P3).
# Direct P%d was WRONG (scrambled archetype labels).
.RELAB <- c("1" = "P4", "2" = "P2", "3" = "P1", "4" = "P3")
long[, Profile := unname(.RELAB[as.character(Cluster)])]
# Add Avg column
avg_rows <- long[, .(Profile = "Avg", delta = mean(delta, na.rm = TRUE)),
                    by = variable]
long_full <- rbind(long[, .(variable, Profile, delta)], avg_rows)
long_full[, Profile := factor(Profile, levels = c(paste0("P", 1:4), "Avg"))]
long_full[, variable := factor(variable, levels = VARS)]
long_full[, var_lab := ifelse(variable == "Hmax", "H[max]",
                                  as.character(variable))]
long_full[, var_lab := factor(var_lab, levels = c("LAI","H[max]","fCover","LAD"))]

max_abs <- max(abs(long_full$delta), na.rm = TRUE) * 1.05
sep_x <- 4.5

p_hm <- ggplot(long_full, aes(x = Profile, y = var_lab, fill = delta)) +
  geom_tile(colour = "white", linewidth = 1) +
  geom_text(aes(label = sprintf("%+.2f", delta)),
             fontface = "bold", size = 5,
             colour = ifelse(abs(long_full$delta) > max_abs * 0.55,
                               "white", "grey15")) +
  geom_vline(xintercept = sep_x, colour = "white", linewidth = 4) +
  geom_vline(xintercept = sep_x, colour = "grey40", linewidth = 1.0) +
  scale_fill_gradient2(low = "#2166AC", mid = "#F7F7F7", high = "#B2182B",
                          midpoint = 0, limits = c(-max_abs, max_abs),
                          name = expression(Delta[v]^"LOO" ~ "(°C)")) +
  scale_y_discrete(limits = rev,
                     labels = function(x) parse(text = as.character(x))) +
  labs(x = NULL, y = NULL) +
  theme_bw(base_size = 15) +
  theme(axis.text.x = element_text(face = "bold", size = 14),
         axis.text.y = element_text(face = "bold", size = 14),
         panel.grid = element_blank(),
         legend.position = "right")
fig_hm <- file.path(OUT, "fig_heatmap_LOO_profiles_gradient.png")
ggsave(fig_hm, p_hm, width = 12, height = 5, dpi = 300)
cli_alert_success("Saved {.path {fig_hm}}")

# ----------------------------------------------------------------------------
# 2) HOBO forward-selection figure
# ----------------------------------------------------------------------------
cli_h1("HOBO forward selection (Hmax fix)")

forward_order <- ah$forward_order
cli_alert("Forward order : {paste(forward_order, collapse=' → ')}")

df_macro <- extract_macro_daily(CFG$forcing_file, CFG$date_seq)
cg <- readRDS(here::here("outputs/lovb/data/DT_cross_gam_targets.rds"))
hobo_obs <- as.data.table(cg$HOBO)[, .(id_plot = id, dT_obs = -dTmax_mean)]

# Bit codes for forward steps
var_pos <- c(LAI = 1, Hmax = 2, fCover = 3, LAD = 4)
bits_seq <- character(5)
bits_seq[1] <- "0000"
bs <- rep("0", 4)
for (i in seq_along(forward_order)) {
  bs[var_pos[forward_order[i]]] <- "1"
  bits_seq[i + 1L] <- paste(bs, collapse = "")
}

# Labels for the figure
build_label <- function(k) {
  if (k == 1) return("Baseline")
  cum_vars <- forward_order[1:(k-1)]
  cum_vars_disp <- ifelse(cum_vars == "Hmax", "H_max", cum_vars)
  paste("+", paste(cum_vars_disp, collapse = " + "))
}
labels <- vapply(seq_len(5), build_label, character(1))
cli_alert("Step labels : {paste(labels, collapse=' | ')}")

# Collect HOBO sims
hobo_root <- here::here("out_files/musica_hobo_hmaxfix")
collect_step <- function(bit) {
  if (bit == "1111") {
    # REF coalition lives in REF/ from the earlier run
    sub <- here::here(file.path(hobo_root, "REF"))
  } else {
    sub <- file.path(hobo_root, bit)
  }
  nc_files <- list.files(sub, pattern = "\\.nc$", full.names = TRUE)
  rows <- list()
  for (f in nc_files) {
    id <- sub("musica_out_HOBO_(.+)\\.nc$", "\\1", basename(f))
    res <- tryCatch(
      extract_deltatmax_one(f, df_macro, CFG$date_seq,
                              z_target = CFG$tair_target_height),
      error = function(e) NULL)
    if (is.null(res) || nrow(res) == 0) next
    rows[[id]] <- data.table(id_plot = id,
                                 Tmax_micro = mean(res$Tmax_micro, na.rm = TRUE),
                                 dT_sim = mean(res$Delta_Tmax, na.rm = TRUE))
  }
  rbindlist(rows)
}

long_fwd <- list()
for (k in seq_along(bits_seq)) {
  bit <- bits_seq[k]
  df_k <- collect_step(bit)
  if (nrow(df_k) == 0) {
    cli_alert_warning("Step {k} (bit={bit}) : 0 sensors collected")
    next
  }
  df_k <- merge(df_k, hobo_obs, by = "id_plot")
  df_k[, step := k]
  df_k[, label := labels[k]]
  long_fwd[[k]] <- df_k
}
long_fwd <- rbindlist(long_fwd, fill = TRUE)
long_fwd[, label := factor(label, levels = labels)]
cli_alert("HOBO collected : {nrow(long_fwd)} rows  ({uniqueN(long_fwd$id_plot)} sensors)")

# Compute per-step r, RMSE, MAE
fits <- long_fwd[, {
  ok <- !is.na(dT_sim) & !is.na(dT_obs)
  x <- dT_obs[ok]; y <- dT_sim[ok]
  r <- suppressWarnings(cor(x, y))
  rmse <- sqrt(mean((y - x)^2))
  mae <- mean(abs(y - x))
  .(n = length(x), r = r, RMSE = rmse, MAE = mae)
}, by = .(step, label)]
cli_alert("Per-step fits :"); print(fits)
fwrite(fits, file.path(OUT, "tab_forward_selection_fits_hmaxfix.csv"))

# Plot
xy_lim <- range(c(long_fwd$dT_sim, long_fwd$dT_obs), na.rm = TRUE)
xy_pad <- diff(xy_lim) * 0.06
xy_lim <- c(xy_lim[1] - xy_pad, xy_lim[2] + xy_pad)

ann <- fits[, .(label,
                  lab_r    = sprintf("italic(r) == %+.2f", r),
                  lab_rmse = sprintf("RMSE == %.2f * ' °C'", RMSE),
                  lab_mae  = sprintf("MAE == %.2f * ' °C'", MAE))]

p_fwd <- ggplot(long_fwd, aes(x = dT_obs, y = dT_sim)) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", colour = "grey55") +
  geom_hline(yintercept = 0, linetype = "dotted", colour = "grey80") +
  geom_vline(xintercept = 0, linetype = "dotted", colour = "grey80") +
  geom_point(colour = "#2C5F2D", alpha = 0.7, size = 1.8) +
  geom_smooth(method = "lm", formula = y ~ x, se = TRUE,
                colour = "#FFB400", fill = "#FFB400",
                linewidth = 0.9, alpha = 0.22) +
  geom_text(data = ann, aes(x = -Inf, y = Inf, label = lab_r),
             parse = TRUE, hjust = -0.08, vjust = 1.4, size = 6,
             inherit.aes = FALSE, colour = "grey15") +
  geom_text(data = ann, aes(x = -Inf, y = Inf, label = lab_rmse),
             parse = TRUE, hjust = -0.08, vjust = 2.8, size = 6,
             inherit.aes = FALSE, colour = "grey15") +
  geom_text(data = ann, aes(x = -Inf, y = Inf, label = lab_mae),
             parse = TRUE, hjust = -0.08, vjust = 4.2, size = 6,
             inherit.aes = FALSE, colour = "grey15") +
  facet_wrap(~ label, nrow = 1) +
  coord_cartesian(xlim = xy_lim, ylim = xy_lim) +
  labs(x = bquote("Field " ~ Delta * T[max] ~ "(°C)"),
        y = bquote("Simulated " ~ Delta * T[max] ~ "(°C)")) +
  theme_bw(base_size = 14) +
  theme(strip.background = element_rect(fill = "grey92", colour = NA),
         strip.text = element_text(face = "bold", size = 13),
         panel.grid.minor = element_blank())
fig_fwd <- file.path(OUT, "fig_forward_selection_HOBO.png")
ggsave(fig_fwd, p_fwd, width = 18, height = 5, dpi = 300)
cli_alert_success("Saved {.path {fig_fwd}}")

cli_h1("Figures done.")
