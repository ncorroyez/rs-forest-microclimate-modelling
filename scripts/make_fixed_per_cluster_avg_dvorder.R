# ==============================================================================
# Per-plot mean AVG fig — same 5×5 grid as fig_forward_HOBO_v10_fixed_per_cluster_plus_all
# BUT each row uses its OWN forward order, derived from the heatmap LOO Δv values:
#
#   Δv = Tmax(REF) − Tmax(LOO_v) per cluster (microclimat convention).
#   Order = Δv DESCENDING : highest Δv (least buffer / most anti-buffer) first,
#   lowest Δv (strongest buffer) last → "weakest buffer → strongest buffer".
#
# Output : fig_forward_HOBO_v10_fixed_per_cluster_plus_all_AVG.{png,pdf}
#          (overwrites the previous fixed-order AVG version).
# ==============================================================================

suppressMessages({
  library(data.table); library(here); library(cli); library(ncdf4)
  library(tidyverse); library(lubridate); library(musica.tools); library(sf)
})

source(here::here("R/config.R"))
source(here::here("R/io.R"))
source(here::here("R/musica.R"))
source(here::here("R/cluster_relabel.R"))
source(here::here("R/h1_shapley_archetypes.R"))

OUT <- here::here("outputs/figs_MEB2026_final")
df_macro <- extract_macro_daily(CFG$forcing_file, CFG$date_seq)

# ---- HOBO cluster assignment -----------------------------------------------
df_forest <- as.data.table(readRDS(here::here(
  "out_files/Sensitivity_Analysis/clhs_sample_floor05_v2.rds")))[
  , .(x, y, Cluster = as.character(Cluster))]
hobo_pts <- sf::st_read(CFG$hobo_geojson, quiet = TRUE) %>%
  filter(!id_plot %in% CFG$ids_to_remove)
hobo_xy <- sf::st_coordinates(hobo_pts)
df_hobo_xy <- data.table(id_plot = hobo_pts$id_plot,
                          x = hobo_xy[, "X"], y = hobo_xy[, "Y"])
df_hobo_clu <- df_hobo_xy[, .(id_plot, Cluster = sapply(seq_len(.N), function(i) {
  d2 <- (df_forest$x - x[i])^2 + (df_forest$y - y[i])^2
  df_forest$Cluster[which.min(d2)]
}))]
df_hobo_clu[, Cluster := relabel_cluster(Cluster)]

hobo_obs <- as.data.table(read_hobo_daily(CFG$hobo_temp_csv, CFG$date_seq,
                                           df_macro, CFG$ids_to_remove))[
  , .(id_plot, date, Delta_obs)]

# ---- HOBO sim per coalition (per-plot mean) --------------------------------
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
    rows[[id]] <- as.data.table(res)[, .(id_plot = id, date,
                                          Delta_sim = Delta_Tmax)]
  }
  rbindlist(rows)
}

per_plot_mean_bit <- function(bit) {
  s <- collect_step_daily(bit)
  if (nrow(s) == 0) return(NULL)
  m <- merge(s, hobo_obs, by = c("id_plot", "date"))
  m <- merge(m, df_hobo_clu, by = "id_plot")
  m[, .(Delta_sim = mean(Delta_sim, na.rm = TRUE),
        Delta_obs = mean(Delta_obs, na.rm = TRUE)),
    by = .(id_plot, Cluster)]
}

# Pre-compute all 16 coalitions at per-plot mean level
cli_h1("Pre-computing all 16 HOBO coalitions (per-plot mean)")
PP <- list()
for (bit in unname(.ARCH_COAL_MAP)) PP[[bit]] <- per_plot_mean_bit(bit)

# ---- Heatmap Δv per cluster (from V10 archetype factorial) -----------------
arch_root  <- here::here("out_files/H1_archetypes_v9")
arch_v10_r <- here::here("out_files/H1_archetypes_v10_fcovmean")
sc_name_of <- setNames(names(.ARCH_COAL_MAP), unname(.ARCH_COAL_MAP))

arch_nc_path <- function(cl, bit) {
  arch_lbl <- sprintf("Arch_C%s", cl)
  v10 <- file.path(arch_v10_r, arch_lbl,
                    sprintf("musica_out_%s_%s.nc", arch_lbl, sc_name_of[bit]))
  v9  <- file.path(arch_root, arch_lbl,
                    sprintf("musica_out_%s_%s.nc", arch_lbl, sc_name_of[bit]))
  if (file.exists(v10)) v10 else v9
}
arch_tmax <- function(cl, bit) {
  res <- tryCatch(extract_deltatmax_one(arch_nc_path(cl, bit), df_macro,
                                         CFG$date_seq,
                                         z_target = CFG$tair_target_height),
                  error = function(e) NULL)
  if (is.null(res) || nrow(res) == 0) return(NA_real_)
  mean(res$Tmax_micro, na.rm = TRUE)
}

VARS     <- c("LAI", "Hmax", "fCover", "LAD")
loo_bits <- c(LAI = "0111", Hmax = "1011", fCover = "1101", LAD = "1110")

cli_h1("Computing Δv per cluster (heatmap convention REF − LOO_v)")
arch_clusters <- c("3", "2", "4", "1")   # old IDs in P1..P4 order (ASC LAI)
new_labels    <- c("P1", "P2", "P3", "P4")
dv_per_cluster <- list()
for (i in seq_along(arch_clusters)) {
  cl_old <- arch_clusters[i]; profile <- new_labels[i]
  tmax_ref <- arch_tmax(cl_old, "1111")
  dv <- sapply(VARS, function(v) tmax_ref - arch_tmax(cl_old, loo_bits[v]))
  dv_per_cluster[[profile]] <- dv
}
dv_avg <- colMeans(do.call(rbind, dv_per_cluster), na.rm = TRUE)
dv_per_cluster[["All"]] <- dv_avg
cli_alert("Δv values per cluster (− = buffers, + = anti-buffers) :")
print(do.call(rbind, dv_per_cluster))

# ---- Per-cluster order : Δv DESCENDING (least buffer first → most buffer last)
cli_h1("Per-cluster forward order (Δv descending)")
order_per_cluster <- lapply(dv_per_cluster, function(dv) names(sort(dv, decreasing = TRUE)))
for (cl in names(order_per_cluster)) {
  cli_alert("{cl} : {paste(order_per_cluster[[cl]], collapse=' → ')}")
}

# Build bit sequence per cluster ; labels are PLOTMATH strings so H_max renders
# as italic H with subscript max.
clusters <- c("P1", "P2", "P3", "P4", "All")
bit_from_set <- function(set) paste(ifelse(VARS %in% set, "1", "0"), collapse = "")
to_pm <- function(v) ifelse(v == "Hmax", "italic(H)[max]", v)
paths <- list()
for (cl in clusters) {
  ord <- order_per_cluster[[cl]]
  bits <- c(bit_from_set(character(0)))
  active <- character(0)
  for (v in ord) { active <- c(active, v); bits <- c(bits, bit_from_set(active)) }
  labels_step <- c("'Baseline'",
                    sapply(seq_along(ord), function(k) {
                      cum_vars <- sapply(ord[1:k], to_pm)
                      joined <- paste(cum_vars, collapse = "*' + '*")
                      if (k == length(ord)) paste0("'REF: '*", joined)
                      else paste0("'+ '*", joined)
                    }))
  paths[[cl]] <- data.table(step = seq_along(bits), Cluster = cl,
                             bit = bits, label = labels_step)
}

# ---- Compute per-(cluster, step) stats and long data ------------------------
cli_h1("Build long data + fits")
long <- list(); ann <- list()
for (cl in clusters) {
  p <- paths[[cl]]
  for (k in seq_len(nrow(p))) {
    bit <- p$bit[k]
    d <- PP[[bit]]
    if (is.null(d)) next
    panel_d <- if (cl == "All") copy(d) else d[Cluster == cl]
    panel_d[, Cluster_row := cl]; panel_d[, step := k]
    panel_d[, label := p$label[k]]
    long[[length(long) + 1L]] <- panel_d
    ok <- !is.na(panel_d$Delta_sim) & !is.na(panel_d$Delta_obs)
    x <- panel_d$Delta_obs[ok]; y <- panel_d$Delta_sim[ok]
    n <- length(x)
    r_val <- if (n < 3) NA_real_ else suppressWarnings(cor(x, y))
    rmse_val <- if (n == 0) NA_real_ else sqrt(mean((y - x)^2))
    ann[[length(ann) + 1L]] <- data.table(
      Cluster_row = cl, step = k, label = p$label[k],
      r_val = r_val, rmse_val = rmse_val, n_val = n)
  }
}
long <- rbindlist(long, fill = TRUE)
ann  <- rbindlist(ann)
long[, Cluster_row := factor(Cluster_row, levels = clusters)]
ann[ , Cluster_row := factor(Cluster_row, levels = clusters)]

ann[, lab_r    := ifelse(is.na(r_val), NA_character_,
                          sprintf(" italic(r)==%+.2f", r_val))]
ann[, lab_rmse := sprintf(" RMSE==%.2f * ' °C'", rmse_val)]
ann[, lab_n    := sprintf(" n==%d", n_val)]

# Per-cell variable label (bottom-right, shows what's added at this step)
ann[, panel_lab := label]

xy_lim <- range(c(long$Delta_sim, long$Delta_obs), na.rm = TRUE)
xy_pad <- diff(xy_lim) * 0.05; xy_lim <- xy_lim + c(-xy_pad, xy_pad)

STAT_SIZE  <- 6
STAT_FACE  <- "bold"
STAT_COLOUR <- "grey15"

p <- ggplot(long, aes(x = Delta_obs, y = Delta_sim,
                       colour = Cluster, shape = Cluster)) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", colour = "grey55") +
  geom_hline(yintercept = 0, linetype = "dotted", colour = "grey80") +
  geom_vline(xintercept = 0, linetype = "dotted", colour = "grey80") +
  geom_point(alpha = 0.85, size = 3.2) +
  geom_smooth(method = "lm", formula = y ~ x, se = TRUE,
              colour = "#FFB400", fill = "#FFB400",
              linewidth = 0.9, alpha = 0.22, aes(group = 1)) +
  geom_text(data = ann[!is.na(lab_r)],
             aes(x = -Inf, y = Inf, label = lab_r),
             parse = TRUE, hjust = -0.10, vjust = 1.4,
             size = STAT_SIZE, fontface = STAT_FACE, colour = STAT_COLOUR,
             inherit.aes = FALSE) +
  geom_text(data = ann, aes(x = -Inf, y = Inf, label = lab_rmse),
             parse = TRUE, hjust = -0.10, vjust = 2.8,
             size = STAT_SIZE, fontface = STAT_FACE, colour = STAT_COLOUR,
             inherit.aes = FALSE) +
  geom_text(data = ann, aes(x = -Inf, y = Inf, label = lab_n),
             parse = TRUE, hjust = -0.10, vjust = 4.2,
             size = STAT_SIZE, fontface = STAT_FACE, colour = STAT_COLOUR,
             inherit.aes = FALSE) +
  geom_text(data = ann, aes(x = Inf, y = -Inf, label = panel_lab),
             hjust = 1.05, vjust = -0.5,
             size = 5, fontface = "italic", parse = TRUE,
             inherit.aes = FALSE, colour = "grey25") +
  facet_grid(Cluster_row ~ step,
              labeller = labeller(step = c("1" = "Baseline",
                                              "2" = "Step 1",
                                              "3" = "Step 2",
                                              "4" = "Step 3",
                                              "5" = "Step 4"))) +
  coord_cartesian(xlim = xy_lim, ylim = xy_lim) +
  scale_colour_manual(values = c(PAL_CLUSTER), na.translate = FALSE) +
  scale_shape_manual(values = c(PAL_SHAPES), na.translate = FALSE) +
  labs(x = expression("Field " ~ Delta * T[max] ~ "(°C)"),
       y = expression("Simulated " ~ Delta * T[max] ~ "(°C)")) +
  theme_bw(base_size = 18) +
  theme(strip.background = element_rect(fill = "grey92", colour = NA),
        strip.text  = element_text(face = "bold", size = 26),
        axis.title  = element_text(face = "bold", size = 28),
        axis.text   = element_text(size = 18),
        legend.text  = element_text(size = 24),
        legend.title = element_text(size = 26),
        legend.position = "bottom",
        legend.key.width = unit(2, "cm"),
        panel.grid.minor = element_blank(),
        aspect.ratio = 0.6)

# Slide 16:9 — 24×13.5 in @ 600 dpi PNG + PDF vector
ggsave(file.path(OUT, "fig_forward_HOBO_v10_fixed_per_cluster_plus_all_AVG.png"),
        p, width = 24, height = 13.5, dpi = 600, bg = "white")
ggsave(file.path(OUT, "fig_forward_HOBO_v10_fixed_per_cluster_plus_all_AVG.pdf"),
        p, width = 24, height = 13.5, device = cairo_pdf)
cli_alert_success("Saved fig_forward_HOBO_v10_fixed_per_cluster_plus_all_AVG (png + pdf, 16:9, per-cluster Δv order)")

# Stats CSV
fwrite(ann[, .(Cluster_row, step, label,
                r = round(r_val, 3), RMSE = round(rmse_val, 3), n = n_val)],
        file.path(OUT, "tab_v10_HOBO_fixed_per_cluster_plus_all_AVG_dvorder.csv"))
cli_alert_success("Saved tab CSV")
