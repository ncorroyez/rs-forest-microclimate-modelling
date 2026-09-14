# ==============================================================================
# V2 (Hmax MAX aggregation) — full pipeline :
#   - 4 NEW archetypes × 16 coalitions = 64 sims  (out_files/H1_archetypes_v2/)
#   - 53 HOBO sensors × 5 forward-selection coalitions  (out_files/musica_hobo_v2/)
#
# HOBO Hmax now extracted via cell-MAX 20 m (no buffer needed — load_lidar_rasters
# already does MAX agg for Hmax). LAI / fCover / LAD remain cell-MEAN, consistent.
# ==============================================================================

suppressMessages({
  library(data.table); library(here); library(cli); library(ncdf4)
  library(sf); library(terra); library(future); library(furrr)
  library(tidyverse); library(musica.tools); library(rmusica); library(lubridate)
})

source(here::here("R/config.R"))
source(here::here("R/io.R"))
source(here::here("R/musica.R"))
source(here::here("R/lad.R"))
source(here::here("R/validation.R"))
source(here::here("R/h1_shapley_archetypes.R"))

N_WORKERS <- 4
plan(multisession, workers = N_WORKERS)

# ----------------------------------------------------------------------------
# Phase 1 : Archetype sims (4 × 16 = 64)
# ----------------------------------------------------------------------------
cli_h1("Phase 1 — Archetype sims (V2)")

df_floor <- as.data.table(readRDS(here::here(
  "out_files/Sensitivity_Analysis/clhs_sample_floor05_v2.rds")))
df_archetypes <- make_synthetic_archetypes(as.data.frame(df_floor))
cli_alert("V2 archetypes (with Hmax MAX) :")
print(df_archetypes[, c("Cluster", "LAI", "Hmax", "fCover")])

fac_scs <- build_factorial_scenarios_archetypes(df_floor)
out_root_arch <- here::here("out_files/H1_archetypes_v2")
dir.create(out_root_arch, recursive = TRUE, showWarnings = FALSE)

jobs <- list()
for (i in seq_len(nrow(df_archetypes))) {
  arch <- df_archetypes[i, , drop = FALSE]
  arch_lbl <- sprintf("Arch_C%s", arch$Cluster)
  out_dir <- file.path(out_root_arch, arch_lbl)
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  for (sc_nm in names(fac_scs)) {
    sc <- fac_scs[[sc_nm]]
    out_nc <- file.path(out_dir, sprintf("musica_out_%s_%s.nc", arch_lbl, sc_nm))
    if (file.exists(out_nc) && file.size(out_nc) > 1e6) next
    jobs[[length(jobs) + 1L]] <- list(arch = arch, sc = sc, out_nc = out_nc,
                                          lbl = paste(arch_lbl, sc_nm))
  }
}
cli_alert("Jobs to run : {length(jobs)}")

t_arch <- Sys.time()
res_arch <- future_walk(jobs, function(job) {
  suppressMessages({
    library(here); library(ncdf4); library(sf); library(terra)
    library(tidyverse); library(musica.tools); library(rmusica); library(lubridate)
  })
  source(here::here("R/config.R")); source(here::here("R/io.R"))
  source(here::here("R/musica.R")); source(here::here("R/lad.R"))
  tryCatch(
    run_musica_one(job$arch, job$sc, job$out_nc,
                    CFG$forcing_file, CFG$musica_cmd),
    error = function(e) cat(sprintf("[%s] ERROR : %s\n", job$lbl, e$message)))
}, .options = furrr_options(seed = TRUE))
cli_alert_success("Archetype sims done in {round(as.numeric(difftime(Sys.time(), t_arch, units='mins')), 1)} min")

# ----------------------------------------------------------------------------
# Phase 2 : Determine forward-selection order from archetype LOO
# ----------------------------------------------------------------------------
cli_h1("Phase 2 — Determine forward order")

df_macro <- extract_macro_daily(CFG$forcing_file, CFG$date_seq)
LADMAP <- .ARCH_COAL_MAP

extract_one <- function(nc_path) {
  if (!file.exists(nc_path) || file.size(nc_path) < 1e6) return(NA_real_)
  res <- tryCatch(extract_deltatmax_one(nc_path, df_macro, CFG$date_seq,
                                            z_target = CFG$tair_target_height),
                    error = function(e) NULL)
  if (is.null(res) || nrow(res) == 0) return(NA_real_)
  mean(res$Tmax_micro, na.rm = TRUE)
}

arch_rows <- list()
for (i in seq_len(nrow(df_archetypes))) {
  arch_lbl <- sprintf("Arch_C%s", df_archetypes$Cluster[i])
  for (sc_nm in names(fac_scs)) {
    bit <- LADMAP[sc_nm]
    nc <- file.path(out_root_arch, arch_lbl,
                      sprintf("musica_out_%s_%s.nc", arch_lbl, sc_nm))
    tair <- extract_one(nc)
    arch_rows[[length(arch_rows) + 1L]] <- data.table(
      Cluster = df_archetypes$Cluster[i],
      bit_code = bit, Tmax_micro = tair)
  }
}
DT_arch <- rbindlist(arch_rows)
DT_arch_w <- dcast(DT_arch, Cluster ~ bit_code, value.var = "Tmax_micro")

VARS <- c("LAI","Hmax","fCover","LAD")
loo_bits <- c(LAI = "0111", Hmax = "1011", fCover = "1101", LAD = "1110")
delta <- list()
for (v in VARS) delta[[v]] <- DT_arch_w[["1111"]] - DT_arch_w[[loo_bits[v]]]
DT_delta <- data.table(Cluster = DT_arch_w$Cluster, do.call(cbind, delta))
DT_delta_long <- melt(DT_delta, id.vars = "Cluster",
                          variable.name = "variable", value.name = "delta")
DT_avg <- DT_delta_long[, .(Avg = mean(delta)), by = variable]
cli_alert("Avg Δ_v per variable :"); print(DT_avg)

setorder(DT_avg, -Avg)
forward_order <- as.character(DT_avg$variable)
cli_alert("Forward order (anti-buffer → strong buffer) : {paste(forward_order, collapse=' → ')}")

saveRDS(list(arch = DT_arch, delta = DT_delta_long, avg = DT_avg,
              forward_order = forward_order),
         here::here("outputs/v2_archetype_heatmap_data.rds"))

# ----------------------------------------------------------------------------
# Phase 3 : HOBO sims (53 × 5 coalitions in forward order)
# ----------------------------------------------------------------------------
cli_h1("Phase 3 — HOBO forward-selection sims")

rasters <- load_lidar_rasters(CFG$in_dir, CFG$agg_factor)
# Default extraction now : cell-MAX Hmax (consistent with cLHS V2)
df_hobo <- as.data.frame(build_hobo_inputs(
  CFG$hobo_geojson, rasters$stack, CFG$ids_to_remove))
df_hobo[df_hobo$fCover < 0.5, "fCover"] <- 0.5
cli_alert("n HOBO sensors : {nrow(df_hobo)}")
cli_alert("HOBO Hmax stats : mean={round(mean(df_hobo$Hmax),2)} m, "
            ~ "range [{round(min(df_hobo$Hmax),2)}, {round(max(df_hobo$Hmax),2)}]")

var_pos <- c(LAI = 1, Hmax = 2, fCover = 3, LAD = 4)
bits_seq <- character(5)
bits_seq[1] <- "0000"
bs <- rep("0", 4)
for (i in seq_along(forward_order)) {
  bs[var_pos[forward_order[i]]] <- "1"
  bits_seq[i + 1L] <- paste(bs, collapse = "")
}
cli_alert("Forward bit sequence : {paste(bits_seq, collapse=' → ')}")

bit2scn <- setNames(names(.ARCH_COAL_MAP), unname(.ARCH_COAL_MAP))
fac_scs_hobo <- build_factorial_scenarios_archetypes(df_floor)

out_root_hobo <- here::here("out_files/musica_hobo_v2")
hobo_jobs <- list()
for (i in seq_len(nrow(df_hobo))) {
  pr <- df_hobo[i, , drop = FALSE]
  for (k in seq_along(bits_seq)) {
    bit <- bits_seq[k]
    sc_name <- bit2scn[bit]
    if (is.na(sc_name)) next
    sc <- fac_scs_hobo[[sc_name]]
    out_dir <- file.path(out_root_hobo, bit)
    dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
    out_nc <- file.path(out_dir, sprintf("musica_out_HOBO_%s.nc", pr$id_plot))
    if (file.exists(out_nc) && file.size(out_nc) > 1e6) next
    hobo_jobs[[length(hobo_jobs) + 1L]] <- list(
      pr = pr, sc = sc, out_nc = out_nc,
      lbl = sprintf("%s/%s", bit, pr$id_plot))
  }
}
cli_alert("HOBO jobs to run : {length(hobo_jobs)}")

t_hobo <- Sys.time()
res_hobo <- future_walk(hobo_jobs, function(job) {
  suppressMessages({
    library(here); library(ncdf4); library(sf); library(terra)
    library(tidyverse); library(musica.tools); library(rmusica); library(lubridate)
  })
  source(here::here("R/config.R")); source(here::here("R/io.R"))
  source(here::here("R/musica.R")); source(here::here("R/lad.R"))
  tryCatch(
    run_musica_one(job$pr, job$sc, job$out_nc,
                    CFG$forcing_file, CFG$musica_cmd),
    error = function(e) cat(sprintf("[%s] ERROR : %s\n", job$lbl, e$message)))
}, .options = furrr_options(seed = TRUE))
cli_alert_success("HOBO sims done in {round(as.numeric(difftime(Sys.time(), t_hobo, units='mins')), 1)} min")

# ----------------------------------------------------------------------------
# Phase 4 : Regenerate figures (heatmap + forward selection)
# ----------------------------------------------------------------------------
cli_h1("Phase 4 — Regenerate figures")

OUT <- here::here("outputs/figs_MEB2026_final")

# Heatmap LOO × profiles (with Avg col)
long <- copy(DT_delta_long)
# Cluster = RAW k-means code; map via cluster_relabel (raw 1->P4, 2->P2, 3->P1, 4->P3).
# Direct P%d was WRONG (scrambled archetype labels).
.RELAB <- c("1" = "P4", "2" = "P2", "3" = "P1", "4" = "P3")
long[, Profile := unname(.RELAB[as.character(Cluster)])]
avg_rows <- long[, .(Profile = "Avg", delta = mean(delta, na.rm = TRUE)),
                    by = variable]
long_full <- rbind(long[, .(variable, Profile, delta)], avg_rows)
long_full[, Profile := factor(Profile, levels = c(paste0("P", 1:4), "Avg"))]
long_full[, variable := factor(variable, levels = VARS)]
long_full[, var_lab := ifelse(variable == "Hmax", "H[max]", as.character(variable))]
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

# Forward selection
cg <- readRDS(here::here("outputs/lovb/data/DT_cross_gam_targets.rds"))
hobo_obs <- as.data.table(cg$HOBO)[, .(id_plot = id, dT_obs = -dTmax_mean)]

build_label <- function(k) {
  if (k == 1) return("Baseline")
  cum_vars <- forward_order[1:(k-1)]
  cum_vars_disp <- ifelse(cum_vars == "Hmax", "H_max", cum_vars)
  paste("+", paste(cum_vars_disp, collapse = " + "))
}
labels <- vapply(seq_len(5), build_label, character(1))

collect_step <- function(bit) {
  sub <- file.path(out_root_hobo, bit)
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
                                 dT_sim = mean(res$Delta_Tmax, na.rm = TRUE))
  }
  rbindlist(rows)
}

long_fwd <- list()
for (k in seq_along(bits_seq)) {
  df_k <- collect_step(bits_seq[k])
  if (nrow(df_k) == 0) next
  df_k <- merge(df_k, hobo_obs, by = "id_plot")
  df_k[, step := k]
  df_k[, label := labels[k]]
  long_fwd[[k]] <- df_k
}
long_fwd <- rbindlist(long_fwd, fill = TRUE)
long_fwd[, label := factor(label, levels = labels)]

fits <- long_fwd[, {
  ok <- !is.na(dT_sim) & !is.na(dT_obs)
  x <- dT_obs[ok]; y <- dT_sim[ok]
  r <- suppressWarnings(cor(x, y))
  rmse <- sqrt(mean((y - x)^2))
  mae <- mean(abs(y - x))
  .(n = length(x), r = r, RMSE = rmse, MAE = mae)
}, by = .(step, label)]
cli_alert("Per-step fits :"); print(fits)
fwrite(fits, file.path(OUT, "tab_forward_selection_fits_v2.csv"))

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

cli_h1("V2 pipeline complete.")
