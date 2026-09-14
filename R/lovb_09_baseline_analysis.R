# ==============================================================================
# LOVB analysis — Module 9 : baseline robustness analysis (Option B)
#
# After Module 8 has run :
#   - 13 plots/cluster x 4 clusters = 52 subset plots
#   - 5 coalitions x 2 new baselines (Q10, Q90) = 10 scenarios
#   - 520 new NCs in out_files/H1_lovb_baselines/BASE_<label>_<bit>/
#
# Plus REF (1111) from existing H1F_r_r_r_r (shared across baselines).
#
# Compute LOVB rankings under each baseline (mean + Q10 + Q90), compute Kendall
# tau between them, generate figure.
# ==============================================================================

suppressMessages({
  library(data.table)
  library(here)
  library(cli)
  library(ggplot2)
  library(ncdf4)
  library(musica.tools); library(rmusica)
  library(tidyverse)
})

source(here::here("R/config.R"))
source(here::here("R/io.R"))
source(here::here("R/lad.R"))
source(here::here("R/musica.R"))

.VARS    <- c("LAI", "Hmax", "fCover", "LAD")
.LOVB_BITS_TO_VAR <- c("0111" = "LAI", "1011" = "Hmax", "1101" = "fCover", "1110" = "LAD")

# ---- 1. Extract daily Delta_Tmax for the subset plots under each baseline ----
lovb_load_daily_baselines <- function(cache_path = here::here("outputs/lovb/data/DT_daily_baselines.rds"),
                                       use_cache = TRUE) {
  if (use_cache && file.exists(cache_path)) {
    return(as.data.table(readRDS(cache_path)))
  }
  cli_h1("Extract daily Delta_Tmax for baseline subset")
  df_macro <- extract_macro_daily(CFG$forcing_file, CFG$date_seq)

  out <- list()
  for (lbl in c("Q10", "Q90")) {
    for (bit in c("0000", "0111", "1011", "1101", "1110")) {
      sc_dir <- here::here("out_files/H1_lovb_baselines", sprintf("BASE_%s_%s", lbl, bit))
      if (!dir.exists(sc_dir)) { cli_alert_warning("Missing : {.path {sc_dir}}"); next }
      nc_files <- list.files(sc_dir, pattern = "\\.nc$", full.names = TRUE)
      nc_files <- nc_files[file.size(nc_files) > 1e6]
      cli_alert("[{lbl} / {bit}] {length(nc_files)} NCs")
      for (f in nc_files) {
        bn <- basename(f)
        x  <- as.integer(sub(".*_X(\\d+)_Y\\d+\\.nc$", "\\1", bn))
        y  <- as.integer(sub(".*_X\\d+_Y(\\d+)\\.nc$", "\\1", bn))
        res <- tryCatch(extract_deltatmax_one(f, df_macro, CFG$date_seq),
                         error = function(e) NULL)
        if (is.null(res) || nrow(res) == 0L) next
        dt <- as.data.table(res)
        dt[, x := x]; dt[, y := y]; dt[, bit_code := bit]; dt[, baseline := lbl]
        out[[length(out) + 1L]] <- dt[, .(x, y, date, bit_code, baseline, Delta_Tmax)]
      }
    }
  }

  # REF (1111) : reuse existing cLHS factorial NCs for the same subset coords
  df_clhs <- as.data.table(readRDS(here::here("outputs/lovb/data/DT_daily_cLHS.rds")))
  subset_coords <- unique(rbindlist(out, fill = TRUE)[, .(x, y)])
  ref_clhs <- df_clhs[bit_code == "1111"][subset_coords, on = c("x", "y"), nomatch = NULL]
  ref_clhs[, baseline := "Q10"]
  ref_clhs2 <- copy(ref_clhs); ref_clhs2[, baseline := "Q90"]
  out[[length(out) + 1L]] <- ref_clhs
  out[[length(out) + 1L]] <- ref_clhs2

  DT <- rbindlist(out, fill = TRUE)
  setkey(DT, baseline, bit_code, x, y, date)

  dir.create(dirname(cache_path), recursive = TRUE, showWarnings = FALSE)
  saveRDS(DT, cache_path)
  cli_alert_success("Cached {.path {cache_path}} ({.val {nrow(DT)}} rows)")
  DT
}

# ---- 2. Compute LOVB metrics per baseline ------------------------------------
lovb_metrics_baselines <- function(DT_daily) {
  cli_h2("LOVB metrics per baseline")
  # Aggregate
  DT_agg <- DT_daily[
    , .(Tmax_mean = mean(Delta_Tmax, na.rm = TRUE),
        Tmax_P90  = quantile(Delta_Tmax, 0.9, type = 7, na.rm = TRUE)),
    by = .(baseline, x, y, bit_code)
  ]

  # Wide pivot per baseline
  out <- list()
  for (lbl in c("Q10", "Q90")) {
    DT_l <- DT_agg[baseline == lbl]
    DT_w <- dcast(DT_l, x + y ~ bit_code,
                   value.var = c("Tmax_mean", "Tmax_P90"))
    if (!all(c("Tmax_mean_1111", "Tmax_mean_0111") %in% names(DT_w))) {
      cli_alert_warning("Baseline {lbl} : missing coalitions")
      next
    }
    for (m in c("mean", "P90")) {
      ref_col <- paste0("Tmax_", m, "_1111")
      for (bit in names(.LOVB_BITS_TO_VAR)) {
        v <- .LOVB_BITS_TO_VAR[bit]
        lovb_col <- paste0("Tmax_", m, "_", bit)
        d <- DT_w[[ref_col]] - DT_w[[lovb_col]]
        out[[length(out) + 1L]] <- data.table(
          Baseline = lbl, Variable = v, Metric = m,
          RMSE = sqrt(mean(d^2, na.rm = TRUE)),
          Bias = mean(d, na.rm = TRUE),
          N    = sum(!is.na(d))
        )
      }
    }
  }
  rbindlist(out)
}

# ---- 3. Add mean-baseline RMSE from existing LOVB (already computed) ----------
lovb_add_mean_baseline <- function(DT_metrics_baselines, n_per_cluster_subset = 13L) {
  cli_h2("Add reference 'mean' baseline (from existing DT_contrib_cLHS)")
  DT_c <- as.data.table(readRDS(here::here("outputs/lovb/data/DT_contrib_cLHS.rds")))
  # For fair comparison, use same subset (13 plots/cluster)
  setDT(DT_c)
  set.seed(123)
  DT_c[, row_id := .I]
  DT_sub <- DT_c[, .SD[sample(.N, min(.N, n_per_cluster_subset))], by = Cluster]

  rows <- list()
  for (m in c("mean", "P90")) {
    for (v in .VARS) {
      d <- DT_sub[[paste0("Tmax_", m, "_REF")]] -
            DT_sub[[paste0("Tmax_", m, "_LOVB_", v)]]
      rows[[length(rows) + 1L]] <- data.table(
        Baseline = "mean", Variable = v, Metric = m,
        RMSE = sqrt(mean(d^2, na.rm = TRUE)),
        Bias = mean(d, na.rm = TRUE),
        N    = sum(!is.na(d))
      )
    }
  }
  rbindlist(list(DT_metrics_baselines, rbindlist(rows)), use.names = TRUE)
}

# ---- 4. Rankings + Kendall tau matrix ----------------------------------------
lovb_baseline_convergence <- function(DT_metrics_full,
                                       out_csv = here::here("outputs/lovb/tables/tab_baseline_robustness.csv")) {
  cli_h1("Baseline robustness — rankings + Kendall tau")
  # Rank by RMSE within each (Baseline, Metric)
  DT_metrics_full[, Rank := frank(-RMSE, ties.method = "average"),
                   by = .(Baseline, Metric)]

  dir.create(dirname(out_csv), recursive = TRUE, showWarnings = FALSE)
  fwrite(DT_metrics_full, out_csv)

  cat("\n### Baseline robustness — ranks per variable per baseline x metric\n\n")
  print(DT_metrics_full)

  # Kendall tau matrix : compare rankings across 3 baselines (mean / Q10 / Q90) per metric
  cat("\n### Kendall tau (between baselines, same metric)\n\n")
  for (m in c("mean", "P90")) {
    sub <- dcast(DT_metrics_full[Metric == m], Variable ~ Baseline, value.var = "Rank")
    base_cols <- intersect(c("mean", "Q10", "Q90"), names(sub))
    K <- matrix(NA_real_, length(base_cols), length(base_cols),
                dimnames = list(base_cols, base_cols))
    for (i in seq_along(base_cols)) for (j in seq_along(base_cols))
      K[i, j] <- cor(sub[[base_cols[i]]], sub[[base_cols[j]]], method = "kendall")
    cat(sprintf("\nMetric = %s :\n", m))
    print(round(K, 3))
  }

  cli_alert_success("Saved {.path {out_csv}}")
  invisible(DT_metrics_full)
}

# ---- 5. Figure : RMSE par variable x baseline -------------------------------
lovb_fig_baselines <- function(DT_metrics_full,
                                out_path = here::here("outputs/lovb/figures/fig_LOVB_baseline_robustness.png")) {
  pal_var <- c(LAI = "#440154", fCover = "#31688e", Hmax = "#35b779", LAD = "#fde725")

  DT_metrics_full[, Baseline := factor(Baseline, levels = c("Q10", "mean", "Q90"))]
  DT_metrics_full[, Variable := factor(Variable, levels = c("LAI", "fCover", "Hmax", "LAD"))]

  p <- ggplot(DT_metrics_full,
               aes(x = Baseline, y = RMSE, fill = Variable, group = Variable)) +
    geom_col(position = position_dodge(width = 0.8), width = 0.7,
              colour = "grey25", linewidth = 0.3) +
    geom_text(aes(label = sprintf("%.2f", RMSE)),
               position = position_dodge(width = 0.8), vjust = -0.4,
               size = 3.0, fontface = "bold") +
    facet_wrap(~ Metric, nrow = 1) +
    scale_fill_manual(values = pal_var) +
    scale_y_continuous(expand = expansion(mult = c(0, 0.2))) +
    labs(
      title    = "LOVB robustness to baseline choice",
      subtitle = "RMSE per LiDAR variable across three baselines (Q10, mean, Q90) on a 52-plot subset",
      x = NULL,
      y = "RMSE of Delta_v  [degree C]"
    ) +
    theme_bw(base_size = 11) +
    theme(strip.text      = element_text(face = "bold"),
          plot.title      = element_text(face = "bold"),
          legend.position = "bottom")

  dir.create(dirname(out_path), recursive = TRUE, showWarnings = FALSE)
  ggsave(out_path, p, width = 11, height = 5.5, dpi = 300)
  cli_alert_success("Saved {.path {out_path}}")
  invisible(p)
}

# ---- Main orchestrator -------------------------------------------------------
lovb_run_baseline_analysis <- function() {
  cli_h1("Module 9 : baseline robustness analysis")
  DT_d <- lovb_load_daily_baselines(use_cache = TRUE)
  DT_m <- lovb_metrics_baselines(DT_d)
  DT_m <- lovb_add_mean_baseline(DT_m)
  res  <- lovb_baseline_convergence(DT_m)
  lovb_fig_baselines(DT_m)
  invisible(res)
}

if (sys.nframe() == 0) {
  lovb_run_baseline_analysis()
}
