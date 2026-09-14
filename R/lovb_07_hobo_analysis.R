# ==============================================================================
# LOVB analysis — Module 7 : HOBO LOVB analysis + comparison with observed
#
# After Module 6 has launched the missing coalitions :
#   1) extract daily Delta_Tmax for the 10 coalitions at 53 HOBO sites
#   2) build DT_contrib_HOBO (wide, 53 rows)
#   3) per-plot LOVB Delta_v
#   4) HOBO observed "buffering sensitivity" : slope(Tmax_obs ~ Tmax_macro) per plot
#      regressed against LAI/Hmax/fCover/FPC1 across the 53 plots  ->  ranking obs
#   5) Kendall tau between LOVB simulated ranking and observed ranking
#   6) figures + CSV
# ==============================================================================

suppressMessages({
  library(data.table)
  library(here)
  library(cli)
  library(ggplot2)
  library(ncdf4)
  library(musica.tools); library(rmusica)
  library(tidyverse); library(sf); library(terra)
  library(mgcv)        # GAM multivarié pour observed buffering vs traits
  library(fda)         # FPCA pour projeter HOBO sur la base cLHS
})

source(here::here("R/config.R"))
source(here::here("R/io.R"))
source(here::here("R/lad.R"))
source(here::here("R/forest.R"))
source(here::here("R/fpca.R"))
source(here::here("R/musica.R"))
source(here::here("R/validation.R"))

# Bit code -> HOBO subdirectory name (folder convention from Module 6 + forward)
.LOVB_HOBO_DIRS <- list(
  "0000" = "H1f_0_Null_baseline",       # existing forward
  "0001" = "HOBO_0001",                  # LVA_LAD             (new from Module 6)
  "0010" = "HOBO_0010",                  # LVA_fCover          (new)
  "0100" = "HOBO_0100",                  # LVA_Hmax            (new)
  "1000" = "H1f_1_LAI_only",             # existing forward    (LVA_LAI)
  "0111" = "HOBO_0111",                  # LOVB_LAI            (new)
  "1011" = "HOBO_1011",                  # LOVB_Hmax           (new)
  "1101" = "HOBO_1101",                  # LOVB_fCover         (new)
  "1110" = "H1f_3_LAI_Hmax_fCover",      # existing forward    (LOVB_LAD)
  "1111" = "REF_all_real"                # existing forward    (REF)
)

# ---- 1. Load daily Delta_Tmax at HOBO for the 10 coalitions ------------------
lovb_load_daily_hobo <- function(cache_path = here::here("outputs/lovb/data/DT_daily_HOBO.rds"),
                                  use_cache  = TRUE) {
  if (use_cache && file.exists(cache_path)) {
    cli_alert_info("Loading cached HOBO daily DT from {.path {cache_path}}")
    return(as.data.table(readRDS(cache_path)))
  }
  cli_h1("Extract daily Delta_Tmax at 53 HOBO sites (10 coalitions)")
  df_macro <- extract_macro_daily(CFG$forcing_file, CFG$date_seq)

  out <- list()
  for (bit in names(.LOVB_HOBO_DIRS)) {
    sc_dir <- here::here("out_files/musica_hobo_validation", .LOVB_HOBO_DIRS[[bit]])
    if (!dir.exists(sc_dir)) { cli_alert_warning("Missing : {.path {sc_dir}}"); next }
    nc_files <- list.files(sc_dir, pattern = "\\.nc$", full.names = TRUE)
    nc_files <- nc_files[file.size(nc_files) > 1e6]
    cli_alert("[{bit}] {basename(sc_dir)} : {length(nc_files)} NCs")
    for (f in nc_files) {
      id_str <- sub(".*HOBO_(.*)\\.nc$", "\\1", basename(f))
      res    <- tryCatch(extract_deltatmax_one(f, df_macro, CFG$date_seq),
                          error = function(e) NULL)
      if (is.null(res) || nrow(res) == 0L) next
      dt <- as.data.table(res)
      dt[, id_plot  := id_str]
      dt[, bit_code := bit]
      out[[length(out) + 1L]] <- dt[, .(id_plot, date, bit_code, Delta_Tmax)]
    }
  }
  DT <- rbindlist(out, fill = TRUE)
  setkey(DT, id_plot, bit_code, date)

  dir.create(dirname(cache_path), recursive = TRUE, showWarnings = FALSE)
  saveRDS(DT, cache_path)
  cli_alert_success("Cached {.path {cache_path}} ({.val {nrow(DT)}} rows)")
  DT
}

# ---- 2. Aggregate + wide pivot (mirror Modules 2/3) --------------------------
lovb_contrib_hobo <- function(DT_daily,
                               cache_path = here::here("outputs/lovb/data/DT_contrib_HOBO.rds"),
                               use_cache  = TRUE) {
  if (use_cache && file.exists(cache_path)) {
    return(as.data.table(readRDS(cache_path)))
  }
  cli_h1("HOBO LOVB : aggregate + contribute")

  # Aggregate mean + P90
  DT_agg <- DT_daily[
    , .(Tmax_mean = mean(Delta_Tmax, na.rm = TRUE),
        Tmax_P90  = quantile(Delta_Tmax, 0.9, type = 7, na.rm = TRUE),
        n_days    = sum(!is.na(Delta_Tmax))),
    by = .(id_plot, bit_code)
  ]

  # Pivot wide
  DT_w <- dcast(DT_agg, id_plot ~ bit_code,
                value.var = c("Tmax_mean", "Tmax_P90"))

  # Rename bit codes -> labels (same convention as Module 3)
  bit2lab <- c("0000" = "NULL", "0001" = "LVA_LAD", "0010" = "LVA_fCover",
                "0100" = "LVA_Hmax", "1000" = "LVA_LAI", "0111" = "LOVB_LAI",
                "1011" = "LOVB_Hmax", "1101" = "LOVB_fCover",
                "1110" = "LOVB_LAD", "1111" = "REF")
  for (b in names(bit2lab)) {
    for (m in c("Tmax_mean", "Tmax_P90")) {
      old <- paste0(m, "_", b); new <- paste0(m, "_", bit2lab[b])
      if (old %in% names(DT_w)) setnames(DT_w, old, new)
    }
  }

  # Derived metrics per variable v
  vars <- c("LAI", "Hmax", "fCover", "LAD")
  for (v in vars) {
    DT_w[, paste0("Delta_", v, "_mean") :=
            get("Tmax_mean_REF") - get(paste0("Tmax_mean_LOVB_", v))]
    DT_w[, paste0("Delta_", v, "_P90") :=
            get("Tmax_P90_REF") - get(paste0("Tmax_P90_LOVB_", v))]
  }

  saveRDS(DT_w, cache_path)
  cli_alert_success("Cached {.path {cache_path}} ({.val {nrow(DT_w)}} rows)")
  DT_w
}

# ---- 3. Observed HOBO buffering sensitivity ----------------------------------
# For each HOBO plot, compute slope(Tmax_obs_subcanopy ~ Tmax_macro). This is
# the buffering sensitivity (~1 = no buffering, <1 = buffered). Then regress
# this slope across plots against LAI/Hmax/fCover/FPC1 to get observed coeffs.
#' Project HOBO LAD profiles onto the cLHS FPCA basis -> get FPC1 scores
.lovb_hobo_fpc1 <- function(df_hobo_inputs, df_sample) {
  cli_h2("Project HOBO LAD on cLHS FPCA basis")
  lad_cols <- grep("^LAD_Layer_", names(df_hobo_inputs), value = TRUE)
  z_breaks <- as.numeric(gsub("LAD_Layer_", "", lad_cols))

  # Same LAD columns in df_sample (cLHS) — check + reorder
  lad_cols_s <- grep("^LAD_Layer_", names(df_sample), value = TRUE)
  z_breaks_s <- as.numeric(gsub("LAD_Layer_", "", lad_cols_s))
  common_z <- intersect(z_breaks, z_breaks_s)
  if (length(common_z) < 5)
    stop("Less than 5 common LAD layers between HOBO and cLHS")
  common_lad <- paste0("LAD_Layer_", common_z)

  mat_clhs <- as.matrix(df_sample[, common_lad])
  mat_clhs[is.na(mat_clhs)] <- 0
  mat_hobo <- as.matrix(df_hobo_inputs[, common_lad])
  mat_hobo[is.na(mat_hobo)] <- 0
  cli_alert("LAD matrix : cLHS {nrow(mat_clhs)} rows, HOBO {nrow(mat_hobo)} rows, {ncol(mat_clhs)} layers")

  hmax_clhs <- as.numeric(df_sample$Hmax)
  hmax_hobo <- as.numeric(df_hobo_inputs$Hmax)

  # Combined FPCA : same harmonics for both, HOBO scores read from bottom rows.
  # Done in a single compute_fpca() call on the concatenated matrix to ensure
  # the SAME basis/eigenfunctions for projection.
  mat_all  <- rbind(mat_clhs, mat_hobo)
  hmax_all <- c(hmax_clhs, hmax_hobo)
  fpca_all <- compute_fpca(mat_all, hmax_all, common_z, n_harm = 3L)

  scores      <- as.data.frame(fpca_all$fpca$scores)
  names(scores) <- paste0("FPC", seq_len(ncol(scores)))
  n_clhs <- nrow(mat_clhs)
  hobo_scores <- scores[(n_clhs + 1L):nrow(scores), , drop = FALSE]
  cli_alert_success("FPC1 projected to {nrow(hobo_scores)} HOBO plots — range [{round(min(hobo_scores$FPC1),3)} ; {round(max(hobo_scores$FPC1),3)}]")

  data.table(
    id_plot = df_hobo_inputs$id_plot,
    FPC1    = hobo_scores$FPC1,
    FPC2    = hobo_scores$FPC2,
    FPC3    = hobo_scores$FPC3
  )
}

lovb_observed_buffering <- function() {
  cli_h2("Compute observed buffering sensitivity from HOBO")
  rasters        <- load_lidar_rasters(CFG$in_dir, CFG$agg_factor)
  df_hobo_inputs <- build_hobo_inputs(CFG$hobo_geojson, rasters$stack, CFG$ids_to_remove)
  df_macro       <- extract_macro_daily(CFG$forcing_file, CFG$date_seq)
  df_hobo_daily  <- read_hobo_daily(CFG$hobo_temp_csv, CFG$date_seq,
                                     df_macro, CFG$ids_to_remove)
  df_sample      <- readRDS(here::here("out_files/Sensitivity_Analysis/clhs_sample.rds"))
  if ("Archetype" %in% names(df_sample) && !"Cluster" %in% names(df_sample)) {
    df_sample <- df_sample %>% dplyr::rename(Cluster = Archetype)
  }

  m <- as.data.table(df_hobo_daily)
  if (!("Tmax_macro" %in% names(m)))
    m <- merge(m, as.data.table(df_macro)[, .(date, Tmax_macro)], by = "date")

  # Slope per plot : Tmax_obs (sub-canopy) ~ Tmax_macro
  slopes <- m[, {
    if (sum(!is.na(Tmax_obs) & !is.na(Tmax_macro)) < 5L) {
      .(slope = NA_real_, intercept = NA_real_)
    } else {
      fit <- tryCatch(lm(Tmax_obs ~ Tmax_macro), error = function(e) NULL)
      if (is.null(fit)) .(slope = NA_real_, intercept = NA_real_)
      else .(slope = unname(coef(fit)[2]), intercept = unname(coef(fit)[1]))
    }
  }, by = id_plot]

  # Traits HOBO scalaires
  traits <- as.data.table(df_hobo_inputs)[, .(id_plot, LAI, Hmax, fCover)]

  # FPC1 par projection sur la base cLHS
  fpc_hobo <- .lovb_hobo_fpc1(df_hobo_inputs, df_sample)

  DT_obs <- merge(slopes, traits,   by = "id_plot")
  DT_obs <- merge(DT_obs, fpc_hobo, by = "id_plot")
  cli_alert("Slopes computed for {sum(!is.na(DT_obs$slope))} / {nrow(DT_obs)} HOBO plots, with FPC1")
  DT_obs
}

# ---- 4. Compare simulated LOVB vs observed sensitivity -----------------------
# Observed sensitivity = partial deviance contribution of each variable in a
# multivariate GAM scat() : slope ~ s(LAI) + s(Hmax) + s(fCover) + s(FPC1).
# Partial deviance = deviance(reduced_v) - deviance(full). Larger = variable
# matters more in the multivariate model. Comparable to LOVB RMSE.
lovb_hobo_convergence <- function(DT_contrib_hobo, DT_obs,
                                   out_csv = here::here("outputs/lovb/tables/tab_HOBO_convergence.csv")) {
  cli_h1("LOVB sim vs observed buffering : multivariate GAM convergence")
  vars      <- c("LAI", "Hmax", "fCover", "LAD")
  trait_for <- c(LAI = "LAI", Hmax = "Hmax", fCover = "fCover", LAD = "FPC1")

  # ---- Simulated LOVB RMSE (already computed pipeline-wise) ----
  rmse_lovb <- list()
  for (m in c("mean", "P90")) {
    rmse_lovb[[m]] <- vapply(vars, function(v) {
      sqrt(mean(DT_contrib_hobo[[paste0("Delta_", v, "_", m)]]^2, na.rm = TRUE))
    }, numeric(1))
    names(rmse_lovb[[m]]) <- vars
  }
  rank_lovb_mean <- rank(-rmse_lovb$mean, ties.method = "average")
  rank_lovb_P90  <- rank(-rmse_lovb$P90,  ties.method = "average")

  # ---- Observed sensitivity : multivariate GAM, partial deviance per var ----
  DT_use <- DT_obs[complete.cases(DT_obs[, .(slope, LAI, Hmax, fCover, FPC1)])]
  cli_alert("GAM fit on {nrow(DT_use)} HOBO plots (complete cases)")
  k_smooth <- 5L
  fit_full <- tryCatch(
    mgcv::gam(slope ~ s(LAI, k = k_smooth) + s(Hmax, k = k_smooth) +
                       s(fCover, k = k_smooth) + s(FPC1, k = k_smooth),
              family = mgcv::scat(), data = DT_use, method = "REML"),
    error = function(e) { cli_alert_warning("Full GAM failed: {e$message}"); NULL })
  if (is.null(fit_full))
    return(invisible(list(table = NULL, tau_mean = NA, tau_P90 = NA)))

  R2_full     <- summary(fit_full)$dev.expl
  dev_full    <- deviance(fit_full)
  cli_alert("GAM full : dev.expl = {round(100*R2_full, 1)}%  | dev = {round(dev_full, 3)}")

  # Reduced models : drop one smooth, re-fit, compare partial deviance
  smooth_terms <- list(
    LAI    = "s(LAI, k = k_smooth)",
    Hmax   = "s(Hmax, k = k_smooth)",
    fCover = "s(fCover, k = k_smooth)",
    LAD    = "s(FPC1, k = k_smooth)"
  )
  partial_dev    <- numeric(length(vars)); names(partial_dev) <- vars
  partial_R2     <- numeric(length(vars)); names(partial_R2)  <- vars
  for (v in vars) {
    terms_keep <- smooth_terms[setdiff(vars, v)]
    fml <- as.formula(paste("slope ~", paste(unlist(terms_keep), collapse = " + ")))
    fit_r <- tryCatch(
      mgcv::gam(fml, family = mgcv::scat(), data = DT_use, method = "REML"),
      error = function(e) NULL)
    if (is.null(fit_r)) {
      partial_dev[v] <- NA_real_; partial_R2[v] <- NA_real_
    } else {
      partial_dev[v] <- deviance(fit_r) - dev_full   # >0 = removing v hurts fit
      partial_R2[v]  <- R2_full - summary(fit_r)$dev.expl
    }
  }

  # Ranking by partial_R2 (= dev.expl_full - dev.expl_reduced). Higher = removing
  # the variable hurts dev.expl more = variable contributes more unique info.
  # NB on small samples (n=53) with REML smooths, partial_dev can be negative
  # due to penalty re-balancing when a smooth is dropped. partial_R2 is the
  # cleaner reading.
  rank_obs <- rank(-partial_R2, ties.method = "average", na.last = "keep")
  rank_obs[is.na(rank_obs)] <- max(rank_obs, na.rm = TRUE) + 1

  DT_out <- data.table(
    Variable          = vars,
    RMSE_LOVB_mean    = rmse_lovb$mean,
    RMSE_LOVB_P90     = rmse_lovb$P90,
    Rank_LOVB_mean    = rank_lovb_mean,
    Rank_LOVB_P90     = rank_lovb_P90,
    Partial_dev_obs   = partial_dev,
    Partial_R2_obs    = partial_R2,
    Rank_obs          = rank_obs
  )

  tau_mean <- suppressWarnings(cor(DT_out$Rank_LOVB_mean, DT_out$Rank_obs,
                                     method = "kendall"))
  tau_P90  <- suppressWarnings(cor(DT_out$Rank_LOVB_P90,  DT_out$Rank_obs,
                                     method = "kendall"))

  dir.create(dirname(out_csv), recursive = TRUE, showWarnings = FALSE)
  fwrite(DT_out, out_csv)

  cat("\n### HOBO convergence — LOVB simulated vs observed (mgcv GAM scat)\n")
  cat(sprintf("Multivariate GAM : slope ~ s(LAI) + s(Hmax) + s(fCover) + s(FPC1)  family = scat()\n"))
  cat(sprintf("Full model dev.expl = %.1f%%   |   n = %d HOBO plots\n", 100 * R2_full, nrow(DT_use)))
  cat("Rank_obs = ranking by partial_R2 (= dev.expl_full - dev.expl_reduced)\n\n")
  print(DT_out)
  cat(sprintf("\nKendall tau (LOVB_mean vs partial_R2_obs) = %.3f\n", tau_mean))
  cat(sprintf("Kendall tau (LOVB_P90  vs partial_R2_obs) = %.3f\n", tau_P90))
  cli_alert_success("Saved {.path {out_csv}}")

  invisible(list(table = DT_out, tau_mean = tau_mean, tau_P90 = tau_P90,
                  R2_full = R2_full, gam_full = fit_full))
}

# ---- 5. Figure : LOVB Delta_v scatter at HOBO + observed slope ---------------
lovb_fig_hobo <- function(DT_contrib_hobo, DT_obs,
                           out_path = here::here("outputs/lovb/figures/fig_LOVB_HOBO_validation.png")) {
  vars <- c("LAI", "Hmax", "fCover", "LAD")
  trait_for <- c(LAI = "LAI", Hmax = "Hmax", fCover = "fCover", LAD = NA)

  # Build a long df : for each var, per-plot Delta_v_mean and observed trait
  long <- rbindlist(lapply(vars, function(v) {
    tr <- trait_for[v]
    trait_vals <- if (is.na(tr)) rep(NA_real_, nrow(DT_contrib_hobo)) else DT_obs[match(DT_contrib_hobo$id_plot, id_plot), get(tr)]
    data.table(
      id_plot    = DT_contrib_hobo$id_plot,
      variable   = v,
      Delta_mean = DT_contrib_hobo[[paste0("Delta_", v, "_mean")]],
      slope_obs  = DT_obs[match(DT_contrib_hobo$id_plot, id_plot), slope],
      trait      = trait_vals
    )
  }))
  long[, variable := factor(variable, levels = vars)]

  long_clean <- long[!is.na(slope_obs) & !is.na(Delta_mean)]
  cli_alert("Figure HOBO : {nrow(long_clean)} valid rows over {nrow(long)} total")

  p <- ggplot(long_clean, aes(x = slope_obs, y = Delta_mean)) +
    geom_hline(yintercept = 0, linetype = "dashed", colour = "grey50") +
    geom_point(alpha = 0.7, size = 2.0, colour = "#31688e")
  if (nrow(long_clean) >= 10) {
    p <- p + geom_smooth(method = "lm", se = TRUE, colour = "grey25", linewidth = 0.8)
  }
  p <- p +
    facet_wrap(~ variable, nrow = 1,
                labeller = labeller(variable = function(x) paste("Delta_", x, " (sim)"))) +
    labs(
      title = "LOVB simulated contribution vs HOBO observed buffering",
      subtitle = "x-axis = slope(Tmax_obs ~ Tmax_macro) per HOBO plot | y-axis = Delta_v_mean from MuSICA",
      x = "Observed buffering slope (HOBO)",
      y = "Delta_v (sim, mean) [degree C]"
    ) +
    theme_bw(base_size = 11) +
    theme(strip.text = element_text(face = "bold"),
          plot.title = element_text(face = "bold"))

  dir.create(dirname(out_path), recursive = TRUE, showWarnings = FALSE)
  ggsave(out_path, p, width = 13, height = 4.5, dpi = 300)
  cli_alert_success("Saved {.path {out_path}}")
  invisible(p)
}

# ---- Main orchestrator for Module 7 ------------------------------------------
lovb_run_hobo_analysis <- function() {
  cli_h1("Module 7 : HOBO LOVB analysis")
  DT_d <- lovb_load_daily_hobo(use_cache = TRUE)
  DT_c <- lovb_contrib_hobo(DT_d, use_cache = FALSE)
  DT_obs <- lovb_observed_buffering()
  res <- lovb_hobo_convergence(DT_c, DT_obs)
  lovb_fig_hobo(DT_c, DT_obs)
  invisible(res)
}

if (sys.nframe() == 0) {
  lovb_run_hobo_analysis()
}
