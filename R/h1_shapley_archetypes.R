# ==============================================================================
# Model-centric Shapley analysis on 4 synthetic archetypes (one per FPCA cluster).
#
# Caracterise la reponse de MuSICA aux traits structuraux d'un archetype typique
# de chaque cluster, independamment de la distribution Blois reelle. Complementaire
# au Shapley cLHS principal (qui mesure "Blois moyen") en isolant l'effet
# purement structural sans variance inter-plot intra-cluster.
#
# Baselines identiques au pipeline cLHS :
#   LAI_b    = mean(df_sample$LAI)        (~3.13)
#   Hmax_b   = round(mean(floor(Hmax)+1)) (~21)
#   fCover_b = 1 (absent)
#   LAD_b    = uniform
# ==============================================================================

#' Construct 4 synthetic archetypes (1 per cluster) from the cLHS sample.
#'
#' @param df_sample cLHS sample with cluster labels, LAD_Layer_* columns.
#' @return data.frame with one row per cluster, ready for run_musica_one().
make_synthetic_archetypes <- function(df_sample) {
  lad_cols <- grep("^LAD_Layer_", names(df_sample), value = TRUE)
  if (length(lad_cols) == 0)
    stop("make_synthetic_archetypes: no LAD_Layer_* columns in df_sample")

  clusters <- sort(unique(as.character(df_sample$Cluster)))
  cat(sprintf("  [archetypes] Building synthetic archetypes for %d clusters\n",
              length(clusters)))

  rows <- lapply(seq_along(clusters), function(i) {
    cl    <- clusters[i]
    df_cl <- df_sample[as.character(df_sample$Cluster) == cl, , drop = FALSE]
    mat   <- as.matrix(df_cl[, lad_cols, drop = FALSE])
    mat[is.na(mat)] <- 0

    lad_mean <- colMeans(mat, na.rm = TRUE)
    lad_df   <- as.data.frame(t(lad_mean))
    names(lad_df) <- lad_cols

    # Coords synthetiques uniques (negatif = sentinelle archetype)
    out <- data.frame(
      x       = -1000 - i,
      y       = -1000 - i,
      Cluster = cl,
      LAI     = mean(df_cl$LAI,    na.rm = TRUE),
      Hmax    = mean(df_cl$Hmax,   na.rm = TRUE),
      fCover  = mean(df_cl$fCover, na.rm = TRUE),
      stringsAsFactors = FALSE
    )
    cbind(out, lad_df)
  })
  df_arch <- do.call(rbind, rows)

  cat("  [archetypes] Resume des 4 archetypes :\n")
  print(df_arch[, c("Cluster", "LAI", "Hmax", "fCover")])

  df_arch
}

#' Build the 16 factorial scenarios applied to archetypes, using baselines from df_sample.
#'
#' @param df_sample cLHS sample dataframe with LAI / Hmax / fCover columns.
#' @param fcov_b Baseline fCover value when bit=0 ("absent"). Default 1
#'   (closed canopy reference). Set to "mean" or a numeric ∈ (0,1] to use a
#'   statistical baseline symmetric with LAI/Hmax (Eva's suggestion).
build_factorial_scenarios_archetypes <- function(df_sample, fcov_b = 1) {
  # Baselines (identiques a scenarios_h1_factorial)
  lai_b   <- mean(df_sample$LAI, na.rm = TRUE)
  hmax_b  <- round(mean(floor(df_sample$Hmax) + 1, na.rm = TRUE))
  if (identical(fcov_b, "mean")) fcov_b <- mean(df_sample$fCover, na.rm = TRUE)
  fcov_b  <- as.numeric(fcov_b)
  cat(sprintf("  [archetypes] Baselines : LAI_b=%.3f | Hmax_b=%d | fCover_b=%g | LAD_b=uniform\n",
              lai_b, hmax_b, fcov_b))

  grid <- expand.grid(
    LAI    = c("real", "mean"),
    Hmax   = c("real", "mean"),
    fCover = c("real", "absent"),
    LAD    = c("real", "uniform"),
    stringsAsFactors = FALSE
  )

  lapply(seq_len(nrow(grid)), function(i) {
    g  <- grid[i, ]
    nm <- sprintf("ARCH_%s_%s_%s_%s",
                  substr(g$LAI, 1, 1), substr(g$Hmax, 1, 1),
                  substr(g$fCover, 1, 1), substr(g$LAD, 1, 1))

    lai_fn    <- if (g$LAI    == "real") function(pr) as.numeric(pr$LAI)
                 else                     function(pr) lai_b
    hmax_fn   <- if (g$Hmax   == "real") function(pr) as.numeric(pr$Hmax)
                 else                     function(pr) hmax_b
    fcover_fn <- if (g$fCover == "real") function(pr) as.numeric(pr$fCover)
                 else                     function(pr) fcov_b
    lad_fn    <- if (g$LAD    == "real") make_lad_real
                 else                     make_lad_uniform

    list(name = nm, grid_row = g,
         lai_fn = lai_fn, hmax_fn = hmax_fn,
         fcover_fn = fcover_fn, lad_fn = lad_fn)
  }) -> fac_scs
  names(fac_scs) <- sapply(fac_scs, `[[`, "name")
  fac_scs
}

#' Run the 4 archetypes x 16 coalitions = 64 MuSICA simulations.
run_archetype_factorial <- function(df_archetypes, fac_scs, out_root,
                                     forcing_file, musica_cmd,
                                     musica_nml, musica_variables) {
  dir.create(out_root, recursive = TRUE, showWarnings = FALSE)
  for (i in seq_len(nrow(df_archetypes))) {
    arch     <- df_archetypes[i, , drop = FALSE]
    arch_lbl <- sprintf("Arch_C%s", arch$Cluster)
    out_dir  <- file.path(out_root, arch_lbl)
    dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
    cat(sprintf("\n[archetype %s] LAI=%.2f Hmax=%.1f fCover=%.2f\n",
                arch_lbl, arch$LAI, arch$Hmax, arch$fCover))

    for (sc_nm in names(fac_scs)) {
      sc     <- fac_scs[[sc_nm]]
      out_nc <- file.path(out_dir, sprintf("musica_out_%s_%s.nc", arch_lbl, sc_nm))
      if (file.exists(out_nc)) next
      cat(sprintf("  [%s] %s\n", arch_lbl, sc_nm))
      run_musica_one(arch, sc, out_nc, forcing_file, musica_cmd,
                      musica_nml       = musica_nml,
                      musica_variables = musica_variables)
    }
  }
  invisible(out_root)
}

#' Extract daily ΔTmax per (archetype, scenario, day).
extract_archetype_deltas <- function(out_root, df_macro, date_seq) {
  archetypes <- list.dirs(out_root, recursive = FALSE, full.names = FALSE)
  out <- list()
  for (arch_lbl in archetypes) {
    arch_dir <- file.path(out_root, arch_lbl)
    nc_files <- list.files(arch_dir, pattern = "\\.nc$", full.names = TRUE)
    nc_files <- nc_files[file.size(nc_files) > 1e6]   # skip truncated
    for (f in nc_files) {
      sc_nm <- sub(sprintf("^musica_out_%s_(.*)\\.nc$", arch_lbl), "\\1", basename(f))
      res   <- extract_deltatmax_one(f, df_macro, date_seq)
      if (is.null(res) || nrow(res) == 0) next
      res$archetype <- arch_lbl
      res$scenario  <- sc_nm
      out[[length(out) + 1]] <- res
    }
  }
  do.call(rbind, out)
}

# Coalition mapping (ARCH_<L>_<H>_<F>_<D>) -> 4-bit code (LAI|Hmax|fCover|LAD)
.ARCH_COAL_MAP <- c(
  "ARCH_m_m_a_u"="0000", "ARCH_r_m_a_u"="1000",
  "ARCH_m_r_a_u"="0100", "ARCH_m_m_r_u"="0010",
  "ARCH_m_m_a_r"="0001", "ARCH_r_r_a_u"="1100",
  "ARCH_r_m_r_u"="1010", "ARCH_r_m_a_r"="1001",
  "ARCH_m_r_r_u"="0110", "ARCH_m_r_a_r"="0101",
  "ARCH_m_m_r_r"="0011", "ARCH_r_r_r_u"="1110",
  "ARCH_r_r_a_r"="1101", "ARCH_r_m_r_r"="1011",
  "ARCH_m_r_r_r"="0111", "ARCH_r_r_r_r"="1111"
)

#' Compute RMSE per scenario for one archetype (vs its full-real reference).
.score_one_archetype <- function(df_a, ref_name = "ARCH_r_r_r_r") {
  df_ref <- df_a[df_a$scenario == ref_name, c("date", "Delta_Tmax")]
  if (nrow(df_ref) == 0)
    stop(sprintf(".score_one_archetype: reference '%s' missing", ref_name))
  names(df_ref)[2] <- "Delta_ref"
  m <- merge(df_a, df_ref, by = "date")
  m$err <- m$Delta_Tmax - m$Delta_ref
  agg <- aggregate(m$err,
                    by = list(scenario = m$scenario),
                    FUN = function(x) sqrt(mean(x^2, na.rm = TRUE)))
  names(agg)[2] <- "rmse"
  agg
}

#' Compute exact Shapley per archetype (point estimate).
compute_archetype_shapley <- function(df_delta) {
  if (is.null(df_delta) || nrow(df_delta) == 0)
    stop("compute_archetype_shapley: df_delta is empty — MuSICA likely failed for all archetypes.")
  archs <- unique(df_delta$archetype)
  variables <- c("LAI", "Hmax", "fCover", "LAD")

  out <- list()
  for (arch in archs) {
    df_a <- df_delta[df_delta$archetype == arch, ]

    scores  <- .score_one_archetype(df_a)
    rmse_v  <- vapply(names(.ARCH_COAL_MAP), function(sc_nm) {
      val <- scores$rmse[scores$scenario == sc_nm]
      if (length(val) == 0) NA_real_ else val[1]
    }, numeric(1))
    names(rmse_v) <- unname(.ARCH_COAL_MAP)

    if (any(is.na(rmse_v))) {
      warning(sprintf("[archetype %s] missing scenarios, skipping", arch))
      next
    }
    attr(rmse_v, "metric") <- "rmse"

    capture.output(shap <- compute_shapley_exact(rmse_v))

    out[[arch]] <- data.frame(
      archetype  = arch,
      variable   = variables,
      phi        = as.numeric(shap$shapley),
      pct        = 100 * as.numeric(shap$shapley) / shap$total_effect,
      v_N        = shap$total_effect,
      check_diff = shap$check_diff,
      stringsAsFactors = FALSE
    )
  }
  do.call(rbind, out)
}

#' Compute pooled Shapley on the 4 archetypes treated as one mini-sample.
#'
#' Different from mean_archetypes (mean of per-archetype phi) :
#' here we POOL all (archetype, day) squared errors first, then compute one
#' RMSE per coalition, then one Shapley. Non-linear interaction of errors
#' before aggregation can change the answer vs averaging phi after the fact.
compute_pooled_archetype_shapley <- function(df_delta) {
  if (is.null(df_delta) || nrow(df_delta) == 0)
    stop("compute_pooled_archetype_shapley: df_delta is empty.")
  variables <- c("LAI", "Hmax", "fCover", "LAD")

  df_ref <- df_delta[df_delta$scenario == "ARCH_r_r_r_r",
                       c("archetype", "date", "Delta_Tmax")]
  if (nrow(df_ref) == 0)
    stop("compute_pooled_archetype_shapley: reference scenario 'ARCH_r_r_r_r' absent — vérifier que MuSICA a tourné pour ce scénario.")
  names(df_ref)[3] <- "Delta_ref"
  m <- merge(df_delta, df_ref, by = c("archetype", "date"))
  m$err <- m$Delta_Tmax - m$Delta_ref

  # Pooled RMSE across all (archetype, day)
  agg <- aggregate(m$err, by = list(scenario = m$scenario),
                    FUN = function(x) sqrt(mean(x^2, na.rm = TRUE)))
  names(agg)[2] <- "rmse"

  rmse_v <- vapply(names(.ARCH_COAL_MAP), function(sc_nm) {
    val <- agg$rmse[agg$scenario == sc_nm]
    if (length(val) == 0) NA_real_ else val[1]
  }, numeric(1))
  names(rmse_v) <- unname(.ARCH_COAL_MAP)
  if (any(is.na(rmse_v)))
    stop("compute_pooled_archetype_shapley: missing scenarios")
  attr(rmse_v, "metric") <- "rmse"

  capture.output(shap <- compute_shapley_exact(rmse_v))

  data.frame(
    archetype  = "pooled_4archetypes",
    variable   = variables,
    phi        = as.numeric(shap$shapley),
    pct        = 100 * as.numeric(shap$shapley) / shap$total_effect,
    v_N        = shap$total_effect,
    check_diff = shap$check_diff,
    stringsAsFactors = FALSE
  )
}

#' Bootstrap pooled archetype Shapley by resampling shared DAYS.
#'
#' Resample 122 days (with replacement), apply SAME resample to all 4
#' archetypes (preserves alignment), recompute pooled RMSE & Shapley.
bootstrap_pooled_archetype_shapley_temporal <- function(df_delta, n_boot = 200) {
  variables <- c("LAI", "Hmax", "fCover", "LAD")
  set.seed(46)

  df_ref <- df_delta[df_delta$scenario == "ARCH_r_r_r_r",
                       c("archetype", "date", "Delta_Tmax")]
  names(df_ref)[3] <- "Delta_ref"
  m <- merge(df_delta, df_ref, by = c("archetype", "date"))
  m$sq_err <- (m$Delta_Tmax - m$Delta_ref)^2

  dates <- unique(m$date)
  n_d   <- length(dates)
  cat(sprintf("  [pooled archetypes] bootstrap %d iter on %d shared days...\n",
              n_boot, n_d))

  boot_recs <- list()
  for (b in seq_len(n_boot)) {
    samp <- sample(dates, n_d, replace = TRUE)
    wt   <- as.data.frame(table(samp), stringsAsFactors = FALSE)
    names(wt) <- c("date", "weight")
    wt$date <- as.Date(wt$date)

    mb <- merge(m, wt, by = "date")
    agg <- aggregate(cbind(mb$weight * mb$sq_err, mb$weight),
                      by = list(scenario = mb$scenario),
                      FUN = sum)
    names(agg) <- c("scenario", "ws", "w")
    agg$rmse <- sqrt(agg$ws / agg$w)

    rmse_v <- vapply(names(.ARCH_COAL_MAP), function(sc_nm) {
      val <- agg$rmse[agg$scenario == sc_nm]
      if (length(val) == 0) NA_real_ else val[1]
    }, numeric(1))
    names(rmse_v) <- unname(.ARCH_COAL_MAP)
    if (any(is.na(rmse_v))) next
    attr(rmse_v, "metric") <- "rmse"

    res <- tryCatch({
      capture.output(sb <- compute_shapley_exact(rmse_v))
      data.frame(boot_id = b, archetype = "pooled_4archetypes",
                  variable = variables,
                  phi = as.numeric(sb$shapley),
                  stringsAsFactors = FALSE)
    }, error = function(e) NULL)
    if (!is.null(res)) boot_recs[[b]] <- res
    if (b %% 50 == 0) cat(sprintf("    pooled boot %d/%d\n", b, n_boot))
  }
  do.call(rbind, boot_recs)
}

#' Bootstrap Shapley per archetype by resampling DAYS (with replacement).
bootstrap_archetype_shapley_temporal <- function(df_delta, n_boot = 200) {
  archs     <- unique(df_delta$archetype)
  variables <- c("LAI", "Hmax", "fCover", "LAD")
  set.seed(45)

  out <- list()
  for (arch in archs) {
    df_a   <- df_delta[df_delta$archetype == arch, ]
    df_ref <- df_a[df_a$scenario == "ARCH_r_r_r_r", c("date", "Delta_Tmax")]
    names(df_ref)[2] <- "Delta_ref"
    m <- merge(df_a, df_ref, by = "date")
    m$sq_err <- (m$Delta_Tmax - m$Delta_ref)^2

    dates <- unique(m$date)
    n_d   <- length(dates)

    cat(sprintf("  [archetype %s] bootstrap %d iter on %d days...\n",
                arch, n_boot, n_d))

    boot_recs <- list()
    for (b in seq_len(n_boot)) {
      samp <- sample(dates, n_d, replace = TRUE)
      wt   <- as.data.frame(table(samp), stringsAsFactors = FALSE)
      names(wt) <- c("date", "weight")
      wt$date <- as.Date(wt$date)

      mb <- merge(m, wt, by = "date")
      agg <- aggregate(cbind(mb$weight * mb$sq_err, mb$weight),
                        by = list(scenario = mb$scenario),
                        FUN = sum)
      names(agg) <- c("scenario", "ws", "w")
      agg$rmse <- sqrt(agg$ws / agg$w)

      rmse_v <- vapply(names(.ARCH_COAL_MAP), function(sc_nm) {
        val <- agg$rmse[agg$scenario == sc_nm]
        if (length(val) == 0) NA_real_ else val[1]
      }, numeric(1))
      names(rmse_v) <- unname(.ARCH_COAL_MAP)
      if (any(is.na(rmse_v))) next
      attr(rmse_v, "metric") <- "rmse"

      res <- tryCatch({
        capture.output(sb <- compute_shapley_exact(rmse_v))
        data.frame(boot_id = b, variable = variables,
                    phi = as.numeric(sb$shapley),
                    stringsAsFactors = FALSE)
      }, error = function(e) NULL)
      if (!is.null(res)) {
        res$archetype <- arch
        boot_recs[[b]] <- res
      }
      if (b %% 50 == 0) cat(sprintf("    boot %d/%d\n", b, n_boot))
    }
    out[[arch]] <- do.call(rbind, boot_recs)
  }
  do.call(rbind, out)
}

# ---- Plots --------------------------------------------------------------------

#' Barplot per archetype (4 panels, 4 bars Shapley + bootstrap CI).
plot_shapley_archetypes_barplot <- function(shap_point, shap_boot, out_path) {
  variables <- c("LAI", "Hmax", "fCover", "LAD")

  ci <- shap_boot %>%
    dplyr::group_by(archetype, variable) %>%
    dplyr::summarise(
      ci_lo = quantile(phi, 0.025, na.rm = TRUE),
      ci_hi = quantile(phi, 0.975, na.rm = TRUE),
      .groups = "drop"
    )

  df <- merge(shap_point, ci, by = c("archetype", "variable"))
  df$variable <- factor(df$variable, levels = variables)
  # Order : Arch_C1, Arch_C2, ..., then pooled_4archetypes at the end
  arch_lvls  <- sort(grep("^Arch_C", unique(df$archetype), value = TRUE))
  other_lvls <- setdiff(unique(df$archetype), arch_lvls)
  df$archetype <- factor(df$archetype, levels = c(arch_lvls, other_lvls))
  pal <- c(LAI = "#440154", Hmax = "#35b779",
            fCover = "#31688e", LAD = "#fde725")

  p <- ggplot(df, aes(x = variable, y = phi, fill = variable)) +
    geom_col(width = 0.7, colour = "grey20", linewidth = 0.3) +
    geom_hline(yintercept = 0, colour = "grey40", linewidth = 0.4) +
    geom_errorbar(aes(ymin = ci_lo, ymax = ci_hi),
                   width = 0.18, linewidth = 0.7) +
    geom_text(aes(label = sprintf("%+.3f", phi),
                   vjust = ifelse(phi >= 0, -1.0, 1.8)),
              size = 3.3, fontface = "bold") +
    facet_wrap(~ archetype, nrow = 1) +
    scale_fill_manual(values = pal, guide = "none") +
    labs(
      title    = "Shapley attribution per synthetic archetype — model-centric view",
      subtitle = paste0("Each archetype = synthetic mean plot per cluster. ",
                        "Error bars = 95% CI from temporal bootstrap (n=200, days)."),
      caption  = "Temporal bootstrap measures sensitivity to seasonal sub-samples for THIS archetype — not structural uncertainty across Blois.",
      x = NULL, y = expression(phi ~ "(" * degree * "C)")
    ) +
    theme_bw(base_size = 11) +
    theme(strip.text   = element_text(face = "bold"),
          plot.title   = element_text(face = "bold"),
          plot.caption = element_text(size = 8, colour = "grey35", hjust = 0,
                                       lineheight = 1.2))

  dir.create(dirname(out_path), recursive = TRUE, showWarnings = FALSE)
  save_plot(p, out_path, width = 13, height = 6)
  cat(sprintf("  [archetype Shapley] PNG : %s\n", out_path))
  invisible(p)
}

#' Heatmap comparant Shapley cLHS vs Shapley archetype par cluster.
plot_shapley_archetypes_comparison <- function(shap_point, df_clhs, out_path) {
  variables <- c("LAI", "Hmax", "fCover", "LAD")

  # cLHS values
  df_clhs_long <- data.frame(
    variable = df_clhs$variable,
    source   = "cLHS_400",
    phi      = df_clhs$shapley,
    stringsAsFactors = FALSE
  )

  # Per-cluster archetypes (rename Arch_C1 -> C1) + pooled separately
  arch_only  <- shap_point[grepl("^Arch_C", shap_point$archetype), ]
  pooled_row <- shap_point[shap_point$archetype == "pooled_4archetypes", ]

  arch_only$source <- sub("^Arch_", "", arch_only$archetype)
  df_arch_long <- arch_only[, c("variable", "source", "phi")]

  # Mean across the 4 per-cluster archetypes (NOT including pooled)
  df_mean <- aggregate(arch_only$phi,
                        by = list(variable = arch_only$variable),
                        FUN = mean, na.rm = TRUE)
  names(df_mean)[2] <- "phi"
  df_mean$source <- "mean_archetypes"
  df_mean <- df_mean[, c("variable", "source", "phi")]

  # Pooled archetype as its own column
  df_pooled <- if (nrow(pooled_row) > 0) {
    pooled_row$source <- "pooled_4archetypes"
    pooled_row[, c("variable", "source", "phi")]
  } else NULL

  df_all <- rbind(df_clhs_long, df_arch_long, df_pooled, df_mean)
  df_all$variable <- factor(df_all$variable, levels = variables)
  source_levels <- c("cLHS_400",
                      sort(unique(df_arch_long$source)),
                      "pooled_4archetypes",
                      "mean_archetypes")
  source_levels <- intersect(source_levels, unique(df_all$source))
  df_all$source <- factor(df_all$source, levels = source_levels)

  # symmetric colour scale
  rng <- max(abs(df_all$phi), na.rm = TRUE)

  p <- ggplot(df_all, aes(x = source, y = variable, fill = phi)) +
    geom_tile(colour = "white", linewidth = 0.7) +
    geom_text(aes(label = sprintf("%+.3f", phi)),
              size = 3.6, fontface = "bold") +
    scale_fill_gradient2(low = "#31688e", mid = "white", high = "#d8576b",
                         midpoint = 0, limits = c(-rng, rng),
                         name = expression(phi ~ "(" * degree * "C)")) +
    labs(
      title    = "Shapley comparison — cLHS sample vs synthetic archetypes",
      subtitle = "Rows : structural variables | Cols : analysis variant. Sign and magnitude reveal whether the hierarchy is universal or cluster-dependent.",
      x = NULL, y = NULL
    ) +
    theme_bw(base_size = 12) +
    theme(plot.title    = element_text(face = "bold"),
          axis.text.x   = element_text(face = "bold"),
          axis.text.y   = element_text(face = "bold"),
          panel.grid    = element_blank())

  dir.create(dirname(out_path), recursive = TRUE, showWarnings = FALSE)
  save_plot(p, out_path, width = 11, height = 6)
  cat(sprintf("  [archetype comparison] PNG : %s\n", out_path))
  invisible(p)
}
