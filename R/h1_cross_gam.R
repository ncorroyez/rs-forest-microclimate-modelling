# ==============================================================================
# H1 cross-trained GAM : compare functional response to LiDAR traits
# between MuSICA simulations (400 cLHS plots) and HOBO observations (53 sensors).
#
# Target dTmax_mean (mean over 122 summer days 2021) :
#   - cLHS  : MuSICA simulated, bit=1111  ->  DT_contrib_cLHS$Tmax_mean_REF
#   - HOBO  : sensor reading -> mean(Tmax_macro - Tmax_obs) per sensor
#
# Same formula on both : dTmax_mean ~ s(LAI,k=5) + s(Hmax,k=5) +
#                                      s(fCover,k=5) + s(FPC1,k=5)
# family = scat()  method = "REML"
# ==============================================================================

suppressMessages({
  library(data.table)
  library(here)
  library(cli)
  library(mgcv)
  library(ncdf4)
  library(tidyverse)   # used only by sourced legacy helpers
  library(sf); library(terra)
  library(musica.tools); library(rmusica)
})

# Sourced helpers (already exist) -----------------------------------------------
source(here::here("R/config.R"))
source(here::here("R/io.R"))
source(here::here("R/lad.R"))
source(here::here("R/forest.R"))
source(here::here("R/scenarios.R"))
source(here::here("R/musica.R"))
source(here::here("R/validation.R"))
source(here::here("R/fpca.R"))
source(here::here("R/lovb_07_hobo_analysis.R"))    # .lovb_hobo_fpc1

# ---- Helper : joint FPCA returning FPC1/2/3 for both cLHS and HOBO -----------
.joint_fpca_scores <- function(df_hobo_inputs, df_sample, n_harm = 3L) {
  cli_h2("Joint FPCA on cLHS (400) + HOBO (53) LAD profiles")
  lad_cols   <- grep("^LAD_Layer_", names(df_hobo_inputs), value = TRUE)
  lad_cols_s <- grep("^LAD_Layer_", names(df_sample),      value = TRUE)
  z_breaks   <- as.numeric(gsub("LAD_Layer_", "", lad_cols))
  z_breaks_s <- as.numeric(gsub("LAD_Layer_", "", lad_cols_s))
  common_z   <- intersect(z_breaks, z_breaks_s)
  common_lad <- paste0("LAD_Layer_", common_z)

  mat_clhs <- as.matrix(df_sample[, common_lad])
  mat_clhs[is.na(mat_clhs)] <- 0
  mat_hobo <- as.matrix(df_hobo_inputs[, common_lad])
  mat_hobo[is.na(mat_hobo)] <- 0
  hmax_all <- c(as.numeric(df_sample$Hmax), as.numeric(df_hobo_inputs$Hmax))
  mat_all  <- rbind(mat_clhs, mat_hobo)
  fpca_all <- compute_fpca(mat_all, hmax_all, common_z, n_harm = n_harm)

  scores <- as.data.frame(fpca_all$fpca$scores)
  names(scores) <- paste0("FPC", seq_len(ncol(scores)))
  n_c <- nrow(mat_clhs)
  list(
    cLHS = data.table(x = df_sample$x, y = df_sample$y,
                       FPC1 = scores$FPC1[1:n_c],
                       FPC2 = scores$FPC2[1:n_c],
                       FPC3 = scores$FPC3[1:n_c]),
    HOBO = data.table(id_plot = df_hobo_inputs$id_plot,
                       FPC1 = scores$FPC1[(n_c + 1L):nrow(scores)],
                       FPC2 = scores$FPC2[(n_c + 1L):nrow(scores)],
                       FPC3 = scores$FPC3[(n_c + 1L):nrow(scores)])
  )
}

# ---- Step 1 : aligned targets ------------------------------------------------
load_cross_targets <- function() {
  cli_h1("Cross-GAM step 1 — aligned dTmax_mean targets")

  # HOBO + cLHS inputs needed for joint FPCA
  rasters        <- load_lidar_rasters(CFG$in_dir, CFG$agg_factor)
  df_hobo_inputs <- build_hobo_inputs(CFG$hobo_geojson, rasters$stack,
                                       CFG$ids_to_remove)
  df_macro       <- extract_macro_daily(CFG$forcing_file, CFG$date_seq)
  df_hobo_daily  <- read_hobo_daily(CFG$hobo_temp_csv, CFG$date_seq,
                                     df_macro, CFG$ids_to_remove)
  df_sample      <- readRDS(here::here("out_files/Sensitivity_Analysis/clhs_sample.rds"))
  if ("Archetype" %in% names(df_sample) && !"Cluster" %in% names(df_sample))
    df_sample <- df_sample %>% dplyr::rename(Cluster = Archetype)

  # Joint FPCA — same basis for cLHS and HOBO, returns FPC1/2/3 both sides
  fpca_scores <- .joint_fpca_scores(df_hobo_inputs, df_sample)

  # cLHS side : target already cached, attach joint FPC1/2/3 via (x,y)
  DT_clhs <- readRDS(here::here("outputs/lovb/data/DT_contrib_cLHS.rds"))
  DT_clhs <- DT_clhs[, .(x, y, LAI, Hmax, fCover, dTmax_mean = Tmax_mean_REF)]
  DT_clhs <- merge(DT_clhs, fpca_scores$cLHS, by = c("x", "y"), all.x = TRUE)
  DT_clhs[, id := paste0("X", x, "_Y", y)]
  DT_clhs <- DT_clhs[, .(id, LAI, Hmax, fCover, FPC1, FPC2, FPC3, dTmax_mean)]
  cli_alert("cLHS  : {nrow(DT_clhs)} plots  |  dTmax_mean from Tmax_mean_REF (MuSICA bit=1111)")

  # Observed dTmax_mean = mean(Tmax_macro - Tmax_obs) over 122 d per sensor
  m <- as.data.table(df_hobo_daily)
  if (!("Tmax_macro" %in% names(m)))
    m <- merge(m, as.data.table(df_macro)[, .(date, Tmax_macro)], by = "date")
  DT_target <- m[, .(dTmax_mean = mean(Tmax_macro - Tmax_obs, na.rm = TRUE),
                      n_days     = sum(!is.na(Tmax_obs) & !is.na(Tmax_macro))),
                  by = id_plot]
  traits <- as.data.table(df_hobo_inputs)[, .(id_plot, LAI, Hmax, fCover)]
  DT_hobo <- merge(merge(traits, fpca_scores$HOBO, by = "id_plot"),
                    DT_target, by = "id_plot")
  DT_hobo[, id := id_plot]
  DT_hobo <- DT_hobo[, .(id, LAI, Hmax, fCover, FPC1, FPC2, FPC3, dTmax_mean, n_days)]

  cli_alert("HOBO  : {nrow(DT_hobo)} sensors |  dTmax_mean = mean(Tmax_macro - Tmax_obs) over n_days")
  cli_alert("        n_days range : [{min(DT_hobo$n_days)} ; {max(DT_hobo$n_days)}]")
  cli_alert("FPC1 range  cLHS [{round(min(DT_clhs$FPC1),3)} ; {round(max(DT_clhs$FPC1),3)}]  HOBO [{round(min(DT_hobo$FPC1),3)} ; {round(max(DT_hobo$FPC1),3)}]")
  cli_alert("FPC2 range  cLHS [{round(min(DT_clhs$FPC2),3)} ; {round(max(DT_clhs$FPC2),3)}]  HOBO [{round(min(DT_hobo$FPC2),3)} ; {round(max(DT_hobo$FPC2),3)}]")
  cli_alert("FPC3 range  cLHS [{round(min(DT_clhs$FPC3),3)} ; {round(max(DT_clhs$FPC3),3)}]  HOBO [{round(min(DT_hobo$FPC3),3)} ; {round(max(DT_hobo$FPC3),3)}]")

  list(cLHS = DT_clhs,
        HOBO = DT_hobo[, .(id, LAI, Hmax, fCover, FPC1, FPC2, FPC3, dTmax_mean)])
}

# ---- Step 2 : symmetric GAM fits ---------------------------------------------
# FPC1 replaced by ti(FPC1, FPC2, FPC3, k=3) = 3D LAD shape interaction.
# k=3 per margin keeps the smooth tractable at n=53 (~27 basis pre-constraints).
CROSS_GAM_FORMULA <- dTmax_mean ~ s(LAI, k = 5) + s(Hmax, k = 5) +
                                    s(fCover, k = 5) +
                                    ti(FPC1, FPC2, FPC3, k = 3)

fit_cross_gams <- function(targets) {
  cli_h1("Cross-GAM step 2 — symmetric fits")
  cli_alert("Formula : {.code {deparse1(CROSS_GAM_FORMULA)}}")
  cli_alert("Family  : scat()   method = REML")

  gam_clhs <- mgcv::gam(CROSS_GAM_FORMULA, data = targets$cLHS,
                         family = mgcv::scat(), method = "REML")
  gam_hobo <- mgcv::gam(CROSS_GAM_FORMULA, data = targets$HOBO,
                         family = mgcv::scat(), method = "REML")
  list(cLHS = gam_clhs, HOBO = gam_hobo)
}

# ---- Step 3 : partial smooths on joint grid ----------------------------------
# For LAI/Hmax/fCover : partial effect from s() term varying along the axis.
# For FPC (3D ti) : 1D slice along FPC1 at joint median FPC2/FPC3.
extract_partial_smooths <- function(targets, gams, n_grid = 100L) {
  cli_h1("Cross-GAM step 3 — partial smooths on joint grid")
  axis_vars <- c("LAI", "Hmax", "fCover", "FPC1")    # plotting axes
  all_pred_vars <- c("LAI", "Hmax", "fCover", "FPC1", "FPC2", "FPC3")

  pooled <- rbind(targets$cLHS[, ..all_pred_vars],
                   targets$HOBO[, ..all_pred_vars])
  joint_median <- pooled[, lapply(.SD, function(x) median(x, na.rm = TRUE))]
  cli_alert("Joint medians : LAI={round(joint_median$LAI,3)}  Hmax={round(joint_median$Hmax,2)}  fCover={round(joint_median$fCover,3)}  FPC1={round(joint_median$FPC1,4)}  FPC2={round(joint_median$FPC2,4)}  FPC3={round(joint_median$FPC3,4)}")

  preds_all <- vector("list", length(axis_vars) * 2L)
  k <- 0L
  ti_term <- "ti(FPC1,FPC2,FPC3)"

  for (v in axis_vars) {
    rng <- range(pooled[[v]], na.rm = TRUE)
    grid_v <- seq(rng[1], rng[2], length.out = n_grid)
    fixed <- setdiff(all_pred_vars, v)
    for (w in c("cLHS", "HOBO")) {
      newd <- as.data.frame(c(setNames(list(grid_v), v),
                                lapply(fixed, function(u)
                                  rep(joint_median[[u]], n_grid))))
      names(newd) <- c(v, fixed)
      pr <- predict(gams[[w]], newdata = newd,
                     type = "terms", se.fit = TRUE)
      term_name <- if (v == "FPC1") ti_term else paste0("s(", v, ")")
      fit_v <- as.numeric(pr$fit[, term_name])
      se_v  <- as.numeric(pr$se.fit[, term_name])
      k <- k + 1L
      preds_all[[k]] <- data.table(
        variable = if (v == "FPC1") "FPC (slice@med)" else v,
        world    = w,
        x        = grid_v,
        fit      = fit_v,
        se       = se_v,
        ci_lo    = fit_v - 1.96 * se_v,
        ci_hi    = fit_v + 1.96 * se_v
      )
    }
  }
  preds <- rbindlist(preds_all)
  var_levels <- c("LAI", "Hmax", "fCover", "FPC (slice@med)")
  preds[, variable := factor(variable, levels = var_levels)]
  preds[, world    := factor(world, levels = c("cLHS", "HOBO"))]
  preds[]
}

# ---- Step 4 : superposition figure -------------------------------------------
plot_cross_smooths <- function(targets, gams, preds, out_path) {
  library(ggplot2)
  cli_h1("Cross-GAM step 4 — superposition figure")

  # annotation per panel : edf + p-value per world
  var_levels <- c("LAI", "Hmax", "fCover", "FPC (slice@med)")
  term_lookup <- c("LAI" = "s(LAI)", "Hmax" = "s(Hmax)",
                    "fCover" = "s(fCover)",
                    "FPC (slice@med)" = "ti(FPC1,FPC2,FPC3)")
  ann_rows <- vector("list", length(var_levels) * 2L); k <- 0L
  for (v in var_levels) {
    for (w in c("cLHS", "HOBO")) {
      sm  <- summary(gams[[w]])
      idx <- match(term_lookup[v], rownames(sm$s.table))
      ed  <- sm$s.table[idx, "edf"]
      pv  <- sm$s.table[idx, "p-value"]
      pv_lab <- if (pv < 1e-3) sprintf("p<0.001") else sprintf("p=%.3f", pv)
      k <- k + 1L
      ann_rows[[k]] <- data.table(variable = v, world = w,
                                    label = sprintf("%s : edf=%.1f, %s", w, ed, pv_lab))
    }
  }
  ann <- rbindlist(ann_rows)
  ann[, variable := factor(variable, levels = var_levels)]
  ann[, world    := factor(world, levels = c("cLHS", "HOBO"))]

  # rug data : raw observation points (FPC slice uses FPC1)
  rug_dt <- rbindlist(list(
    targets$cLHS[, .(x = LAI,    variable = "LAI",              world = "cLHS")],
    targets$cLHS[, .(x = Hmax,   variable = "Hmax",             world = "cLHS")],
    targets$cLHS[, .(x = fCover, variable = "fCover",           world = "cLHS")],
    targets$cLHS[, .(x = FPC1,   variable = "FPC (slice@med)",  world = "cLHS")],
    targets$HOBO[, .(x = LAI,    variable = "LAI",              world = "HOBO")],
    targets$HOBO[, .(x = Hmax,   variable = "Hmax",             world = "HOBO")],
    targets$HOBO[, .(x = fCover, variable = "fCover",           world = "HOBO")],
    targets$HOBO[, .(x = FPC1,   variable = "FPC (slice@med)",  world = "HOBO")]
  ))
  rug_dt[, variable := factor(variable, levels = var_levels)]
  rug_dt[, world    := factor(world, levels = c("cLHS", "HOBO"))]

  cols  <- c(cLHS = "#2C5F2D", HOBO = "#FFB400")
  fills <- c(cLHS = "#2C5F2D", HOBO = "#FFB400")

  p <- ggplot(preds, aes(x = x, y = fit, colour = world, fill = world)) +
    geom_hline(yintercept = 0, linetype = "dotted", colour = "grey50") +
    geom_ribbon(aes(ymin = ci_lo, ymax = ci_hi), alpha = 0.25, colour = NA) +
    geom_line(linewidth = 1.1) +
    geom_rug(data = rug_dt, aes(x = x, colour = world), inherit.aes = FALSE,
              sides = "b", alpha = 0.45, length = unit(0.02, "npc")) +
    geom_text(data = ann[world == "cLHS"],
               aes(x = -Inf, y = Inf, label = label),
               hjust = -0.05, vjust = 1.4, size = 3.0, colour = cols["cLHS"],
               inherit.aes = FALSE, fontface = "bold") +
    geom_text(data = ann[world == "HOBO"],
               aes(x = -Inf, y = Inf, label = label),
               hjust = -0.05, vjust = 2.9, size = 3.0, colour = cols["HOBO"],
               inherit.aes = FALSE, fontface = "bold") +
    facet_wrap(~ variable, scales = "free_x", nrow = 1L) +
    scale_colour_manual(values = cols, name = NULL,
                         labels = c("cLHS sim (n=400)", "HOBO obs (n=53)")) +
    scale_fill_manual(values = fills, guide = "none") +
    labs(
      title    = "GAM partial response of dTmax_mean to LiDAR traits — sim vs obs",
      subtitle = "Formula  dTmax_mean ~ s(LAI) + s(Hmax) + s(fCover) + ti(FPC1,FPC2,FPC3, k=3), family=scat(), REML  |  others fixed at JOINT median",
      x = "Trait value (joint cLHS + HOBO range)",
      y = "Partial effect on dTmax_mean (K)",
      caption = paste0(
        "FPC panel : 1D slice along FPC1 at joint median of FPC2 and FPC3. ti() captures 3D LAD shape interaction.\n",
        "CI95 = fit ± 1.96 × SE from predict(type='terms', se.fit=TRUE). Rug : raw observations per world."
      )
    ) +
    theme_bw(base_size = 11) +
    theme(
      plot.title    = element_text(face = "bold"),
      plot.subtitle = element_text(colour = "grey25", size = 9),
      plot.caption  = element_text(size = 7.5, colour = "grey35", hjust = 0,
                                    lineheight = 1.2),
      strip.background = element_rect(fill = "grey92", colour = NA),
      strip.text       = element_text(face = "bold", size = 11),
      panel.grid.minor = element_blank(),
      legend.position  = "bottom"
    )

  dir.create(dirname(out_path), recursive = TRUE, showWarnings = FALSE)
  ggsave(out_path, p, width = 11.4, height = 4.6, dpi = 300)
  cli_alert_success("Saved {.path {out_path}}")
  invisible(out_path)
}

# ---- Step 5 : variance partial explained -------------------------------------
# Variables : LAI, Hmax, fCover dropped one at a time via s(); FPC block
# dropped via removing the whole ti(FPC1, FPC2, FPC3) term.
compute_partial_varexp <- function(targets, gams) {
  cli_h1("Cross-GAM step 5 — partial deviance explained")
  vars <- c("LAI", "Hmax", "fCover", "FPC")
  rows <- vector("list", 2L * length(vars)); k <- 0L

  build_reduced_fml <- function(drop_var) {
    parts <- c(if (drop_var != "LAI")    "s(LAI, k=5)",
                if (drop_var != "Hmax")   "s(Hmax, k=5)",
                if (drop_var != "fCover") "s(fCover, k=5)",
                if (drop_var != "FPC")    "ti(FPC1, FPC2, FPC3, k=3)")
    as.formula(paste("dTmax_mean ~", paste(parts, collapse = " + ")))
  }

  # Use response-scale R² for partial drop comparison : invariant to scat()/REML
  # quirks across refits with tensor smooths.
  r2_resp <- function(g, y) {
    yhat <- as.numeric(predict(g, type = "response"))
    1 - sum((y - yhat)^2, na.rm = TRUE) / sum((y - mean(y, na.rm = TRUE))^2, na.rm = TRUE)
  }

  for (w in c("cLHS", "HOBO")) {
    g_full   <- gams[[w]]
    y        <- targets[[w]]$dTmax_mean
    R2_full  <- r2_resp(g_full, y)
    for (v in vars) {
      fml_red <- build_reduced_fml(v)
      g_red <- mgcv::gam(fml_red, data = targets[[w]],
                          family = mgcv::scat(), method = "REML")
      R2_red <- r2_resp(g_red, y)
      k <- k + 1L
      rows[[k]] <- data.table(world = w, variable = v,
                                R2_full = R2_full, R2_red = R2_red,
                                VarExp = max(R2_full - R2_red, 0))
    }
  }
  ve <- rbindlist(rows)
  ve[, VarExp_norm := VarExp / sum(VarExp), by = world]
  ve[, Rank := frank(-VarExp, ties.method = "min"), by = world]
  ve[, variable := factor(variable, levels = vars)]
  ve[, world    := factor(world, levels = c("cLHS", "HOBO"))]
  ve[]
}

plot_cross_varexp <- function(ve, out_path) {
  library(ggplot2)
  cli_h1("Cross-GAM step 5b — varexp barplot")
  cols <- c(cLHS = "#2C5F2D", HOBO = "#FFB400")
  p <- ggplot(ve, aes(x = variable, y = VarExp_norm, fill = world)) +
    geom_col(position = position_dodge(width = 0.7), width = 0.65,
              colour = "white") +
    geom_text(aes(label = sprintf("#%d\n%.1f%%", Rank, 100 * VarExp_norm)),
               position = position_dodge(width = 0.7),
               vjust = -0.2, size = 3.1, fontface = "bold",
               colour = "grey20") +
    scale_fill_manual(values = cols, name = NULL,
                       labels = c("cLHS sim", "HOBO obs")) +
    scale_y_continuous(labels = scales::percent_format(accuracy = 1),
                        expand = expansion(mult = c(0, 0.18))) +
    labs(
      title    = "Partial R² explained per trait — sim vs obs",
      subtitle = "Normalised per world (bars sum to 100%). Rank annotated.",
      x = NULL, y = "Share of explained R²",
      caption = paste0(
        "VarExp(v) = R²(full) - R²(no v) on response scale.  FPC = ti(FPC1, FPC2, FPC3, k=3) 3D LAD shape.\n",
        "Hmax in HOBO clipped to 0 : refit without Hmax gave marginally HIGHER R² (Hmax NS, p=0.43)."
      )
    ) +
    theme_bw(base_size = 11) +
    theme(plot.title = element_text(face = "bold"),
           plot.subtitle = element_text(colour = "grey25"),
           plot.caption = element_text(size = 8, colour = "grey35", hjust = 0),
           legend.position = "bottom")
  dir.create(dirname(out_path), recursive = TRUE, showWarnings = FALSE)
  ggsave(out_path, p, width = 7.5, height = 4.6, dpi = 300)
  cli_alert_success("Saved {.path {out_path}}")
}

# ---- Step 6 : recap table + classification -----------------------------------
build_recap_table <- function(gams, preds, ve, out_path) {
  cli_h1("Cross-GAM step 6 — recap table")
  var_axes  <- c("LAI", "Hmax", "fCover", "FPC (slice@med)")
  ve_keys   <- c("LAI", "Hmax", "fCover", "FPC")
  term_keys <- c("s(LAI)", "s(Hmax)", "s(fCover)", "ti(FPC1,FPC2,FPC3)")

  rows <- vector("list", length(var_axes)); k <- 0L
  for (i in seq_along(var_axes)) {
    v <- var_axes[i]; vk <- ve_keys[i]; tk <- term_keys[i]
    pc <- preds[variable == v & world == "cLHS"]
    ph <- preds[variable == v & world == "HOBO"]
    setorder(pc, x); setorder(ph, x)
    tau <- cor(pc$fit, ph$fit, method = "kendall")
    overlap <- mean( (pc$ci_lo <= ph$ci_hi) & (ph$ci_lo <= pc$ci_hi) ) * 100

    sm_c <- summary(gams$cLHS)$s.table
    sm_h <- summary(gams$HOBO)$s.table
    edf_c <- sm_c[tk, "edf"]; edf_h <- sm_h[tk, "edf"]
    p_c   <- sm_c[tk, "p-value"]; p_h <- sm_h[tk, "p-value"]

    ve_c <- ve[world == "cLHS" & variable == vk]
    ve_h <- ve[world == "HOBO" & variable == vk]

    convergence <- if (tau > 0.8 & overlap > 70) "CONVERGE"
      else if (tau < 0.5 | overlap < 40) "DIVERGE"
      else "PARTIAL"

    k <- k + 1L
    rows[[k]] <- data.table(
      Variable        = vk,
      edf_cLHS        = round(edf_c, 2),
      edf_HOBO        = round(edf_h, 2),
      p_cLHS          = signif(p_c, 3),
      p_HOBO          = signif(p_h, 3),
      VarExp_cLHS_pct = round(100 * ve_c$VarExp_norm, 1),
      VarExp_HOBO_pct = round(100 * ve_h$VarExp_norm, 1),
      Rank_cLHS       = ve_c$Rank,
      Rank_HOBO       = ve_h$Rank,
      Kendall_tau     = round(tau, 3),
      CI_overlap_pct  = round(overlap, 1),
      Convergence     = convergence
    )
  }
  tab <- rbindlist(rows)

  dir.create(dirname(out_path), recursive = TRUE, showWarnings = FALSE)
  fwrite(tab, out_path)
  cli_alert_success("Saved {.path {out_path}}")
  print(tab)

  # Console synthesis
  cli_h2("Classification convergence sim vs obs")
  for (i in seq_len(nrow(tab))) {
    r <- tab[i]
    tag <- switch(r$Convergence,
                    "CONVERGE" = "[CONVERGE]",
                    "PARTIAL"  = "[PARTIAL ]",
                    "DIVERGE"  = "[DIVERGE ]")
    cli_alert(sprintf("  %s  %-7s  tau=%+.2f  overlap=%5.1f%%  |  Rank cLHS=%d  HOBO=%d",
                       tag, r$Variable, r$Kendall_tau, r$CI_overlap_pct,
                       r$Rank_cLHS, r$Rank_HOBO))
  }
  invisible(tab)
}

print_gam_diagnostics <- function(gams) {
  for (w in c("cLHS", "HOBO")) {
    g <- gams[[w]]
    cli_h2(sprintf("Summary GAM_%s   (n = %d)", w, length(g$residuals)))
    print(summary(g))
    cli_h2(sprintf("Concurvity matrix (full) — GAM_%s", w))
    print(round(mgcv::concurvity(g, full = TRUE), 3))
  }
}

# ---- script entry : run step depending on STEP argument ----------------------
if (sys.nframe() == 0) {
  args <- commandArgs(trailingOnly = TRUE)
  step <- if (length(args) >= 1L) args[[1L]] else "1"

  if (step == "1") {
  out <- load_cross_targets()

  cli_h2("HEAD DT_clhs")
  print(head(out$cLHS, 4))
  cli_h2("SUMMARY dTmax_mean (cLHS)")
  print(summary(out$cLHS$dTmax_mean))

  cli_h2("HEAD DT_hobo")
  print(head(out$HOBO, 4))
  cli_h2("SUMMARY dTmax_mean (HOBO observed)")
  print(summary(out$HOBO$dTmax_mean))

  cli_h2("Coverage range per trait (joint)")
  for (v in c("LAI", "Hmax", "fCover", "FPC1")) {
    r_c <- range(out$cLHS[[v]], na.rm = TRUE)
    r_h <- range(out$HOBO[[v]], na.rm = TRUE)
    cli_alert(sprintf("%-7s  cLHS [%.3f ; %.3f]  HOBO [%.3f ; %.3f]",
                       v, r_c[1], r_c[2], r_h[1], r_h[2]))
  }

  # Cache step 1 output for downstream steps
  out_path <- here::here("outputs/lovb/data/DT_cross_gam_targets.rds")
  saveRDS(out, out_path)
  cli_alert_success("Saved {.path {out_path}}")

  } else if (step == "2") {
    targets <- readRDS(here::here("outputs/lovb/data/DT_cross_gam_targets.rds"))
    gams <- fit_cross_gams(targets)
    print_gam_diagnostics(gams)
    saveRDS(gams, here::here("outputs/lovb/data/cross_gam_fits.rds"))
    cli_alert_success("Saved cross_gam_fits.rds")
  } else if (step == "6") {
    gams  <- readRDS(here::here("outputs/lovb/data/cross_gam_fits.rds"))
    preds <- readRDS(here::here("outputs/lovb/data/cross_gam_predictions.rds"))
    ve    <- readRDS(here::here("outputs/lovb/data/cross_gam_varexp.rds"))
    build_recap_table(gams, preds, ve,
                       here::here("outputs/lovb/tables/tab_GAM_cross_validation.csv"))
  } else if (step == "4") {
    targets <- readRDS(here::here("outputs/lovb/data/DT_cross_gam_targets.rds"))
    gams    <- readRDS(here::here("outputs/lovb/data/cross_gam_fits.rds"))
    preds   <- readRDS(here::here("outputs/lovb/data/cross_gam_predictions.rds"))
    plot_cross_smooths(targets, gams, preds,
                        here::here("outputs/lovb/figures/fig_GAM_cross_smooths.png"))
    ve <- compute_partial_varexp(targets, gams)
    cli_h2("Variance explained table")
    print(ve)
    saveRDS(ve, here::here("outputs/lovb/data/cross_gam_varexp.rds"))
    plot_cross_varexp(ve,
                       here::here("outputs/lovb/figures/fig_GAM_cross_varexp.png"))
  } else if (step == "3") {
    targets <- readRDS(here::here("outputs/lovb/data/DT_cross_gam_targets.rds"))
    gams    <- readRDS(here::here("outputs/lovb/data/cross_gam_fits.rds"))
    preds   <- extract_partial_smooths(targets, gams)
    cli_h2("HEAD predictions table (long)")
    print(head(preds, 6))
    cli_h2("Predictions per variable x world")
    print(preds[, .N, by = .(variable, world)])
    saveRDS(preds, here::here("outputs/lovb/data/cross_gam_predictions.rds"))
    cli_alert_success("Saved cross_gam_predictions.rds")
  } else {
    stop("Unknown step : ", step)
  }
}
