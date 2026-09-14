# ==============================================================================
# MEB 2026 v9 — refresh after meeting feedback (2026-05-22 AM)
#
# Changes vs v8 :
#   - rename "Archetype" → "Profile" in heatmap and elsewhere
#   - baselines figure : 4 panels showing cLHS distributions with baselines
#   - LOO-only narrative figure : simple horizontal schema
# ==============================================================================

suppressMessages({
  library(data.table); library(here); library(cli); library(ggplot2)
  library(patchwork); library(cowplot)
})

DATA <- here::here("outputs/lovb_floor05/data")
OUT  <- here::here("outputs/figs_MEB2026_final")
dir.create(OUT, showWarnings = FALSE, recursive = TRUE)

VARS <- c("LAI","Hmax","fCover","LAD")
.PAL_PROFILE <- c("P1"="#E69F00","P2"="#0072B2","P3"="#009E73","P4"="#CC79A7")

# ============================================================================
# (1) Heatmap LOO × 4 profiles (gradient) — "Profile" labeling
# ============================================================================
build_heatmap_profiles <- function() {
  cli_h1("Heatmap LOO × 4 profiles (gradient, 'Profile' labeling)")
  arch <- readRDS(file.path(DATA, "DT_contrib_archetypes_floor05.rds"))
  long_rows <- list()
  for (v in VARS) {
    col_lovb <- paste0("Tmax_mean_LOVB_", v)
    if (col_lovb %in% names(arch) && "Tmax_mean_REF" %in% names(arch)) {
      # Convention : Δ_v = T_max(REF) − T_max(LOO_v) ; negative ⇒ v buffers
      long_rows[[length(long_rows) + 1L]] <- data.table(
        Cluster = arch$Cluster, variable = v,
        delta = arch$Tmax_mean_REF - arch[[col_lovb]])
    }
  }
  long <- rbindlist(long_rows)
  # Cluster = RAW k-means code; map via cluster_relabel (raw 1->P4, 2->P2, 3->P1, 4->P3).
  # Direct P%d was WRONG (scrambled archetype labels).
  .RELAB <- c("1" = "P4", "2" = "P2", "3" = "P1", "4" = "P3")
  long[, Profile := unname(.RELAB[as.character(Cluster)])]

  # Add "Avg" column = mean of |Δ_v| across the 4 profiles, signed by majority
  avg_rows <- long[, .(Profile = "Avg",
                          delta = mean(delta, na.rm = TRUE)),
                      by = .(variable)]
  long <- rbind(long[, .(variable, Profile, delta)], avg_rows)
  long[, Profile := factor(Profile, levels = c(paste0("P", 1:4), "Avg"))]
  long[, variable := factor(variable, levels = VARS)]
  long[, var_lab := ifelse(variable == "Hmax", "H[max]", as.character(variable))]
  long[, var_lab := factor(var_lab, levels = c("LAI","H[max]","fCover","LAD"))]

  max_abs <- max(abs(long$delta), na.rm = TRUE) * 1.05
  # Separator between profile columns and average column
  sep_x <- 4.5

  p <- ggplot(long, aes(x = Profile, y = var_lab, fill = delta)) +
    geom_tile(colour = "white", linewidth = 1) +
    geom_text(aes(label = sprintf("%+.2f", delta)),
               fontface = "bold", size = 5,
               colour = ifelse(abs(long$delta) > max_abs * 0.55,
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
  fig_path <- file.path(OUT, "fig_heatmap_LOO_profiles_gradient.png")
  ggsave(fig_path, p, width = 12, height = 5, dpi = 300)
  cli_alert_success("Saved {.path {fig_path}}")
}

# ============================================================================
# (2) Profile structural profiles — relabel C1-C4 → P1-P4
# ============================================================================
build_profiles_structural <- function() {
  cli_h1("Profiles structural figure (LAD + traits, P1-P4)")
  df_sample <- as.data.table(readRDS("out_files/Sensitivity_Analysis/clhs_sample.rds"))
  if ("Archetype" %in% names(df_sample) && !"Cluster" %in% names(df_sample))
    df_sample[, Cluster := Archetype]
  df_sample[, fCover_floored := pmax(fCover, 0.5)]
  lad_cols <- grep("^LAD_Layer_", names(df_sample), value = TRUE)
  z_breaks <- as.numeric(gsub("LAD_Layer_", "", lad_cols))

  profiles <- list()
  for (cl in c("1","2","3","4")) {
    sub <- df_sample[as.character(Cluster) == cl, ..lad_cols]
    sub_mat <- as.matrix(sub); sub_mat[is.na(sub_mat)] <- 0
    profiles[[cl]] <- data.table(
      Profile = paste0("P", cl),
      height  = z_breaks,
      lad_mean = colMeans(sub_mat),
      lad_sd   = apply(sub_mat, 2, sd))
  }
  prof_long <- rbindlist(profiles)
  prof_long[, Profile := factor(Profile, levels = paste0("P", 1:4))]

  stats <- df_sample[, .(LAI_mean   = mean(LAI, na.rm = TRUE),
                            LAI_min   = min(LAI, na.rm = TRUE),
                            LAI_max   = max(LAI, na.rm = TRUE),
                            Hmax_mean = mean(Hmax, na.rm = TRUE),
                            Hmax_min  = min(Hmax, na.rm = TRUE),
                            Hmax_max  = max(Hmax, na.rm = TRUE),
                            fCover_mean = mean(fCover_floored, na.rm = TRUE),
                            fCover_min  = min(fCover_floored, na.rm = TRUE),
                            fCover_max  = max(fCover_floored, na.rm = TRUE)),
                        by = Cluster]
  # raw k-means code -> display label (raw 1->P4, 2->P2, 3->P1, 4->P3); direct P%d was WRONG
  .RELAB <- c("1" = "P4", "2" = "P2", "3" = "P1", "4" = "P3")
  stats[, Profile := unname(.RELAB[as.character(Cluster)])]
  stats[, Profile := factor(Profile, levels = paste0("P", 1:4))]
  setorder(stats, Profile)
  # explicit MuSICA archetype value (cluster mean, Hmax rounded as in sims)
  # + sampling range [min, max] across the 100 plots of the cluster
  stats[, Hmax_used := round(Hmax_mean)]
  stats[, lab_LAI    := sprintf("LAI*':  '*'%.2f'", LAI_mean)]
  stats[, lab_LAI_r  := sprintf("'[%.2f - %.2f]'", LAI_min, LAI_max)]
  stats[, lab_Hmax   := sprintf("H[max]*':  '*%d*' m'", as.integer(Hmax_used))]
  stats[, lab_Hmax_r := sprintf("'[%.0f - %.0f] m'", round(Hmax_min), round(Hmax_max))]
  stats[, lab_fCover := sprintf("fCover*':  '*%.2f", fCover_mean)]
  stats[, lab_fCover_r := sprintf("'[%.2f - %.2f]'", fCover_min, fCover_max)]

  pal <- setNames(unname(.PAL_PROFILE), paste0("P", 1:4))
  max_lad <- max(prof_long$lad_mean + prof_long$lad_sd, na.rm = TRUE) * 1.05
  x_ann <- max_lad * 0.18
  p <- ggplot(prof_long, aes(x = lad_mean, y = height,
                                colour = Profile, fill = Profile)) +
    geom_ribbon(aes(xmin = pmax(0, lad_mean - lad_sd),
                      xmax = lad_mean + lad_sd),
                  alpha = 0.25, colour = NA) +
    geom_path(linewidth = 1.2) +
    geom_hline(data = stats, aes(yintercept = Hmax_used, colour = Profile),
                  linetype = "dashed", linewidth = 0.8,
                  inherit.aes = FALSE, show.legend = FALSE) +
    # Header value labels (centroid value used in MuSICA)
    geom_text(data = stats, aes(x = x_ann, y = 41, label = lab_LAI,
                                    colour = Profile),
               parse = TRUE, inherit.aes = FALSE, hjust = 0, vjust = 1,
               size = 5.2, fontface = "bold", show.legend = FALSE) +
    geom_text(data = stats, aes(x = x_ann, y = 38, label = lab_Hmax,
                                    colour = Profile),
               parse = TRUE, inherit.aes = FALSE, hjust = 0, vjust = 1,
               size = 5.2, fontface = "bold", show.legend = FALSE) +
    geom_text(data = stats, aes(x = x_ann, y = 35, label = lab_fCover,
                                    colour = Profile),
               parse = TRUE, inherit.aes = FALSE, hjust = 0, vjust = 1,
               size = 5.2, fontface = "bold", show.legend = FALSE) +
    # Range subtitle in lighter weight just below
    geom_text(data = stats, aes(x = x_ann, y = 39.2, label = lab_LAI_r,
                                    colour = Profile),
               parse = TRUE, inherit.aes = FALSE, hjust = 0, vjust = 1,
               size = 3.8, show.legend = FALSE) +
    geom_text(data = stats, aes(x = x_ann, y = 36.2, label = lab_Hmax_r,
                                    colour = Profile),
               parse = TRUE, inherit.aes = FALSE, hjust = 0, vjust = 1,
               size = 3.8, show.legend = FALSE) +
    geom_text(data = stats, aes(x = x_ann, y = 33.2, label = lab_fCover_r,
                                    colour = Profile),
               parse = TRUE, inherit.aes = FALSE, hjust = 0, vjust = 1,
               size = 3.8, show.legend = FALSE) +
    facet_wrap(~ Profile, nrow = 1L) +
    scale_colour_manual(values = pal, name = NULL, drop = FALSE) +
    scale_fill_manual(values = pal, name = NULL, drop = FALSE) +
    scale_x_continuous(limits = c(0, max_lad * 0.95),
                        expand = expansion(mult = c(0, 0.02))) +
    scale_y_continuous(limits = c(0, 45), breaks = c(0, 10, 20, 30, 40)) +
    labs(x = bquote("LAD" ~ "(m"^"2"~"m"^"-3"*")"),
          y = "Height (m)") +
    theme_bw(base_size = 15) +
    theme(strip.background = element_rect(fill = "grey92", colour = NA),
           strip.text = element_text(face = "bold", size = 14),
           legend.position = "none",
           panel.grid.minor = element_blank())
  fig_path <- file.path(OUT, "fig_profiles_structural.png")
  ggsave(fig_path, p, width = 14, height = 6, dpi = 300)
  cli_alert_success("Saved {.path {fig_path}}")
}

# ============================================================================
# (3) Baselines figure : 4 panels showing cLHS distribution + baseline lines
# ============================================================================
build_baselines_figure <- function() {
  cli_h1("Baselines distribution figure")
  clhs <- as.data.table(readRDS("out_files/Sensitivity_Analysis/clhs_sample.rds"))
  clhs[, fCover_floored := pmax(fCover, 0.5)]

  base_pal <- "#FFB400"   # baseline line colour
  hist_fill <- "#56B4E9"

  panel_hist <- function(values, baseline, title_expr, xlab,
                            extra_v = NULL, extra_label = NULL) {
    df <- data.frame(val = values)
    p <- ggplot(df, aes(x = val)) +
      geom_histogram(fill = hist_fill, colour = "white",
                       bins = 30, alpha = 0.85) +
      geom_vline(xintercept = baseline,
                   colour = base_pal, linewidth = 1.3, linetype = "dashed") +
      annotate("text", x = baseline, y = Inf,
                label = sprintf("baseline = %.2f", baseline),
                colour = base_pal, fontface = "bold",
                hjust = -0.05, vjust = 1.5, size = 4.5) +
      labs(x = xlab, y = "n plots", title = title_expr) +
      theme_bw(base_size = 13) +
      theme(plot.title = element_text(face = "bold", hjust = 0.5, size = 14),
             panel.grid.minor = element_blank())
    if (!is.null(extra_v)) {
      p <- p + geom_vline(xintercept = extra_v, colour = "grey30",
                              linewidth = 0.8, linetype = "dotted") +
        annotate("text", x = extra_v, y = Inf,
                  label = extra_label, colour = "grey30",
                  hjust = -0.05, vjust = 3, size = 3.6)
    }
    p
  }

  p_lai <- panel_hist(clhs$LAI, baseline = mean(clhs$LAI),
                          title_expr = expression(bold(LAI)), xlab = bquote(m^2 ~ m^{-2}))
  p_hmax <- panel_hist(clhs$Hmax, baseline = mean(clhs$Hmax),
                          title_expr = expression(bold(H[max])), xlab = "m")
  p_fc <- panel_hist(clhs$fCover, baseline = 1,
                        title_expr = expression(bold(fCover)),
                        xlab = "–")

  # LAD baseline = uniform shape ; envelope shown as min/max across all plots
  lad_cols <- grep("^LAD_Layer_", names(clhs), value = TRUE)
  z_breaks <- as.numeric(gsub("LAD_Layer_", "", lad_cols))
  lad_mat <- as.matrix(clhs[, ..lad_cols]); lad_mat[is.na(lad_mat)] <- 0
  lad_long <- data.table(
    height = rep(z_breaks, each = nrow(clhs)),
    plot_id = rep(seq_len(nrow(clhs)), times = length(z_breaks)),
    lad = as.vector(lad_mat))
  lad_summary <- lad_long[, .(lad_mean = mean(lad), lad_p10 = quantile(lad, 0.1),
                                lad_p90 = quantile(lad, 0.9)), by = height]
  # Baseline uniform : flat profile equal to global mean LAI / canopy_depth
  baseline_uniform <- sum(lad_summary$lad_mean) / nrow(lad_summary)

  p_lad <- ggplot(lad_summary, aes(x = lad_mean, y = height)) +
    geom_ribbon(aes(xmin = lad_p10, xmax = lad_p90),
                  fill = hist_fill, alpha = 0.4) +
    geom_path(linewidth = 1.1, colour = hist_fill) +
    geom_vline(xintercept = baseline_uniform,
                 colour = base_pal, linewidth = 1.3, linetype = "dashed") +
    annotate("text", x = baseline_uniform, y = max(z_breaks),
              label = sprintf("baseline = %.2f", baseline_uniform),
              colour = base_pal, fontface = "bold",
              hjust = -0.05, vjust = -0.5, size = 4.3) +
    labs(x = bquote("LAD" ~ "(m"^2~"m"^{-3}*")"), y = "Height (m)",
          title = "LAD vertical profile") +
    theme_bw(base_size = 13) +
    theme(plot.title = element_text(face = "bold", hjust = 0.5, size = 14),
           panel.grid.minor = element_blank())

  combined <- (p_lai | p_hmax | p_fc | p_lad) + plot_layout(nrow = 1)
  fig_path <- file.path(OUT, "fig_baselines_distributions.png")
  ggsave(fig_path, combined, width = 18, height = 5, dpi = 300)
  cli_alert_success("Saved {.path {fig_path}}")
}

# ============================================================================
# (4) LOO illustration schema (horizontal)
# ============================================================================
build_loo_schema <- function() {
  cli_h1("LOO horizontal illustration")
  # Conceptual schematic : REF (4 traits real) -> remove v -> measure ΔTmax
  steps <- data.table(
    step = factor(c("REF", "LOVB(LAI)", "LOVB(H_max)",
                       "LOVB(fCover)", "LOVB(LAD)"),
                    levels = c("REF","LOVB(LAI)","LOVB(H_max)",
                                 "LOVB(fCover)","LOVB(LAD)")),
    LAI    = c(1,0,1,1,1),
    Hmax   = c(1,1,0,1,1),
    fCover = c(1,1,1,0,1),
    LAD    = c(1,1,1,1,0))
  long <- melt(steps, id.vars = "step", variable.name = "trait", value.name = "state")
  long[, trait := factor(trait, levels = c("LAI","Hmax","fCover","LAD"))]
  long[, label := ifelse(state == 1, "real", "baseline")]
  long[, fillc := ifelse(state == 1, "#2C5F2D", "#D55E00")]
  long[, trait_lab := ifelse(trait == "Hmax", "H[max]", as.character(trait))]
  long[, trait_lab := factor(trait_lab, levels = c("LAI","H[max]","fCover","LAD"))]

  p <- ggplot(long, aes(x = trait_lab, y = step, fill = label)) +
    geom_tile(colour = "white", linewidth = 2) +
    geom_text(aes(label = label), fontface = "bold", size = 4.5,
                colour = "white") +
    scale_fill_manual(values = c("real" = "#2C5F2D", "baseline" = "#D55E00"),
                         name = NULL) +
    scale_x_discrete(labels = function(x) parse(text = as.character(x))) +
    scale_y_discrete(limits = rev) +
    labs(x = NULL, y = NULL,
          title = "LOO design : remove one trait at a time, compare to REF",
          subtitle = expression(Delta[v] == T[max]^{LOO(v)} - T[max]^{REF})) +
    theme_minimal(base_size = 15) +
    theme(plot.title    = element_text(face = "bold", size = 15),
           plot.subtitle = element_text(size = 14, hjust = 0),
           axis.text.x   = element_text(face = "bold", size = 14),
           axis.text.y   = element_text(face = "bold", size = 13),
           panel.grid    = element_blank(),
           legend.position = "bottom",
           legend.text = element_text(size = 12))
  fig_path <- file.path(OUT, "fig_LOO_schema.png")
  ggsave(fig_path, p, width = 11, height = 6, dpi = 300)
  cli_alert_success("Saved {.path {fig_path}}")
}

if (sys.nframe() == 0) {
  build_heatmap_profiles()
  build_profiles_structural()
  build_baselines_figure()
  build_loo_schema()
  cli_alert_success("v9 figures generated.")
}
