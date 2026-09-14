# ==============================================================================
# Recompute attribution with 3 unified metrics :
#   1. MAE = mean(|Δ_v|) for LOVB, LVA, Shapley per-plot (harmonized)
#   2. Shapley RMSE-based : v(S) = RMSE(sim_S vs sim_REF) on 400 plots, global
#   3. Shapley R²-based : v(S) = R²(sim_S vs sim_REF), global
# All using floor05 data.
# Outputs : tab_attribution_*_floor05_<metric>.csv and corresponding heatmaps.
# ==============================================================================

suppressMessages({
  library(data.table); library(here); library(cli); library(ggplot2)
})

DATA <- here::here("outputs/lovb_floor05/data")
TAB  <- here::here("outputs/lovb_floor05/tables")
FIG  <- here::here("outputs/lovb_floor05/figures")

VARS <- c("LAI","fCover","Hmax","LAD")
all_bits <- c("0000","0001","0010","0011","0100","0101","0110","0111",
                "1000","1001","1010","1011","1100","1101","1110","1111")

# ============================================================================
# Shapley exact (16-coalition, generic over scalar v(S))
# ============================================================================
shap_one <- function(v_func) {
  # v_func : function taking bit_string -> scalar value
  mb <- function(bo, n=4L) { b <- rep("0", n); b[bo] <- "1"; paste(b, collapse="") }
  out <- numeric(4L)
  for (vb in 1:4) {
    others <- setdiff(1:4, vb); total <- 0
    for (s in 0:3) {
      ss <- if (s == 0L) list(integer(0)) else combn(others, s, simplify = FALSE)
      w <- factorial(s) * factorial(4-s-1) / factorial(4)
      for (S in ss) total <- total + w * (v_func(mb(c(S, vb))) - v_func(mb(S)))
    }
    out[vb] <- total
  }
  setNames(out, c("LAI","Hmax","fCover","LAD"))
}

# ============================================================================
# METRIC 1 : MAE harmonized
# ============================================================================
build_MAE_tables <- function() {
  cli_h1("METRIC 1 — MAE harmonized")
  DT_c   <- readRDS(file.path(DATA, "DT_contrib_cLHS_floor05.rds"))
  DT_a   <- readRDS(file.path(DATA, "DT_contrib_archetypes_floor05.rds"))
  DT_phi <- readRDS(file.path(DATA, "DT_shapley_per_plot_floor05.rds"))

  # LOVB cLHS : mean(|Delta_v_mean|) over 400 plots
  lovb_clhs <- data.table(variable = VARS,
    LOVB_cLHS_value = c(mean(abs(DT_c$Delta_LAI_mean),    na.rm = TRUE),
                          mean(abs(DT_c$Delta_Hmax_mean),   na.rm = TRUE),
                          mean(abs(DT_c$Delta_fCover_mean), na.rm = TRUE),
                          mean(abs(DT_c$Delta_LAD_mean),    na.rm = TRUE)))
  lovb_clhs[, LOVB_cLHS_rank := frank(-LOVB_cLHS_value, ties.method = "min")]

  # LVA cLHS
  lva_clhs <- data.table(variable = VARS,
    LVA_cLHS_value = c(mean(abs(DT_c$LVA_LAI_mean),    na.rm = TRUE),
                         mean(abs(DT_c$LVA_Hmax_mean),   na.rm = TRUE),
                         mean(abs(DT_c$LVA_fCover_mean), na.rm = TRUE),
                         mean(abs(DT_c$LVA_LAD_mean),    na.rm = TRUE)))
  lva_clhs[, LVA_cLHS_rank := frank(-LVA_cLHS_value, ties.method = "min")]

  # Shapley cLHS : MAE = mean(|phi_v|) across 400 plots
  shap_clhs <- DT_phi[, .(Shapley_cLHS_value = mean(abs(phi), na.rm = TRUE)),
                        by = variable]
  shap_clhs[, Shapley_cLHS_rank := frank(-Shapley_cLHS_value, ties.method = "min")]

  tab_clhs <- merge(merge(shap_clhs, lovb_clhs, by = "variable"),
                      lva_clhs, by = "variable")
  tab_clhs[, variable := factor(variable, levels = VARS)]
  setorder(tab_clhs, variable)
  fwrite(tab_clhs, file.path(TAB, "tab_attribution_cLHS_mean_MAE.csv"))
  cli_alert("cLHS MAE :"); print(tab_clhs)

  # Archetypes
  lovb_arch <- data.table(variable = VARS,
    LOVB_arch_value = c(mean(abs(DT_a$Delta_LAI_mean),    na.rm = TRUE),
                          mean(abs(DT_a$Delta_Hmax_mean),   na.rm = TRUE),
                          mean(abs(DT_a$Delta_fCover_mean), na.rm = TRUE),
                          mean(abs(DT_a$Delta_LAD_mean),    na.rm = TRUE)))
  lovb_arch[, LOVB_arch_rank := frank(-LOVB_arch_value, ties.method = "min")]
  lva_arch <- data.table(variable = VARS,
    LVA_arch_value = c(mean(abs(DT_a$LVA_LAI_mean),    na.rm = TRUE),
                         mean(abs(DT_a$LVA_Hmax_mean),   na.rm = TRUE),
                         mean(abs(DT_a$LVA_fCover_mean), na.rm = TRUE),
                         mean(abs(DT_a$LVA_LAD_mean),    na.rm = TRUE)))
  lva_arch[, LVA_arch_rank := frank(-LVA_arch_value, ties.method = "min")]

  # Shapley archetypes : need 16-coal data
  DT_daily_a <- readRDS(file.path(DATA, "DT_daily_archetypes_floor05.rds"))
  # Load missing 6 from H1_archetypes if needed (C2/C3/C4)
  if (length(unique(DT_daily_a$bit_code)) < 16) {
    suppressMessages({library(ncdf4); library(tidyverse); library(musica.tools)})
    source(here::here("R/config.R")); source(here::here("R/musica.R"))
    source(here::here("R/io.R")); source(here::here("R/lovb_01_load.R"))
    MISSING <- c("0011","0101","0110","1001","1010","1100")
    df_macro <- extract_macro_daily(CFG$forcing_file, CFG$date_seq)
    extra <- list()
    for (arch in c("Arch_C2","Arch_C3","Arch_C4")) {
      for (bit in MISSING) {
        nc <- here::here("out_files/H1_archetypes", arch,
                            sprintf("musica_out_%s_ARCH_%s.nc", arch, lovb_bit_to_suffix(bit)))
        if (!file.exists(nc) || file.size(nc) < 1e5) next
        res <- tryCatch(extract_deltatmax_one(nc, df_macro, CFG$date_seq),
                         error = function(e) NULL)
        if (is.null(res) || nrow(res) == 0L) next
        dt <- as.data.table(res); dt[, archetype := arch]; dt[, bit_code := bit]
        extra[[length(extra)+1L]] <- dt[, .(archetype, date, bit_code, Delta_Tmax)]
      }
    }
    DT_daily_a <- rbind(DT_daily_a, rbindlist(extra, fill = TRUE), fill = TRUE)
  }
  agg_a <- DT_daily_a[, .(Tmax_mean = mean(Delta_Tmax, na.rm = TRUE)),
                        by = .(archetype, bit_code)]
  Tw_a <- dcast(agg_a, archetype ~ bit_code, value.var = "Tmax_mean")
  setnames(Tw_a, all_bits, paste0("T_", all_bits))
  phi_a <- list()
  for (i in seq_len(nrow(Tw_a))) {
    tvec <- as.numeric(unlist(Tw_a[i, paste0("T_", all_bits), with = FALSE]))
    names(tvec) <- all_bits
    if (any(is.na(tvec))) next
    v_func <- function(b) as.numeric(tvec[b]) - as.numeric(tvec["0000"])
    phi <- shap_one(v_func)
    phi_a[[Tw_a$archetype[i]]] <- data.table(archetype = Tw_a$archetype[i],
                                                variable = names(phi),
                                                phi = unname(phi))
  }
  DT_phi_a <- rbindlist(phi_a)
  shap_arch <- DT_phi_a[, .(Shapley_arch_value = mean(abs(phi), na.rm = TRUE)),
                          by = variable]
  shap_arch[, Shapley_arch_rank := frank(-Shapley_arch_value, ties.method = "min")]

  tab_arch <- merge(merge(shap_arch, lovb_arch, by = "variable"),
                      lva_arch, by = "variable")
  tab_arch[, variable := factor(variable, levels = VARS)]
  setorder(tab_arch, variable)
  fwrite(tab_arch, file.path(TAB, "tab_attribution_archetypes_pooled_MAE.csv"))
  cli_alert("Archetype MAE :"); print(tab_arch)

  list(clhs = tab_clhs, arch = tab_arch)
}

# ============================================================================
# METRIC 2 : Shapley RMSE-based (global)
#   v(S) = RMSE(sim_S vs sim_REF) on the cLHS plots, where the residuals
#   are taken at the plot level.
# ============================================================================
build_shapley_RMSE <- function() {
  cli_h1("METRIC 2 — Shapley RMSE-based (global, cLHS)")
  DT_full <- readRDS(file.path(DATA, "DT_contrib_cLHS_floor05_16coalitions.rds"))

  # RMSE per coalition vs REF (using Tmax_mean cols)
  Tmax_REF <- DT_full$Tmax_mean_1111
  v_loss <- numeric(length(all_bits)); names(v_loss) <- all_bits
  for (b in all_bits) {
    col <- paste0("Tmax_mean_", b)
    v_loss[b] <- sqrt(mean((DT_full[[col]] - Tmax_REF)^2, na.rm = TRUE))
  }
  # By definition v_loss["1111"] = 0 (RMSE of REF vs REF)
  cli_alert("RMSE per coalition (vs REF) :"); print(round(v_loss, 4))

  # Convert to Shapley : v(S) = v_loss_empty - v_loss(S)
  # Shapley value of v = sum over coalitions of marginal contribution
  v_empty <- v_loss["0000"]
  v_func <- function(b) v_empty - v_loss[b]
  phi <- shap_one(v_func)
  tab <- data.table(variable = names(phi), Shapley_RMSE_value = unname(phi))
  tab[, Shapley_RMSE_rank := frank(-abs(Shapley_RMSE_value), ties.method = "min")]
  tab[, variable := factor(variable, levels = VARS)]
  setorder(tab, variable)
  fwrite(tab, file.path(TAB, "tab_attribution_cLHS_Shapley_RMSE.csv"))
  cli_alert("Shapley RMSE-based :"); print(tab)
  tab
}

# ============================================================================
# METRIC 3 : Shapley R²-based (global)
#   v(S) = R²(sim_S, sim_REF) = correlation² between simulated values
# ============================================================================
build_shapley_R2 <- function() {
  cli_h1("METRIC 3 — Shapley R²-based (global, cLHS)")
  DT_full <- readRDS(file.path(DATA, "DT_contrib_cLHS_floor05_16coalitions.rds"))
  Tmax_REF <- DT_full$Tmax_mean_1111
  v_R2 <- numeric(length(all_bits)); names(v_R2) <- all_bits
  for (b in all_bits) {
    col <- paste0("Tmax_mean_", b)
    sim_S <- DT_full[[col]]
    ok <- !is.na(sim_S) & !is.na(Tmax_REF)
    if (sum(ok) < 3) { v_R2[b] <- NA; next }
    # R² = 1 - SS_res/SS_tot where SS_res = sum((sim_S - REF)²), SS_tot = sum((REF - mean(REF))²)
    ss_res <- sum((sim_S[ok] - Tmax_REF[ok])^2)
    ss_tot <- sum((Tmax_REF[ok] - mean(Tmax_REF[ok]))^2)
    v_R2[b] <- 1 - ss_res / ss_tot
  }
  # By definition v_R2["1111"] should be 1 (perfect match REF vs REF)
  cli_alert("R² per coalition (vs REF) :"); print(round(v_R2, 4))

  # Shapley : convert to "loss" form (1 - R²) since lower=better in loss convention
  v_loss_r2 <- 1 - v_R2
  v_empty <- v_loss_r2["0000"]
  v_func <- function(b) v_empty - v_loss_r2[b]
  phi <- shap_one(v_func)
  # Sum should equal v_loss_r2["0000"] - v_loss_r2["1111"] = (1 - R²_NULL) - 0
  tab <- data.table(variable = names(phi), Shapley_R2_value = unname(phi))
  tab[, Shapley_R2_rank := frank(-abs(Shapley_R2_value), ties.method = "min")]
  tab[, variable := factor(variable, levels = VARS)]
  setorder(tab, variable)
  fwrite(tab, file.path(TAB, "tab_attribution_cLHS_Shapley_R2.csv"))
  cli_alert("Shapley R²-based :"); print(tab)
  tab
}

# ============================================================================
# Main : all 3 metrics + heatmaps
# ============================================================================
if (sys.nframe() == 0) {
  m1 <- build_MAE_tables()
  m2 <- build_shapley_RMSE()
  m3 <- build_shapley_R2()

  # --- Comparative heatmap : 3 columns (MAE, RMSE Shapley, R² Shapley) for cLHS Shapley
  cli_h1("Comparative table : 3 Shapley metrics (cLHS, per-variable)")
  comp <- data.table(
    variable = VARS,
    MAE_rank  = m1$clhs[match(VARS, variable)]$Shapley_cLHS_rank,
    RMSE_rank = m2[match(VARS, variable)]$Shapley_RMSE_rank,
    R2_rank   = m3[match(VARS, variable)]$Shapley_R2_rank,
    MAE_value  = m1$clhs[match(VARS, variable)]$Shapley_cLHS_value,
    RMSE_value = m2[match(VARS, variable)]$Shapley_RMSE_value,
    R2_value   = m3[match(VARS, variable)]$Shapley_R2_value
  )
  fwrite(comp, file.path(TAB, "tab_Shapley_3metrics_comparison.csv"))
  cli_alert("Comparison :"); print(comp)

  # Heatmap of 3 metric rankings
  long <- melt(comp, id.vars = "variable",
                 measure.vars = c("MAE_rank","RMSE_rank","R2_rank"),
                 variable.name = "metric", value.name = "Rank")
  long[, metric := factor(metric,
                            levels = c("MAE_rank","RMSE_rank","R2_rank"),
                            labels = c("MAE\n(per-plot)","RMSE\n(global)","R²\n(global)"))]
  long[, variable := factor(variable, levels = VARS)]
  pal <- c("1"="#440154","2"="#3B528B","3"="#5DC863","4"="#FDE725")
  p <- ggplot(long, aes(x = metric, y = variable, fill = factor(Rank))) +
    geom_tile(colour = "white", linewidth = 1) +
    geom_text(aes(label = Rank), colour = ifelse(long$Rank <= 2, "white", "grey15"),
               fontface = "bold", size = 6) +
    scale_fill_manual(values = pal, name = "Rank", drop = FALSE) +
    scale_y_discrete(limits = rev) +
    labs(
      title    = "Shapley attribution ranking — 3 metrics (cLHS, n=400)",
      subtitle = "MAE = mean(|φ_v|) per plot     |     RMSE = global Shapley on residuals     |     R² = global Shapley on explained variance",
      x = NULL, y = NULL,
      caption = "All 3 metrics use exact 16-coalition Shapley. Differences reveal how the choice of metric affects ranking."
    ) +
    theme_bw(base_size = 14) +
    theme(plot.title = element_text(face="bold", size=14),
           plot.subtitle = element_text(colour="grey25", size=11),
           plot.caption = element_text(size = 9, colour = "grey35", hjust = 0),
           axis.text.x = element_text(face="bold", size=11),
           axis.text.y = element_text(face="bold", size=13),
           panel.grid = element_blank())
  ggsave(file.path(FIG, "fig_Shapley_3metrics_comparison.png"),
          p, width = 9, height = 5, dpi = 300)
  cli_alert_success("Saved fig_Shapley_3metrics_comparison.png")
}
