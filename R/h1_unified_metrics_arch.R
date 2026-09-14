# ==============================================================================
# Same 3-metric Shapley analysis but at the archetype level (n=4).
# Outputs : tab_attribution_archetypes_Shapley_{RMSE,R2}.csv + comparison figure.
# ==============================================================================

suppressMessages({
  library(data.table); library(here); library(cli); library(ggplot2)
  library(ncdf4); library(tidyverse); library(musica.tools)
})
source(here::here("R/config.R"))
source(here::here("R/io.R"))
source(here::here("R/musica.R"))
source(here::here("R/lovb_01_load.R"))

DATA <- here::here("outputs/lovb_floor05/data")
TAB  <- here::here("outputs/lovb_floor05/tables")
FIG  <- here::here("outputs/lovb_floor05/figures")

VARS <- c("LAI","fCover","Hmax","LAD")
all_bits <- c("0000","0001","0010","0011","0100","0101","0110","0111",
                "1000","1001","1010","1011","1100","1101","1110","1111")

shap_one <- function(v_func) {
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

# ---- Load daily archetype data with all 16 coalitions ------------------------
cli_h1("Archetype 3-metric Shapley (n=4)")
DT_daily_a <- readRDS(file.path(DATA, "DT_daily_archetypes_floor05.rds"))
# Check per-archetype : if ANY arch < 16, augment
n_per_arch <- DT_daily_a[, length(unique(bit_code)), by = archetype]
if (min(n_per_arch$V1) < 16) {
  MISSING <- c("0011","0101","0110","1001","1010","1100")
  df_macro <- extract_macro_daily(CFG$forcing_file, CFG$date_seq)
  extra <- list()
  for (arch in c("Arch_C2","Arch_C3","Arch_C4")) {
    for (bit in MISSING) {
      nc <- here::here("out_files/H1_archetypes", arch,
                          sprintf("musica_out_%s_ARCH_%s.nc", arch,
                                    lovb_bit_to_suffix(bit)))
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

# Aggregate mean per (archetype, coalition)
DT_agg <- DT_daily_a[, .(Tmax_mean = mean(Delta_Tmax, na.rm = TRUE)),
                      by = .(archetype, bit_code)]
DT_w <- dcast(DT_agg, archetype ~ bit_code, value.var = "Tmax_mean")
have_cols <- intersect(all_bits, names(DT_w))
setnames(DT_w, have_cols, paste0("T_", have_cols))
cli_alert("Wide DT : {nrow(DT_w)} archetypes × {length(have_cols)} coalitions")

# ============================================================================
# METRIC 2 : Shapley RMSE-based across 4 archetypes
# ============================================================================
cli_h2("Shapley RMSE-based (n=4 archetypes)")
T_REF <- DT_w$T_1111
v_RMSE <- numeric(length(all_bits)); names(v_RMSE) <- all_bits
for (b in all_bits) {
  col <- paste0("T_", b)
  if (!(col %in% names(DT_w))) { v_RMSE[b] <- NA; next }
  v_RMSE[b] <- sqrt(mean((DT_w[[col]] - T_REF)^2, na.rm = TRUE))
}
cli_alert("RMSE per coalition (n=4) :"); print(round(v_RMSE, 4))
v_empty <- v_RMSE["0000"]
v_func <- function(b) v_empty - v_RMSE[b]
phi_RMSE <- shap_one(v_func)
tab_RMSE <- data.table(variable = names(phi_RMSE),
                          Shapley_arch_RMSE_value = unname(phi_RMSE))
tab_RMSE[, Shapley_arch_RMSE_rank := frank(-abs(Shapley_arch_RMSE_value),
                                                ties.method = "min")]
tab_RMSE[, variable := factor(variable, levels = VARS)]
setorder(tab_RMSE, variable)
fwrite(tab_RMSE, file.path(TAB, "tab_attribution_archetypes_Shapley_RMSE.csv"))
cli_alert("Shapley RMSE archetypes :"); print(tab_RMSE)

# ============================================================================
# METRIC 3 : Shapley R²-based across 4 archetypes
# ============================================================================
cli_h2("Shapley R²-based (n=4 archetypes)")
v_R2 <- numeric(length(all_bits)); names(v_R2) <- all_bits
for (b in all_bits) {
  col <- paste0("T_", b)
  if (!(col %in% names(DT_w))) { v_R2[b] <- NA; next }
  sim_S <- DT_w[[col]]
  ok <- !is.na(sim_S) & !is.na(T_REF)
  if (sum(ok) < 3) { v_R2[b] <- NA; next }
  ss_res <- sum((sim_S[ok] - T_REF[ok])^2)
  ss_tot <- sum((T_REF[ok] - mean(T_REF[ok]))^2)
  v_R2[b] <- 1 - ss_res / ss_tot
}
cli_alert("R² per coalition (n=4) :"); print(round(v_R2, 4))
v_loss_r2 <- 1 - v_R2
v_empty_r2 <- v_loss_r2["0000"]
v_func_r2 <- function(b) v_empty_r2 - v_loss_r2[b]
phi_R2 <- shap_one(v_func_r2)
tab_R2 <- data.table(variable = names(phi_R2),
                        Shapley_arch_R2_value = unname(phi_R2))
tab_R2[, Shapley_arch_R2_rank := frank(-abs(Shapley_arch_R2_value),
                                            ties.method = "min")]
tab_R2[, variable := factor(variable, levels = VARS)]
setorder(tab_R2, variable)
fwrite(tab_R2, file.path(TAB, "tab_attribution_archetypes_Shapley_R2.csv"))
cli_alert("Shapley R² archetypes :"); print(tab_R2)

# ============================================================================
# Comparison heatmap (3 metrics × 4 vars) at archetype level
# ============================================================================
mae_arch <- fread(file.path(TAB, "tab_attribution_archetypes_pooled_MAE.csv"))
comp <- data.table(
  variable = VARS,
  MAE_rank  = mae_arch[match(VARS, variable)]$Shapley_arch_rank,
  RMSE_rank = tab_RMSE[match(VARS, variable)]$Shapley_arch_RMSE_rank,
  R2_rank   = tab_R2[match(VARS, variable)]$Shapley_arch_R2_rank,
  MAE_value  = mae_arch[match(VARS, variable)]$Shapley_arch_value,
  RMSE_value = tab_RMSE[match(VARS, variable)]$Shapley_arch_RMSE_value,
  R2_value   = tab_R2[match(VARS, variable)]$Shapley_arch_R2_value
)
fwrite(comp, file.path(TAB, "tab_Shapley_3metrics_archetypes_comparison.csv"))
cli_alert("Archetype comparison :"); print(comp)

long <- melt(comp, id.vars = "variable",
               measure.vars = c("MAE_rank","RMSE_rank","R2_rank"),
               variable.name = "metric", value.name = "Rank")
long[, metric := factor(metric,
                          levels = c("MAE_rank","RMSE_rank","R2_rank"),
                          labels = c("MAE\n(per-archetype)","RMSE\n(n=4)","R²\n(n=4)"))]
long[, variable := factor(variable, levels = VARS)]
pal <- c("1"="#440154","2"="#3B528B","3"="#5DC863","4"="#FDE725")
p <- ggplot(long, aes(x = metric, y = variable, fill = factor(Rank))) +
  geom_tile(colour = "white", linewidth = 1) +
  geom_text(aes(label = Rank), colour = ifelse(long$Rank <= 2, "white", "grey15"),
             fontface = "bold", size = 6) +
  scale_fill_manual(values = pal, name = "Rank", drop = FALSE) +
  scale_y_discrete(limits = rev) +
  labs(
    title    = "Shapley attribution ranking — 3 metrics (archetypes, n=4)",
    subtitle = "MAE = mean(|φ_v|) per archetype     |     RMSE / R² = global Shapley over 4 archetypes",
    x = NULL, y = NULL,
    caption = "With n=4 archetypes, RMSE and R² metrics are sensitive to outliers. Compare with cLHS (n=400) for stability."
  ) +
  theme_bw(base_size = 14) +
  theme(plot.title = element_text(face="bold", size=14),
         plot.subtitle = element_text(colour="grey25", size=11),
         plot.caption = element_text(size = 9, colour = "grey35", hjust = 0),
         axis.text.x = element_text(face="bold", size=11),
         axis.text.y = element_text(face="bold", size=13),
         panel.grid = element_blank())
ggsave(file.path(FIG, "fig_Shapley_3metrics_archetypes.png"),
        p, width = 9, height = 5, dpi = 300)
cli_alert_success("Saved fig_Shapley_3metrics_archetypes.png")

# ============================================================================
# Combined figure : cLHS + archetypes (6 columns)
# ============================================================================
comp_clhs <- fread(file.path(TAB, "tab_Shapley_3metrics_comparison.csv"))
comp_clhs[, scope := "cLHS\n(n=400)"]
comp[, scope := "Archetypes\n(n=4)"]
both <- rbind(comp_clhs, comp)
both_long <- melt(both, id.vars = c("variable", "scope"),
                    measure.vars = c("MAE_rank","RMSE_rank","R2_rank"),
                    variable.name = "metric", value.name = "Rank")
both_long[, metric := factor(metric,
                                levels = c("MAE_rank","RMSE_rank","R2_rank"),
                                labels = c("MAE","RMSE","R²"))]
both_long[, scope := factor(scope, levels = c("cLHS\n(n=400)","Archetypes\n(n=4)"))]
both_long[, variable := factor(variable, levels = VARS)]

p2 <- ggplot(both_long, aes(x = metric, y = variable, fill = factor(Rank))) +
  geom_tile(colour = "white", linewidth = 1) +
  geom_text(aes(label = Rank), colour = ifelse(both_long$Rank <= 2, "white", "grey15"),
             fontface = "bold", size = 5.5) +
  facet_grid(. ~ scope, switch = "x") +
  scale_fill_manual(values = pal, name = "Rank", drop = FALSE) +
  scale_y_discrete(limits = rev) +
  labs(
    title    = "Shapley attribution ranking across 3 metrics × 2 scales",
    subtitle = "MAE = mean(|φ_v|)     |     RMSE = sqrt(mean(Δ²)) global Shapley     |     R² = explained variance global Shapley",
    x = NULL, y = NULL,
    caption = "Consensus across 6 cells per variable indicates the robust rank. LAD #4 across all 6. LAI/fCover swap depending on metric scale."
  ) +
  theme_bw(base_size = 14) +
  theme(strip.background = element_rect(fill = "grey92", colour = NA),
         strip.text = element_text(face = "bold", size = 12),
         strip.placement = "outside",
         plot.title = element_text(face="bold", size=14),
         plot.subtitle = element_text(colour="grey25", size=10),
         plot.caption = element_text(size = 9, colour = "grey35", hjust = 0),
         axis.text.x = element_text(face="bold", size=11),
         axis.text.y = element_text(face="bold", size=13),
         panel.grid = element_blank())
ggsave(file.path(FIG, "fig_Shapley_3metrics_combined.png"),
        p2, width = 12, height = 5.5, dpi = 300)
cli_alert_success("Saved fig_Shapley_3metrics_combined.png")
