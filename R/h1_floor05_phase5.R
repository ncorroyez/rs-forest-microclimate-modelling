# ==============================================================================
# Phase 5 — Regenerate floor05 figures (heatmap, Type A/B, HOBO).
# Outputs in outputs/lovb_floor05/figures/.
# Plus a comparison table showing rank shifts vs original (non-floored).
# ==============================================================================

suppressMessages({
  library(data.table); library(here); library(cli); library(ggplot2)
  library(patchwork); library(cowplot)
})

OUT_DATA <- here::here("outputs/lovb_floor05/data")
OUT_TAB  <- here::here("outputs/lovb_floor05/tables")
OUT_FIG  <- here::here("outputs/lovb_floor05/figures")
dir.create(OUT_FIG, recursive = TRUE, showWarnings = FALSE)

VARS         <- c("LAI", "fCover", "Hmax", "LAD")
.CLUSTER_LABS <- c("1" = "C1 Open", "2" = "C2 Dense",
                    "3" = "C3 Low cover", "4" = "C4 Inter")
.PAL_CLUSTER  <- c("C1 Open" = "#E69F00", "C2 Dense" = "#0072B2",
                    "C3 Low cover" = "#009E73", "C4 Inter" = "#CC79A7")

cli_h1("Phase 5 — regenerate floor05 figures")
t_start <- Sys.time()

# ---- 1. Comparison table : ranks original vs floor05 -------------------------
cli_h2("Compare original ranks vs floor05 ranks")
orig_clhs <- fread(here::here("outputs/lovb/tables/tab_attribution_cLHS_mean.csv"))
new_clhs  <- fread(file.path(OUT_TAB, "tab_attribution_cLHS_mean_floor05.csv"))
orig_arch <- fread(here::here("outputs/lovb/tables/tab_attribution_archetypes_pooled.csv"))
new_arch  <- fread(file.path(OUT_TAB, "tab_attribution_archetypes_pooled_floor05.csv"))

cmp <- merge(
  orig_clhs[, .(variable, Shapley_cLHS = Shapley_cLHS_rank,
                  LOVB_cLHS = LOVB_cLHS_rank, LVA_cLHS = LVA_cLHS_rank)],
  new_clhs[, .(variable,
                 Shapley_cLHS_f05 = Shapley_cLHS_rank,
                 LOVB_cLHS_f05 = LOVB_cLHS_rank,
                 LVA_cLHS_f05 = LVA_cLHS_rank)],
  by = "variable"
)
cmp_arch <- merge(
  orig_arch[, .(variable, Shapley_arch = Shapley_arch_rank,
                  LOVB_arch = LOVB_arch_rank, LVA_arch = LVA_arch_rank)],
  new_arch[, .(variable,
                 Shapley_arch_f05 = Shapley_arch_rank,
                 LOVB_arch_f05 = LOVB_arch_rank,
                 LVA_arch_f05 = LVA_arch_rank)],
  by = "variable"
)
full_cmp <- merge(cmp_arch, cmp, by = "variable")
full_cmp[, variable := factor(variable, levels = VARS)]
setorder(full_cmp, variable)
fwrite(full_cmp, file.path(OUT_TAB, "tab_rank_shift_floor05_vs_orig.csv"))
cli_alert("Rank shifts (orig → floor05) :"); print(full_cmp)

# ---- 2. Heatmap 2 levels floor05 ---------------------------------------------
cli_h2("Heatmap 2 levels floor05")
arch <- fread(file.path(OUT_TAB, "tab_attribution_archetypes_pooled_floor05.csv"))
clhs <- fread(file.path(OUT_TAB, "tab_attribution_cLHS_mean_floor05.csv"))
arch[, variable := factor(variable, levels = VARS)]
clhs[, variable := factor(variable, levels = VARS)]

methods <- c("LOVB_arch", "LVA_arch", "Shapley_arch",
              "LOVB_cLHS", "LVA_cLHS", "Shapley_cLHS")
labels  <- c("LOVB", "LVA", "Shapley", "LOVB", "LVA", "Shapley")
long_list <- vector("list", length(methods))
for (i in seq_along(methods)) {
  m <- methods[i]
  src <- if (grepl("_arch$", m)) arch else clhs
  rank_col <- paste0(m, "_rank")
  long_list[[i]] <- src[, .(variable, method = labels[i], col_idx = i,
                              Rank = get(rank_col))]
}
long <- rbindlist(long_list)
cons <- long[, .(Rank = mean(Rank), col_idx = 7L), by = variable]
cons[, method := "Consensus"]
long_full <- rbind(long[, .(variable, method, col_idx, Rank)],
                    cons[, .(variable, method, col_idx, Rank)])
long_full[, fill_rank := as.character(round(Rank))]

pal <- c("1"="#440154","2"="#3B528B","3"="#5DC863","4"="#FDE725")
x_labels <- c("LOVB","LVA","Shapley","LOVB","LVA","Shapley","Consensus")
top_labels_df <- data.table(
  x = c(2, 5, 7),
  label = c("Level 1 — archetypes (4 K-means centroids, floor05)",
              "Level 2 — cLHS (400 plots, floor05)",
              "Consensus"))
sep_xs <- c(3.5, 6.5)

p_hm <- ggplot(long_full, aes(x = col_idx, y = variable, fill = fill_rank)) +
  geom_tile(colour = "white", linewidth = 1.0) +
  geom_text(aes(label = ifelse(method == "Consensus",
                                 sprintf("%.1f", Rank),
                                 sprintf("%g", Rank))),
             colour = ifelse(as.numeric(long_full$fill_rank) <= 2,
                               "white", "grey15"),
             fontface = "bold", size = 5.5) +
  geom_vline(xintercept = sep_xs, colour = "white", linewidth = 3) +
  geom_text(data = top_labels_df,
             aes(x = x, y = 4.85, label = label),
             inherit.aes = FALSE, fontface = "bold", size = 3.3,
             colour = "grey20") +
  scale_fill_manual(values = pal, name = "Rank") +
  scale_x_continuous(breaks = 1:7, labels = x_labels, expand = c(0,0)) +
  scale_y_discrete(limits = rev, expand = c(0,0)) +
  coord_cartesian(ylim = c(0.5, 5.1), clip = "off") +
  labs(
    title    = "Attribution ranking — floor05 (fCover floored to 0.5 where real < 0.5)",
    subtitle = "Dark = rank 1 (strongest)  |  light = rank 4 (weakest)  |  54 cLHS plots + 4 HOBO + C1 archetype affected",
    x = NULL, y = NULL,
    caption = "Floor05 sensitivity test : forçage fCover >= 0.5 partout. Comparer avec heatmap_2levels_v3 original pour mesurer impact."
  ) +
  theme_bw(base_size = 12) +
  theme(plot.title    = element_text(face = "bold", size = 14),
         plot.subtitle = element_text(colour = "grey25", size = 11),
         plot.caption  = element_text(size = 8, colour = "grey35", hjust = 0),
         axis.text.x   = element_text(face = "bold", size = 10),
         axis.text.y   = element_text(face = "bold", size = 13),
         plot.margin   = margin(15, 10, 10, 10),
         panel.grid    = element_blank())
ggsave(file.path(OUT_FIG, "fig_heatmap_2levels_floor05.png"),
        p_hm, width = 12, height = 5.5, dpi = 300)
cli_alert_success("Saved fig_heatmap_2levels_floor05.png")

# ---- 3. HOBO Spearman facet floor05 (3 methods x 4 vars) --------------------
cli_h2("HOBO Spearman facet floor05")
DT_h_f <- readRDS(file.path(OUT_DATA, "DT_contrib_HOBO_floor05.rds"))
cg     <- readRDS(here::here("outputs/lovb/data/DT_cross_gam_targets.rds"))
hobo   <- as.data.table(cg$HOBO)[, .(id_plot = id, dT = -dTmax_mean)]
DT_h_f <- merge(DT_h_f, hobo, by = "id_plot")
for (v in VARS) {
  DT_h_f[, paste0("LVA_", v) := get(paste0("Tmax_mean_LVA_", v)) - get("Tmax_mean_NULL")]
  DT_h_f[, paste0("LOVB_", v) := get(paste0("Delta_", v, "_mean"))]
  DT_h_f[, paste0("Shapley_", v) := (get(paste0("LVA_", v)) + get(paste0("LOVB_", v))) / 2]
}
long_hobo <- rbindlist(lapply(c("LOVB","LVA","Shapley"), function(m)
  rbindlist(lapply(VARS, function(v) data.table(
    method = m, Variable = v,
    dT = DT_h_f$dT,
    y  = DT_h_f[[paste0(m, "_", v)]])))))
long_hobo[, method   := factor(method, levels = c("LOVB","LVA","Shapley"))]
long_hobo[, Variable := factor(Variable, levels = VARS)]

tab_sp <- fread(file.path(OUT_TAB, "tab_HOBO_spearman_dTmax_floor05.csv"))
tab_sp[, Variable := factor(Variable, levels = VARS)]
tab_sp[, method   := factor(method, levels = c("LOVB","LVA","Shapley"))]
tab_sp[, rho_lab := sprintf("rho == %+.2f", rho_S)]
tab_sp[, p_lab   := ifelse(p_value < 1e-3, "italic(p) < 0.001",
                              sprintf("italic(p) == %.3f", p_value))]
tab_sp[, ci_lab  := sprintf("CI[95] *' '* '[' * %+.2f * '; ' * %+.2f * ']'",
                                ci_lo, ci_hi)]

p_hobo <- ggplot(long_hobo, aes(x = dT, y = y)) +
  geom_hline(yintercept = 0, linetype = "dotted", colour = "grey60") +
  geom_vline(xintercept = 0, linetype = "dotted", colour = "grey60") +
  geom_point(colour = "#2C5F2D", alpha = 0.7, size = 1.7) +
  geom_smooth(method = "loess", formula = y ~ x, se = TRUE,
                colour = "#FFB400", fill = "#FFB400",
                linewidth = 0.7, alpha = 0.18, span = 0.9) +
  geom_text(data = tab_sp, aes(x = -Inf, y = Inf, label = rho_lab),
             parse = TRUE, hjust = -0.08, vjust = 1.3, size = 4.5,
             fontface = "bold", inherit.aes = FALSE, colour = "grey15") +
  geom_text(data = tab_sp, aes(x = -Inf, y = Inf, label = p_lab),
             parse = TRUE, hjust = -0.08, vjust = 3.0, size = 3.8,
             inherit.aes = FALSE, colour = "grey25") +
  geom_text(data = tab_sp, aes(x = -Inf, y = Inf, label = ci_lab),
             parse = TRUE, hjust = -0.06, vjust = 5.0, size = 3.3,
             inherit.aes = FALSE, colour = "grey35") +
  facet_grid(method ~ Variable, scales = "free_y", switch = "y") +
  labs(
    title    = bquote("HOBO validation (n = 53, floor05)" ~ " — " ~ rho[Spearman]
                        ~ "vs observed" ~ Delta * T[max]),
    subtitle = "fCover floored to 0.5 where real < 0.5     |     Affects 4 HOBO sensors and the 8 fCover-real coalitions",
    x = bquote(Delta * T[max]^"obs" ~ "(°C)"),
    y = bquote("Simulated  " * Delta[v] ~ "(°C)"),
    caption = "Compare with fig_HOBO_facet_v3quater (original) to measure floor05 sensitivity. Sign convention : positive Δ_v ⇒ v contributes to buffering."
  ) +
  theme_bw(base_size = 14) +
  theme(strip.background = element_rect(fill = "grey92", colour = NA),
         strip.text       = element_text(face = "bold", size = 13),
         strip.placement  = "outside",
         plot.title       = element_text(face = "bold", size = 16),
         plot.subtitle    = element_text(colour = "grey25", size = 12),
         plot.caption     = element_text(size = 10, colour = "grey35", hjust = 0),
         panel.grid.minor = element_blank())
ggsave(file.path(OUT_FIG, "fig_HOBO_facet_floor05.png"),
        p_hobo, width = 16, height = 13, dpi = 300)
cli_alert_success("Saved fig_HOBO_facet_floor05.png")

# ---- 4. Quick LOVB attribution per cluster floor05 ---------------------------
cli_h2("LOVB attribution per cluster floor05")
DT_c <- readRDS(file.path(OUT_DATA, "DT_contrib_cLHS_floor05.rds"))
long_c <- rbindlist(lapply(VARS, function(v) data.table(
  Cluster = DT_c$Cluster, variable = v,
  delta_LOVB = DT_c[[paste0("Delta_", v, "_mean")]],
  LVA = DT_c[[paste0("LVA_", v, "_mean")]]
)))
long_c[, variable := factor(variable, levels = VARS)]
long_c[, Cluster_lab := .CLUSTER_LABS[as.character(Cluster)]]
long_c[, Cluster_lab := factor(Cluster_lab, levels = .CLUSTER_LABS)]
# Median per cluster x variable for both LOVB and LVA
med_summary <- long_c[, .(LOVB_median = median(delta_LOVB, na.rm = TRUE),
                            LVA_median = median(LVA, na.rm = TRUE)),
                       by = .(Cluster_lab, variable)]
fwrite(med_summary, file.path(OUT_TAB, "tab_LOVB_LVA_median_per_cluster_floor05.csv"))
cli_alert("Median LOVB/LVA per cluster x variable :"); print(med_summary)

cli_alert_success("Phase 5 done in {round(as.numeric(difftime(Sys.time(), t_start, units='mins')),1)} min")
