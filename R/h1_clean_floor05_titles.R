# Regenerate the 2 existing figures (heatmap_2levels, HOBO_facet) WITHOUT
# any mention of "floor05" in titles/subtitles/captions.
suppressMessages({
  library(data.table); library(here); library(cli); library(ggplot2)
})

DATA <- here::here("outputs/lovb_floor05/data")
TAB  <- here::here("outputs/lovb_floor05/tables")
FIG  <- here::here("outputs/lovb_floor05/figures")

VARS <- c("LAI","fCover","Hmax","LAD")

# ----- Heatmap 2 levels -----
cli_h2("Heatmap 2 levels (no floor05 in titles)")
arch <- fread(file.path(TAB, "tab_attribution_archetypes_pooled_floor05.csv"))
clhs <- fread(file.path(TAB, "tab_attribution_cLHS_mean_floor05.csv"))
arch[, variable := factor(variable, levels = VARS)]
clhs[, variable := factor(variable, levels = VARS)]
methods <- c("LOVB_arch","LVA_arch","Shapley_arch",
              "LOVB_cLHS","LVA_cLHS","Shapley_cLHS")
labels  <- c("LOVB","LVA","Shapley","LOVB","LVA","Shapley")
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
  label = c("Level 1 — archetypes (4 K-means centroids)",
              "Level 2 — cLHS (400 plots)",
              "Consensus"))
p_hm <- ggplot(long_full, aes(x = col_idx, y = variable, fill = fill_rank)) +
  geom_tile(colour = "white", linewidth = 1.0) +
  geom_text(aes(label = ifelse(method == "Consensus",
                                 sprintf("%.1f", Rank),
                                 sprintf("%g", Rank))),
             colour = ifelse(as.numeric(long_full$fill_rank) <= 2,
                               "white", "grey15"),
             fontface = "bold", size = 5.5) +
  geom_vline(xintercept = c(3.5, 6.5), colour = "white", linewidth = 3) +
  geom_text(data = top_labels_df,
             aes(x = x, y = 4.85, label = label),
             inherit.aes = FALSE, fontface = "bold", size = 3.3,
             colour = "grey20") +
  scale_fill_manual(values = pal, name = "Rank") +
  scale_x_continuous(breaks = 1:7, labels = x_labels, expand = c(0,0)) +
  scale_y_discrete(limits = rev, expand = c(0,0)) +
  coord_cartesian(ylim = c(0.5, 5.1), clip = "off") +
  labs(
    title    = "Attribution ranking — 2 aggregation levels × 3 methods",
    subtitle = "Dark = rank 1 (strongest)  |  light = rank 4 (weakest)",
    x = NULL, y = NULL
  ) +
  theme_bw(base_size = 13) +
  theme(plot.title    = element_text(face = "bold", size = 14),
         plot.subtitle = element_text(colour = "grey25", size = 11),
         axis.text.x   = element_text(face = "bold", size = 10),
         axis.text.y   = element_text(face = "bold", size = 13),
         plot.margin   = margin(15, 10, 10, 10),
         panel.grid    = element_blank())
ggsave(file.path(FIG, "fig_heatmap_2levels.png"),
        p_hm, width = 12, height = 5.5, dpi = 300)
cli_alert_success("Refreshed fig_heatmap_2levels.png (no floor05 in title)")

# ----- HOBO facet -----
cli_h2("HOBO facet (no floor05 in titles)")
DT_h_f <- readRDS(file.path(DATA, "DT_contrib_HOBO_floor05.rds"))
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
tab_sp <- fread(file.path(TAB, "tab_HOBO_spearman_dTmax_floor05.csv"))
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
    title    = bquote("HOBO sensor-level validation (n = 53) — "
                        ~ rho[Spearman] ~ "vs observed" ~ Delta * T[max]),
    subtitle = "Rows : LOVB / LVA / Shapley     |     LOESS span = 0.9     |     Spearman ρ + bootstrap CI95 (1000 reps)",
    x = bquote(Delta * T[max]^"obs" ~ "(°C)"),
    y = bquote("Simulated  " * Delta[v] ~ "(°C)")
  ) +
  theme_bw(base_size = 14) +
  theme(strip.background = element_rect(fill = "grey92", colour = NA),
         strip.text       = element_text(face = "bold", size = 13),
         strip.placement  = "outside",
         plot.title       = element_text(face = "bold", size = 16),
         plot.subtitle    = element_text(colour = "grey25", size = 12),
         panel.grid.minor = element_blank())
ggsave(file.path(FIG, "fig_HOBO_facet.png"),
        p_hobo, width = 16, height = 13, dpi = 300)
cli_alert_success("Refreshed fig_HOBO_facet.png (no floor05 in title)")
