# ==============================================================================
# 8 figures Type A/B × cLHS/archetypes × mean/daily.
# Naming convention :
#   fig_typeA_cLHS_mean.png            scatter Tmax_REF vs Tmax_LOVB_v (400 plots)
#   fig_typeA_cLHS_daily.png           idem but 122 days per plot (hexbin)
#   fig_typeA_archetypes_mean.png      idem on 4 archetypes
#   fig_typeA_archetypes_daily.png     4 archetypes × 122 days
#   fig_typeB_cLHS_mean.png            scatter Δ_v vs trait value (400 plots)
#   fig_typeB_cLHS_daily.png           Δ_v daily vs trait value
#   fig_typeB_archetypes_mean.png      4 archetypes
#   fig_typeB_archetypes_daily.png     4 archetypes × 122 days
#
# No mention of "floor05" in any title/subtitle/caption.
# Data from outputs/lovb_floor05/data/.
# ==============================================================================

suppressMessages({
  library(data.table); library(here); library(cli); library(ggplot2)
})

DATA <- here::here("outputs/lovb_floor05/data")
FIG  <- here::here("outputs/lovb_floor05/figures")
dir.create(FIG, recursive = TRUE, showWarnings = FALSE)

VARS         <- c("LAI", "fCover", "Hmax", "LAD")
.CLUSTER_LABS <- c("1" = "C1 Open", "2" = "C2 Dense",
                    "3" = "C3 Low cover", "4" = "C4 Inter")
.PAL_CLUSTER  <- c("C1 Open" = "#E69F00", "C2 Dense" = "#0072B2",
                    "C3 Low cover" = "#009E73", "C4 Inter" = "#CC79A7")
.CLUSTER_ARCH <- c("Arch_C1" = "C1 Open", "Arch_C2" = "C2 Dense",
                    "Arch_C3" = "C3 Low cover", "Arch_C4" = "C4 Inter")

trait_map <- c(LAI = "LAI", fCover = "fCover", Hmax = "Hmax", LAD = "FPC1")

# ============================================================================
# Helpers : load data
# ============================================================================
load_clhs_mean <- function() {
  DT <- readRDS(file.path(DATA, "DT_contrib_cLHS_floor05.rds"))
  DT
}
load_clhs_daily <- function() {
  DT <- readRDS(file.path(DATA, "DT_daily_cLHS_floor05.rds"))
  DT
}
load_arch_mean <- function() {
  DT <- readRDS(file.path(DATA, "DT_contrib_archetypes_floor05.rds"))
  DT
}
load_arch_daily <- function() {
  DT <- readRDS(file.path(DATA, "DT_daily_archetypes_floor05.rds"))
  DT
}

# ============================================================================
# Type A : Tmax_REF vs Tmax_LOVB_v
# ============================================================================
build_typeA_clhs_mean <- function() {
  cli_h2("Type A cLHS mean")
  DT <- load_clhs_mean()
  long <- rbindlist(lapply(VARS, function(v) data.table(
    Cluster = DT$Cluster, Variable = v,
    Tmax_REF  = DT$Tmax_mean_REF,
    Tmax_LOVB = DT[[paste0("Tmax_mean_LOVB_", v)]])))
  long[, Variable := factor(Variable, levels = VARS)]
  long[, Cluster_lab := .CLUSTER_LABS[as.character(Cluster)]]
  long[, Cluster_lab := factor(Cluster_lab, levels = .CLUSTER_LABS)]
  pal <- setNames(unname(.PAL_CLUSTER), .CLUSTER_LABS)
  ann <- long[, .(rmse = round(sqrt(mean((Tmax_LOVB - Tmax_REF)^2, na.rm = TRUE)), 2)),
                by = Variable]
  ann[, label := paste0("RMSE == ", sprintf("%.2f", rmse), " * '°C'")]

  p <- ggplot(long, aes(x = Tmax_REF, y = Tmax_LOVB, colour = Cluster_lab)) +
    geom_abline(slope = 1, intercept = 0, linetype = "dashed", colour = "grey50") +
    geom_point(alpha = 0.55, size = 1.6) +
    geom_text(data = ann, aes(x = -Inf, y = Inf, label = label),
               parse = TRUE, hjust = -0.06, vjust = 1.4, size = 4,
               fontface = "bold", inherit.aes = FALSE, colour = "grey15") +
    facet_wrap(~ Variable, nrow = 1L, scales = "free") +
    scale_colour_manual(values = pal, name = NULL, drop = FALSE) +
    labs(
      title    = bquote("Type A — " ~ T[max]^"LOVB,v" ~ "vs" ~ T[max]^"REF" ~ "(cLHS, MEAN over 122 d)"),
      subtitle = "n = 400 plots     |     Dashed line = 1:1 identity     |     Deviation = v's contribution",
      x = bquote(T[max] ~ "buffering REF" ~ "(°C)"),
      y = bquote(T[max] ~ "buffering LOVB_v" ~ "(°C)")
    ) +
    theme_bw(base_size = 14) +
    theme(strip.background = element_rect(fill="grey92", colour=NA),
           strip.text = element_text(face="bold", size=13),
           plot.title = element_text(face="bold", size=15),
           plot.subtitle = element_text(colour="grey25", size=11),
           legend.position = "bottom",
           panel.grid.minor = element_blank())
  ggsave(file.path(FIG, "fig_typeA_cLHS_mean.png"), p,
          width = 16, height = 5.5, dpi = 300)
  cli_alert_success("Saved fig_typeA_cLHS_mean.png")
}

build_typeA_clhs_daily <- function() {
  cli_h2("Type A cLHS daily")
  DT_d <- load_clhs_daily()
  ref <- DT_d[bit_code == "1111", .(x, y, date, Tmax_REF = Delta_Tmax)]
  bit_map <- c(LAI = "0111", Hmax = "1011", fCover = "1101", LAD = "1110")
  long_list <- list()
  for (v in VARS) {
    lovb <- DT_d[bit_code == bit_map[v],
                   .(x, y, date, Tmax_LOVB = Delta_Tmax)]
    m <- merge(ref, lovb, by = c("x","y","date"))
    m[, Variable := v]
    long_list[[v]] <- m
  }
  long <- rbindlist(long_list)
  long[, Variable := factor(Variable, levels = VARS)]
  # Add cluster info
  cl <- as.data.table(readRDS(here::here("out_files/Sensitivity_Analysis/clhs_sample_floor05.rds")))[, .(x, y, Cluster)]
  long <- merge(long, cl, by = c("x","y"))
  long[, Cluster_lab := .CLUSTER_LABS[as.character(Cluster)]]
  long[, Cluster_lab := factor(Cluster_lab, levels = .CLUSTER_LABS)]
  pal <- setNames(unname(.PAL_CLUSTER), .CLUSTER_LABS)
  ann <- long[, .(rmse = round(sqrt(mean((Tmax_LOVB - Tmax_REF)^2, na.rm = TRUE)), 2),
                    n_obs = .N),
                by = Variable]
  ann[, label := paste0("RMSE == ", sprintf("%.2f", rmse), " * '°C   n='",
                          " * '", n_obs, "'")]

  p <- ggplot(long, aes(x = Tmax_REF, y = Tmax_LOVB)) +
    geom_abline(slope = 1, intercept = 0, linetype = "dashed", colour = "grey50") +
    geom_point(aes(colour = Cluster_lab), alpha = 0.06, size = 0.4) +
    geom_text(data = ann, aes(x = -Inf, y = Inf, label = label),
               parse = TRUE, hjust = -0.06, vjust = 1.4, size = 4,
               fontface = "bold", inherit.aes = FALSE, colour = "grey15") +
    facet_wrap(~ Variable, nrow = 1L, scales = "free") +
    scale_colour_manual(values = pal, name = NULL, drop = FALSE) +
    guides(colour = guide_legend(override.aes = list(alpha = 1, size = 3))) +
    labs(
      title    = bquote("Type A — " ~ T[max]^"LOVB,v" ~ "vs" ~ T[max]^"REF" ~ "(cLHS, DAILY)"),
      subtitle = "n = 400 plots × 122 days = ~48,800 pairs per panel     |     Daily values, no averaging     |     Dashed line = 1:1 identity",
      x = bquote(Delta * T[max] ~ "buffering REF, daily (°C)"),
      y = bquote(Delta * T[max] ~ "buffering LOVB_v, daily (°C)")
    ) +
    theme_bw(base_size = 14) +
    theme(strip.background = element_rect(fill="grey92", colour=NA),
           strip.text = element_text(face="bold", size=13),
           plot.title = element_text(face="bold", size=15),
           plot.subtitle = element_text(colour="grey25", size=11),
           legend.position = "bottom",
           panel.grid.minor = element_blank())
  ggsave(file.path(FIG, "fig_typeA_cLHS_daily.png"), p,
          width = 16, height = 5.5, dpi = 300)
  cli_alert_success("Saved fig_typeA_cLHS_daily.png")
}

build_typeA_arch_mean <- function() {
  cli_h2("Type A archetypes mean")
  DT <- load_arch_mean()
  long <- rbindlist(lapply(VARS, function(v) data.table(
    archetype = DT$archetype, Variable = v,
    Tmax_REF  = DT$Tmax_mean_REF,
    Tmax_LOVB = DT[[paste0("Tmax_mean_LOVB_", v)]])))
  long[, Variable := factor(Variable, levels = VARS)]
  long[, archetype_lab := .CLUSTER_ARCH[archetype]]
  long[, archetype_lab := factor(archetype_lab, levels = .CLUSTER_ARCH)]
  pal <- setNames(unname(.PAL_CLUSTER), .CLUSTER_ARCH)
  ann <- long[, .(rmse = round(sqrt(mean((Tmax_LOVB - Tmax_REF)^2, na.rm = TRUE)), 2)),
                by = Variable]
  ann[, label := paste0("RMSE == ", sprintf("%.2f", rmse), " * '°C'")]

  p <- ggplot(long, aes(x = Tmax_REF, y = Tmax_LOVB, colour = archetype_lab)) +
    geom_abline(slope = 1, intercept = 0, linetype = "dashed", colour = "grey50") +
    geom_point(size = 4) +
    geom_text(aes(label = archetype_lab), nudge_y = 0.05, size = 3, vjust = 0,
                colour = "grey20", fontface = "bold", show.legend = FALSE) +
    geom_text(data = ann, aes(x = -Inf, y = Inf, label = label),
               parse = TRUE, hjust = -0.06, vjust = 1.4, size = 4,
               fontface = "bold", inherit.aes = FALSE, colour = "grey15") +
    facet_wrap(~ Variable, nrow = 1L, scales = "free") +
    scale_colour_manual(values = pal, name = NULL, drop = FALSE) +
    labs(
      title    = bquote("Type A — " ~ T[max]^"LOVB,v" ~ "vs" ~ T[max]^"REF" ~ "(archetypes, MEAN over 122 d)"),
      subtitle = "n = 4 K-means centroids     |     Dashed line = 1:1 identity",
      x = bquote(T[max] ~ "buffering REF" ~ "(°C)"),
      y = bquote(T[max] ~ "buffering LOVB_v" ~ "(°C)")
    ) +
    theme_bw(base_size = 14) +
    theme(strip.background = element_rect(fill="grey92", colour=NA),
           strip.text = element_text(face="bold", size=13),
           plot.title = element_text(face="bold", size=15),
           plot.subtitle = element_text(colour="grey25", size=11),
           legend.position = "bottom",
           panel.grid.minor = element_blank())
  ggsave(file.path(FIG, "fig_typeA_archetypes_mean.png"), p,
          width = 16, height = 5.5, dpi = 300)
  cli_alert_success("Saved fig_typeA_archetypes_mean.png")
}

build_typeA_arch_daily <- function() {
  cli_h2("Type A archetypes daily")
  DT_d <- load_arch_daily()
  ref <- DT_d[bit_code == "1111", .(archetype, date, Tmax_REF = Delta_Tmax)]
  bit_map <- c(LAI = "0111", Hmax = "1011", fCover = "1101", LAD = "1110")
  long_list <- list()
  for (v in VARS) {
    lovb <- DT_d[bit_code == bit_map[v],
                   .(archetype, date, Tmax_LOVB = Delta_Tmax)]
    m <- merge(ref, lovb, by = c("archetype","date"))
    m[, Variable := v]
    long_list[[v]] <- m
  }
  long <- rbindlist(long_list)
  long[, Variable := factor(Variable, levels = VARS)]
  long[, archetype_lab := .CLUSTER_ARCH[archetype]]
  long[, archetype_lab := factor(archetype_lab, levels = .CLUSTER_ARCH)]
  pal <- setNames(unname(.PAL_CLUSTER), .CLUSTER_ARCH)
  ann <- long[, .(rmse = round(sqrt(mean((Tmax_LOVB - Tmax_REF)^2, na.rm = TRUE)), 2),
                    n_obs = .N),
                by = Variable]
  ann[, label := paste0("RMSE == ", sprintf("%.2f", rmse), " * '°C   n='",
                          " * '", n_obs, "'")]

  p <- ggplot(long, aes(x = Tmax_REF, y = Tmax_LOVB, colour = archetype_lab)) +
    geom_abline(slope = 1, intercept = 0, linetype = "dashed", colour = "grey50") +
    geom_point(alpha = 0.45, size = 1.2) +
    geom_text(data = ann, aes(x = -Inf, y = Inf, label = label),
               parse = TRUE, hjust = -0.06, vjust = 1.4, size = 4,
               fontface = "bold", inherit.aes = FALSE, colour = "grey15") +
    facet_wrap(~ Variable, nrow = 1L, scales = "free") +
    scale_colour_manual(values = pal, name = NULL, drop = FALSE) +
    labs(
      title    = bquote("Type A — " ~ T[max]^"LOVB,v" ~ "vs" ~ T[max]^"REF" ~ "(archetypes, DAILY)"),
      subtitle = "n = 4 archetypes × 122 days = 488 daily pairs per panel     |     Daily values, no averaging",
      x = bquote(Delta * T[max] ~ "buffering REF, daily (°C)"),
      y = bquote(Delta * T[max] ~ "buffering LOVB_v, daily (°C)")
    ) +
    theme_bw(base_size = 14) +
    theme(strip.background = element_rect(fill="grey92", colour=NA),
           strip.text = element_text(face="bold", size=13),
           plot.title = element_text(face="bold", size=15),
           plot.subtitle = element_text(colour="grey25", size=11),
           legend.position = "bottom",
           panel.grid.minor = element_blank())
  ggsave(file.path(FIG, "fig_typeA_archetypes_daily.png"), p,
          width = 16, height = 5.5, dpi = 300)
  cli_alert_success("Saved fig_typeA_archetypes_daily.png")
}

# ============================================================================
# Type B : Δ_v vs trait value
# ============================================================================
build_typeB_clhs_mean <- function() {
  cli_h2("Type B cLHS mean")
  DT <- load_clhs_mean()
  long <- rbindlist(lapply(VARS, function(v) data.table(
    Cluster = DT$Cluster, Variable = v,
    trait = DT[[trait_map[v]]],
    delta = DT[[paste0("Delta_", v, "_mean")]])))
  long[, Variable := factor(Variable, levels = VARS)]
  long[, Cluster_lab := .CLUSTER_LABS[as.character(Cluster)]]
  long[, Cluster_lab := factor(Cluster_lab, levels = .CLUSTER_LABS)]
  pal <- setNames(unname(.PAL_CLUSTER), .CLUSTER_LABS)

  p <- ggplot(long, aes(x = trait, y = delta)) +
    geom_hline(yintercept = 0, linetype = "dashed", colour = "grey40") +
    geom_point(aes(colour = Cluster_lab), alpha = 0.55, size = 1.7) +
    geom_smooth(method = "loess", span = 0.75, se = TRUE,
                  colour = "grey20", fill = "grey70", alpha = 0.3) +
    facet_wrap(~ Variable, nrow = 1L, scales = "free_x") +
    scale_colour_manual(values = pal, name = NULL, drop = FALSE) +
    labs(
      title    = bquote("Type B — " ~ Delta[v]^"LOVB" ~ "vs trait value (cLHS, MEAN over 122 d)"),
      subtitle = "n = 400 plots     |     LOESS span = 0.75 with 95%CI     |     LAD trait = FPC1 score",
      x = "Trait value",
      y = bquote(Delta[v]^"LOVB" ~ "(°C)")
    ) +
    theme_bw(base_size = 14) +
    theme(strip.background = element_rect(fill="grey92", colour=NA),
           strip.text = element_text(face="bold", size=13),
           plot.title = element_text(face="bold", size=15),
           plot.subtitle = element_text(colour="grey25", size=11),
           legend.position = "bottom",
           panel.grid.minor = element_blank())
  ggsave(file.path(FIG, "fig_typeB_cLHS_mean.png"), p,
          width = 16, height = 5.5, dpi = 300)
  cli_alert_success("Saved fig_typeB_cLHS_mean.png")
}

build_typeB_clhs_daily <- function() {
  cli_h2("Type B cLHS daily")
  DT_d <- load_clhs_daily()
  ref <- DT_d[bit_code == "1111", .(x, y, date, Tmax_REF = Delta_Tmax)]
  bit_map <- c(LAI = "0111", Hmax = "1011", fCover = "1101", LAD = "1110")
  long_list <- list()
  for (v in VARS) {
    lovb <- DT_d[bit_code == bit_map[v],
                   .(x, y, date, Tmax_LOVB = Delta_Tmax)]
    m <- merge(ref, lovb, by = c("x","y","date"))
    m[, delta := Tmax_REF - Tmax_LOVB]
    m[, Variable := v]
    long_list[[v]] <- m[, .(x, y, date, delta, Variable)]
  }
  long <- rbindlist(long_list)
  # Merge with traits (constant per plot)
  cl <- as.data.table(readRDS(here::here("out_files/Sensitivity_Analysis/clhs_sample_floor05.rds")))[, .(x, y, Cluster, LAI, Hmax, fCover, FPC1)]
  long <- merge(long, cl, by = c("x","y"))
  long[, trait := NA_real_]
  for (v in VARS) long[Variable == v, trait := get(trait_map[v])]
  long[, Variable := factor(Variable, levels = VARS)]
  long[, Cluster_lab := .CLUSTER_LABS[as.character(Cluster)]]
  long[, Cluster_lab := factor(Cluster_lab, levels = .CLUSTER_LABS)]
  pal <- setNames(unname(.PAL_CLUSTER), .CLUSTER_LABS)

  p <- ggplot(long, aes(x = trait, y = delta)) +
    geom_hline(yintercept = 0, linetype = "dashed", colour = "grey40") +
    geom_point(aes(colour = Cluster_lab), alpha = 0.06, size = 0.4) +
    facet_wrap(~ Variable, nrow = 1L, scales = "free_x") +
    scale_colour_manual(values = pal, name = NULL, drop = FALSE) +
    guides(colour = guide_legend(override.aes = list(alpha = 1, size = 3))) +
    labs(
      title    = bquote("Type B — " ~ Delta[v]^"LOVB" ~ "vs trait value (cLHS, DAILY)"),
      subtitle = "n = 400 plots × 122 days = ~48,800 daily values per panel     |     Trait value constant per plot, Δ_v varies day-by-day",
      x = "Trait value (constant per plot)",
      y = bquote(Delta[v]^"LOVB" ~ "daily (°C)")
    ) +
    theme_bw(base_size = 14) +
    theme(strip.background = element_rect(fill="grey92", colour=NA),
           strip.text = element_text(face="bold", size=13),
           plot.title = element_text(face="bold", size=15),
           plot.subtitle = element_text(colour="grey25", size=11),
           legend.position = "bottom",
           panel.grid.minor = element_blank())
  ggsave(file.path(FIG, "fig_typeB_cLHS_daily.png"), p,
          width = 16, height = 5.5, dpi = 300)
  cli_alert_success("Saved fig_typeB_cLHS_daily.png")
}

build_typeB_arch_mean <- function() {
  cli_h2("Type B archetypes mean")
  DT <- load_arch_mean()
  # Need traits per archetype (centroid means from df_sample_floor)
  df_sample <- as.data.table(readRDS(here::here("out_files/Sensitivity_Analysis/clhs_sample_floor05.rds")))
  arch_traits <- df_sample[, .(LAI = mean(LAI), Hmax = mean(Hmax),
                                  fCover = mean(fCover), FPC1 = mean(FPC1)),
                              by = Cluster]
  arch_traits[, archetype := paste0("Arch_C", Cluster)]
  setkey(arch_traits, archetype)

  long_list <- list()
  for (v in VARS) {
    DT_v <- DT[, .(archetype, delta = get(paste0("Delta_", v, "_mean")))]
    DT_v <- merge(DT_v, arch_traits, by = "archetype")
    DT_v[, trait := get(trait_map[v])]
    DT_v[, Variable := v]
    long_list[[v]] <- DT_v[, .(archetype, Variable, trait, delta)]
  }
  long <- rbindlist(long_list)
  long[, Variable := factor(Variable, levels = VARS)]
  long[, archetype_lab := .CLUSTER_ARCH[archetype]]
  long[, archetype_lab := factor(archetype_lab, levels = .CLUSTER_ARCH)]
  pal <- setNames(unname(.PAL_CLUSTER), .CLUSTER_ARCH)

  p <- ggplot(long, aes(x = trait, y = delta, colour = archetype_lab)) +
    geom_hline(yintercept = 0, linetype = "dashed", colour = "grey40") +
    geom_point(size = 4) +
    geom_text(aes(label = archetype_lab), nudge_y = 0.05, size = 3, vjust = 0,
                colour = "grey20", fontface = "bold", show.legend = FALSE) +
    facet_wrap(~ Variable, nrow = 1L, scales = "free_x") +
    scale_colour_manual(values = pal, name = NULL, drop = FALSE) +
    labs(
      title    = bquote("Type B — " ~ Delta[v]^"LOVB" ~ "vs trait value (archetypes, MEAN over 122 d)"),
      subtitle = "n = 4 K-means centroids",
      x = "Trait value (cluster mean)",
      y = bquote(Delta[v]^"LOVB" ~ "(°C)")
    ) +
    theme_bw(base_size = 14) +
    theme(strip.background = element_rect(fill="grey92", colour=NA),
           strip.text = element_text(face="bold", size=13),
           plot.title = element_text(face="bold", size=15),
           plot.subtitle = element_text(colour="grey25", size=11),
           legend.position = "bottom",
           panel.grid.minor = element_blank())
  ggsave(file.path(FIG, "fig_typeB_archetypes_mean.png"), p,
          width = 16, height = 5.5, dpi = 300)
  cli_alert_success("Saved fig_typeB_archetypes_mean.png")
}

build_typeB_arch_daily <- function() {
  cli_h2("Type B archetypes daily")
  DT_d <- load_arch_daily()
  df_sample <- as.data.table(readRDS(here::here("out_files/Sensitivity_Analysis/clhs_sample_floor05.rds")))
  arch_traits <- df_sample[, .(LAI = mean(LAI), Hmax = mean(Hmax),
                                  fCover = mean(fCover), FPC1 = mean(FPC1)),
                              by = Cluster]
  arch_traits[, archetype := paste0("Arch_C", Cluster)]

  ref <- DT_d[bit_code == "1111", .(archetype, date, Tmax_REF = Delta_Tmax)]
  bit_map <- c(LAI = "0111", Hmax = "1011", fCover = "1101", LAD = "1110")
  long_list <- list()
  for (v in VARS) {
    lovb <- DT_d[bit_code == bit_map[v],
                   .(archetype, date, Tmax_LOVB = Delta_Tmax)]
    m <- merge(ref, lovb, by = c("archetype","date"))
    m[, delta := Tmax_REF - Tmax_LOVB]
    m[, Variable := v]
    m <- merge(m, arch_traits, by = "archetype")
    m[, trait := get(trait_map[v])]
    long_list[[v]] <- m[, .(archetype, date, Variable, trait, delta)]
  }
  long <- rbindlist(long_list)
  long[, Variable := factor(Variable, levels = VARS)]
  long[, archetype_lab := .CLUSTER_ARCH[archetype]]
  long[, archetype_lab := factor(archetype_lab, levels = .CLUSTER_ARCH)]
  pal <- setNames(unname(.PAL_CLUSTER), .CLUSTER_ARCH)

  p <- ggplot(long, aes(x = trait, y = delta, colour = archetype_lab)) +
    geom_hline(yintercept = 0, linetype = "dashed", colour = "grey40") +
    geom_point(alpha = 0.4, size = 1.3) +
    facet_wrap(~ Variable, nrow = 1L, scales = "free_x") +
    scale_colour_manual(values = pal, name = NULL, drop = FALSE) +
    labs(
      title    = bquote("Type B — " ~ Delta[v]^"LOVB" ~ "vs trait value (archetypes, DAILY)"),
      subtitle = "n = 4 archetypes × 122 days = 488 daily values per panel",
      x = "Trait value (constant per archetype)",
      y = bquote(Delta[v]^"LOVB" ~ "daily (°C)")
    ) +
    theme_bw(base_size = 14) +
    theme(strip.background = element_rect(fill="grey92", colour=NA),
           strip.text = element_text(face="bold", size=13),
           plot.title = element_text(face="bold", size=15),
           plot.subtitle = element_text(colour="grey25", size=11),
           legend.position = "bottom",
           panel.grid.minor = element_blank())
  ggsave(file.path(FIG, "fig_typeB_archetypes_daily.png"), p,
          width = 16, height = 5.5, dpi = 300)
  cli_alert_success("Saved fig_typeB_archetypes_daily.png")
}

if (sys.nframe() == 0) {
  cli_h1("8 figures Type A/B × cLHS/archetypes × mean/daily")
  build_typeA_clhs_mean()
  build_typeA_clhs_daily()
  build_typeA_arch_mean()
  build_typeA_arch_daily()
  build_typeB_clhs_mean()
  build_typeB_clhs_daily()
  build_typeB_arch_mean()
  build_typeB_arch_daily()
  cli_alert_success("All 8 figures in {.path {FIG}}")
}
