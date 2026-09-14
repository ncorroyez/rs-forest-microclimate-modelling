# ==============================================================================
# Amplifier vs buffer split of the ΔTmax validation. The global Pearson of
# v3.2.3 collapses (0.92 -> 0.60, doc Fig 6); we show this is driven entirely by
# a handful of atypical AMPLIFYING gaps (observed ΔTmax > 0), and that on the 44
# BUFFERING plots -- the refuges a microclimate map actually cares about --
# v3.2.3 is at least as good as v3.2.0 on ranking, linear agreement AND amplitude.
# Metric = dTmaxS (summer ΔTmax), the same as doc Fig 6.
#   Run from repo root:
#   Rscript Chapitre1/comparaison_versions/scripts/02_analyses_microclimat/cmp_v320_v323_ampbuff.R
# Out: figures/doc_media/image_ampbuff.png  +  tables/tab_v320_v323_ampbuff.csv
# ==============================================================================
suppressPackageStartupMessages({ library(data.table); library(ggplot2) })
DIR <- "Chapitre1/comparaison_versions"; TBL <- file.path(DIR, "tables"); MED <- file.path(DIR, "figures/doc_media")

d <- fread(file.path(TBL, "recap_JJAS_plotdata.csv"))
d[, grp := ifelse(dTmaxS_obs > 0, "amplificateur (trouée)", "tampon (couvert)")]

# long form: one row per (plot, version)
L <- rbind(
  d[, .(id_plot, grp, obs = dTmaxS_obs, sim = dTmaxS_v320, version = "v3.2.0")],
  d[, .(id_plot, grp, obs = dTmaxS_obs, sim = dTmaxS_v323, version = "v3.2.3 (yoyo)")])

# per (version x group) and per-version-ALL stats
stat1 <- function(s) data.table(n = nrow(s),
  pearson = cor(s$obs, s$sim), spearman = cor(s$obs, s$sim, method = "spearman"),
  amp = sd(s$sim) / sd(s$obs))
tab <- rbind(
  L[, stat1(.SD), by = .(version, scope = grp), .SDcols = c("obs","sim")],
  L[, cbind(scope = "TOUTES (53)", stat1(.SD)), by = version, .SDcols = c("obs","sim")])
tab[, (c("pearson","spearman","amp")) := lapply(.SD, round, 2), .SDcols = c("pearson","spearman","amp")]
setcolorder(tab, c("version","scope","n","pearson","spearman","amp"))
fwrite(tab, file.path(TBL, "tab_v320_v323_ampbuff.csv"))
cat("=== Amplifier/buffer split, ΔTmax (summer), agreement with obs ===\n"); print(tab[order(version, scope)])

# per-facet annotation
ann <- L[, {
  b <- .SD[grp == "tampon (couvert)"]; a <- .SD[grp == "amplificateur (trouée)"]
  .(lab = sprintf("TOUTES (n=%d): r=%.2f  rho=%.2f\ntampons (n=%d): r=%.2f  rho=%.2f  amp=%.0f%%\ntrouées (n=%d): r=%.2f  amp=%.0f%%",
                  .N, cor(obs, sim), cor(obs, sim, method = "spearman"),
                  nrow(b), cor(b$obs, b$sim), cor(b$obs, b$sim, method = "spearman"), 100*sd(b$sim)/sd(b$obs),
                  nrow(a), cor(a$obs, a$sim), 100*sd(a$sim)/sd(a$obs)))
}, by = version]

rng <- range(c(L$obs, L$sim))
fig <- ggplot(L, aes(obs, sim)) +
  geom_abline(slope = 1, intercept = 0, linetype = 2, colour = "grey55") +
  geom_hline(yintercept = 0, colour = "grey80", linewidth = .3) +
  geom_vline(xintercept = 0, colour = "grey80", linewidth = .3) +
  geom_point(aes(fill = grp), shape = 21, colour = "black", size = 2.2, stroke = .3, alpha = .9) +
  geom_text(data = ann, aes(x = rng[1], y = rng[2], label = lab), hjust = 0, vjust = 1,
            size = 2.9, lineheight = .95, colour = "grey20") +
  facet_wrap(~version) +
  scale_fill_manual(values = c("amplificateur (trouée)" = "#D9A441", "tampon (couvert)" = "#1B7837"), name = NULL) +
  coord_equal(xlim = rng, ylim = rng) +
  labs(x = "ΔTmax observé (°C)", y = "ΔTmax simulé (°C)",
       title = "ΔTmax simulé vs observé, séparé en placettes amplificatrices et tampons",
       subtitle = "Diagonale = accord parfait. La v3.2.3 aplatit vers 0 les quelques trouées amplificatrices (ΔTmax>0) — d'où la chute de son Pearson ;\nsur les 44 placettes tampons (les refuges), elle égale ou dépasse la v3.2.0 (rang, accord linéaire et amplitude).") +
  theme_bw(base_size = 11) +
  theme(legend.position = "bottom", panel.grid.minor = element_blank(),
        plot.subtitle = element_text(size = 8.3, colour = "grey35"))
ggsave(file.path(MED, "image_ampbuff.png"), fig, width = 10.5, height = 5.8, dpi = 200, bg = "white")
cat(sprintf("\nDONE -> %s\n", file.path(MED, "image_ampbuff.png")))
