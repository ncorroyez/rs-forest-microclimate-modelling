# ==============================================================================
# A13 — Matrice de colinéarité des traits structuraux (53 placettes).
# Fonde le discours méthodo : tous les traits sont corrélés → il faut le Shapley
# (toutes coalitions) et non des partial plots ; explique pourquoi le COM proxy
# fCover/Hmax (A2) et le confondement cluster↔LAI.
# Traits : LAI, Hmax, fCover, COM (position verticale), VCI (complexité verticale).
#
# Sortie : outputs/figures_pipeline/annex/fig_trait_collinearity.{png,pdf}
# ==============================================================================
suppressMessages({ library(here); library(data.table); library(terra); library(tidyverse) })
source(here::here("scripts/_article_style.R"))
OUT <- here::here("outputs/figures_pipeline/annex"); dir.create(OUT, recursive = TRUE, showWarnings = FALSE)
OAK <- path.expand("~/Documents/NC_Full/03_RESULTS/Blois/Metrics/Oak_Only")

lad <- fread(here::here("in_files/lad_z05/Blois_lad_z05_r25.csv"))
ladc <- grep("^LAD_Layer_", names(lad), value = TRUE); hts <- as.numeric(gsub("LAD_Layer_","",ladc))
lad[, COM := apply(.SD, 1, function(r){ d<-as.numeric(r[ladc]); d[is.na(d)]<-0
  if(sum(d)<=0) NA_real_ else (sum(d*hts)/sum(d))/as.numeric(r[["Hmax"]]) }), .SDcols=c(ladc,"Hmax")]
fc <- fread(here::here("transfer/musica_version_benchmark/run_simulations/inputs/fcover_per_plot.csv"))
D  <- merge(lad[, .(id_plot, LAI, Hmax, COM)], fc, by = "id_plot")

# VCI extrait au point (complexité verticale)
vci_f <- file.path(OAK, "vci_res_10_m.tif")
if (file.exists(vci_f)) {
  pts <- vect(as.data.frame(lad[, .(x, y, id_plot)]), geom = c("x","y"), crs = "EPSG:32631")
  r <- rast(vci_f); p <- if (crs(r)!=crs(pts)) project(pts, crs(r)) else pts
  vci <- data.table(id_plot = lad$id_plot, VCI = terra::extract(r, p, ID = FALSE)[[1]])
  D <- merge(D, vci, by = "id_plot")
}
traits <- intersect(c("LAI","Hmax","fCover","COM","VCI"), names(D))
M <- cor(D[, ..traits], use = "pairwise.complete.obs")

# matrice complète, valeurs annotées (aucune paire perdue)
ML <- as.data.table(reshape2::melt(M, varnames = c("x","y"), value.name = "r"))
ML[, `:=`(x = factor(x, levels = traits), y = factor(y, levels = rev(traits)))]
p <- ggplot(ML, aes(x, y, fill = r)) +
  geom_tile(colour = "white", linewidth = 1) +
  geom_text(aes(label = sprintf("%.2f", r)), size = 4) +
  scale_fill_gradient2(low = "#2166AC", mid = "white", high = "#B2182B",
                       midpoint = 0, limits = c(-1, 1), name = "r") +
  coord_equal() +
  labs(x = NULL, y = NULL,
       subtitle = "Structural-trait collinearity (53 plots) — why exact Shapley, not partial plots") +
  theme_article() +
  theme(panel.border = element_blank(), axis.ticks = element_blank())
ggsave_article(file.path(OUT, "fig_trait_collinearity"), p, 6.2, 5.4)
cat("\n=== matrice de corrélation des traits ===\n"); print(round(M, 2))
cli::cli_alert_success("Saved fig_trait_collinearity")
