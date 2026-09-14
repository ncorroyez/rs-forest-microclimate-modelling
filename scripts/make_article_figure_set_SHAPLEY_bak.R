# ==============================================================================
# Jeu de figures NUMÉROTÉ façon article (Chap 1) — copie les figures déjà
# produites vers outputs/figures_article_<branch>/ avec la numérotation narrative
# du plan (review/plan_article_chap1.md). Réordonner = éditer MAIN / ANNEX.
#
# Source = branche z05 (v3.2.0, article) + annexes partagées (figures_pipeline/annex).
# Sortie = Fig01_*, … (corps) et A1_*, … (annexes) en png + pdf si dispo,
#          + manifeste tab_article_figures.csv, + liste des manquantes.
#   Rscript chapter1/ledger/make_article_figure_set.R
# ==============================================================================

suppressMessages({ library(here) })
BRANCH <- Sys.getenv("PIPE_BRANCH", "z05")
P  <- here::here(sprintf("outputs/figures_pipeline_%s", BRANCH))   # figures de la branche
PA <- here::here("outputs/figures_pipeline/annex")                  # annexes partagées
OUT <- here::here(sprintf("outputs/figures_article_%s", BRANCH))
unlink(OUT, recursive = TRUE); dir.create(OUT, recursive = TRUE, showWarnings = FALSE)

# ---- ORDRE (éditable) : list(id, titre, source_sans_extension) ---------------
MAIN <- list(
  list("Fig1", "typology_archetype_profiles", here::here("outputs/figs_MEB2026_final/fig_archetypes_profiles")),  # mean+ribbon+traits (style MEB)
  list("Fig2", "shapley_ranking",          file.path(P,  "fig_shapley_ranking")),
  list("Fig3", "shapley_by_cluster",       file.path(P,  "fig_shapley_by_cluster")),
  list("Fig4", "validation_glass_ceiling", file.path(P,  "fig_validation_perplot")),
  list("Fig5", "forward_validation",       file.path(P,  "fig_validation_forward_shapley")),
  list("Fig6", "vertical_gradient",        file.path(PA, "fig_vertical_profiles_legacy"))
)
## Annexes consolidées : A6+A7 fusionnées (A6a/b/c), A10+A11 fusionnées (A8a/b),
## ancienne A8 (LAD z05 vs z1) supprimée. Renumérotation propre A1→A10.
ANNEX <- list(
  list("A1a", "fpca_loadings",             file.path(P,  "fig_fpca_loadings")),     # FPCs du LAD normalisé x&y -> clustering
  list("A1b", "fpca_harmonics",            file.path(P,  "fig_fpca_harmonics")),
  list("A1c", "kmeans_elbow",              file.path(P,  "fig_kmeans_elbow")),
  list("A1d", "cluster_map",               here::here("outputs/figs_MEB2026_final/fig_map_clusters_blois")),
  list("A2",  "h2_field_validation",       file.path(PA, "fig_h2_field_validation")),
  list("A3",  "h2_controlled_insilico",    file.path(PA, "fig_h2_controlled_topheavy")),
  list("A4",  "conditional_shapley_LAD",   file.path(PA, "fig_annex_conditional_shapley_LAD")),
  list("A4b", "conditional_shapley_LAI_fCover", file.path(P, "annex/fig_annex_conditional_shapley_lai_fcover")),
  list("A5",  "shapley_extra_metrics",     file.path(PA, "fig_annex_shapley_extra_by_cluster")),  # inclut ΔTmin (cf. texte)
  list("A6a", "version_benchmark_scatter", file.path(PA, "fig_musica_version_benchmark_scatter")),  # A6 = robustesse version (A6+A7 fusionnées)
  list("A6b", "version_benchmark_heatmap", file.path(PA, "fig_musica_version_benchmark_heatmap")),
  list("A6c", "version_sim_vs_obs",        file.path(PA, "fig_sim_vs_obs_v320_v323_1m")),
  list("A7",  "sentinel2_optical_illusion", file.path(PA, "fig_annex_s2_optical_illusion")),  # pont Ch.2
  list("A8",  "glass_ceiling",             file.path(PA, "fig_A8_glass_ceiling")),           # composite (a topo + b pente)
  list("A9",  "height_sampling_robustness", here::here("outputs/figs_MEB2026_final/fig_loo_height_compare_tmax")),
  list("A10", "trait_collinearity",         file.path(PA, "fig_trait_collinearity"))
)

# ---- copie numérotée ---------------------------------------------------------
copy_set <- function(items, kind) {
  do.call(rbind, lapply(items, function(it) {
    id <- it[[1]]; title <- it[[2]]; src <- it[[3]]
    got <- character(0)
    for (ext in c("png")) {                  # PNG only pour l'instant (pas de PDF)
      s <- paste0(src, ".", ext)
      if (file.exists(s)) {
        d <- file.path(OUT, sprintf("%s_%s.%s", id, title, ext))
        file.copy(s, d, overwrite = TRUE); got <- c(got, ext)
      }
    }
    data.frame(kind = kind, id = id, title = title,
               source = sub(paste0(here::here(), "/"), "", paste0(src, ".png")),
               status = if (length(got)) paste(got, collapse="+") else "MISSING")
  }))
}
man <- rbind(copy_set(MAIN, "main"), copy_set(ANNEX, "annex"))
write.csv(man, file.path(OUT, "tab_article_figures.csv"), row.names = FALSE)

cat(sprintf("\n=== Jeu de figures article (branche %s) -> %s ===\n", BRANCH, basename(OUT)))
print(man, row.names = FALSE)
miss <- man[man$status == "MISSING", ]
if (nrow(miss)) cat(sprintf("\n⚠ %d figure(s) MANQUANTE(S) :\n  %s\n",
                            nrow(miss), paste(sprintf("%s (%s)", miss$id, miss$title), collapse = "\n  ")))
