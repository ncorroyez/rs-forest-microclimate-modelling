suppressPackageStartupMessages({
  library(png)
  library(grid)
})

out_pdf <- "outputs/figures/Chap1_figures.pdf"

main_figs <- c(
  "outputs/figures/01_h1_forward_curve.png",
  "outputs/figures/02_h1_ranking.png",
  "outputs/figures/03_h1_loo_scatters_2x2.png",
  "outputs/figures/04_h1_scatters_vs_ref.png",
  "outputs/figures/05_h2_paired_real_vs_uniform.png",
  "outputs/figures/06_h2_vertical_delta_by_cluster.png",
  "outputs/figures/07_h2_vertical_3way_canicule.png",
  "outputs/figures/08_clusters_lad_profiles.png",
  "outputs/figures/09_fpca_harmonics.png",
  "outputs/figures/10_fpca_reconstruction.png"
)

annex_figs <- sort(list.files("outputs/figures/annex", pattern = "\\.png$",
                               full.names = TRUE))

all_figs <- c(main_figs, annex_figs)
all_figs <- all_figs[file.exists(all_figs)]

cat(sprintf("Assemblage de %d figures → %s\n", length(all_figs), out_pdf))

cairo_pdf(out_pdf, width = 14, height = 9, onefile = TRUE)
for (f in all_figs) {
  img <- readPNG(f)
  grid.newpage()
  grid.raster(img, width = unit(1, "npc"), height = unit(1, "npc"))
  grid.text(basename(f), x = 0.01, y = 0.01, just = c("left", "bottom"),
            gp = gpar(fontsize = 7, col = "grey50"))
}
dev.off()

cat(sprintf("PDF généré : %s\n", out_pdf))
