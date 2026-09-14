# ==============================================================================
# Figure 2 — analysis pipeline schematic, rendered from a graphviz source.
#
# Needs the graphviz `dot` binary on PATH; the script stops if it is absent.
#
# Reads : scripts/c1_fig2_methodo.dot
# Writes: out_files/Chapter1/figures/Fig2_methodo.png   (-Gdpi=200)
#         out_files/Chapter1/figures/Fig2_methodo.pdf
#   Rscript scripts/c1_fig2_methodo.R
# ==============================================================================
# Figure 2 — analysis pipeline schematic. Renders scripts/c1_fig2_methodo.dot with
# graphviz. This was the last figure in the chapter with no generating script, and it
# had drifted from the analysis (it showed "10 m" traits and a "±1 SD perturbation",
# both superseded). Source of truth is now the .dot file next to this script.
DOT <- "scripts/c1_fig2_methodo.dot"
OUT <- "out_files/Chapter1/figures/Fig2_methodo"
stopifnot(file.exists(DOT))
if (nzchar(Sys.which("dot")) == FALSE) stop("graphviz `dot` not found on PATH")
dir.create(dirname(OUT), recursive = TRUE, showWarnings = FALSE)
for (fmt in c("png", "pdf")) {
  args <- c(sprintf("-T%s", fmt), if (fmt == "png") "-Gdpi=200", DOT, "-o", paste0(OUT, ".", fmt))
  st <- system2("dot", args)
  if (st != 0) stop("dot failed for format ", fmt)
}
cat(sprintf("DONE -> %s.png (%.0f kB) + .pdf\n", OUT, file.size(paste0(OUT, ".png"))/1024))
