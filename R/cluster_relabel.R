# ==============================================================================
# Cluster relabel — single source of truth for cluster display labels.
#
# Underlying cLHS data uses Cluster IDs 1, 2, 3, 4.
# We DISPLAY them re-ordered by ASCENDING LAI, so :
#   old Cluster 3 (LAI ≈ 2.4, sparsest) → displayed as "P1"
#   old Cluster 2 (LAI ≈ 6.0)            → displayed as "P2"
#   old Cluster 4 (LAI ≈ 7.7)            → displayed as "P3"
#   old Cluster 1 (LAI ≈ 10.8, densest)  → displayed as "P4"
#
# The palette is then attached to the NEW labels :
#   P1 (smallest LAI)  = "#1B9E77" (green)
#   P2                  = "#D95F02" (orange)
#   P3                  = "#7570B3" (violet)
#   P4 (largest LAI)    = "#E7298A" (pink)
# ==============================================================================

# Mapping : original Cluster ID (character) → display label
CLUSTER_RELABEL <- c("3" = "P1",
                      "2" = "P2",
                      "4" = "P3",
                      "1" = "P4")

# Unified palette keyed by the NEW labels — thermal red→green (P1 sparse/warm →
# P4 dense/cool), aligned with scripts/_article_style.R and the amplify/buffer
# convention (single source of truth for cluster colours across the article).
PAL_CLUSTER <- c(P1 = "#D7191C",
                  P2 = "#FDAE61",
                  P3 = "#74C476",
                  P4 = "#1A9850")

PAL_SHAPES  <- c(P1 = 16, P2 = 17, P3 = 15, P4 = 18)

#' Apply the cluster remap to a vector / column.
#'
#' @param x character/integer/factor of old Cluster IDs (1..4).
#' @return character vector of new labels ("P1".."P4"), with factor levels
#'   ordered P1..P4 (left to right in plots).
relabel_cluster <- function(x) {
  out <- CLUSTER_RELABEL[as.character(x)]
  factor(out, levels = paste0("P", 1:4))
}
