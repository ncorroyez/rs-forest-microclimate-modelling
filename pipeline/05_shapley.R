# ==============================================================================
# PIPELINE STAGE 05 — Exact Shapley attribution per plot (ΔTmax & ΔVPDmax, 1 m).
# Weight = factorial(s)·factorial(4-s-1)/factorial(4). Additivity: Σφ = REF − Base.
# Regression assert (constraint #6): median Σφ on ΔTmax(all) ≈ -0.311 °C.
# Output: tab_shapley_perplot.csv (+ rds)
# ==============================================================================

if (!exists("PIPE")) source(here::here("pipeline/00_config.R"))
cli_h1("STAGE 05 — exact Shapley")

coal <- as.data.table(readRDS(file.path(PIPE$OUT_DATA, "coal_metrics.rds")))
clu  <- as.data.table(readRDS(file.path(PIPE$OUT_DATA, "clusters.rds")))
bits <- PIPE$BITS

shap_one <- function(vec) {                         # vec named by bit
  v0 <- as.numeric(vec[PIPE$BASE_BIT]); v <- function(b) as.numeric(vec[b]) - v0
  mb <- function(bo,n=4L){ b<-rep("0",n); b[bo]<-"1"; paste(b,collapse="") }
  out <- numeric(4L)
  for (vb in 1:4){ others<-setdiff(1:4,vb); tot<-0
    for (s in 0:3){ subs <- if(s==0) list(integer(0)) else combn(others,s,simplify=FALSE)
      w <- factorial(s)*factorial(4-s-1)/factorial(4)
      for (S in subs) tot <- tot + w*(v(mb(c(S,vb)))-v(mb(S))) }
    out[vb] <- tot }
  setNames(out, c("LAI","Hmax","fCover","LAD"))
}
phi_for <- function(col, period_lbl) {
  W <- dcast(coal, id_plot ~ bit, value.var = col)
  rbindlist(lapply(seq_len(nrow(W)), function(i){
    vec <- as.numeric(unlist(W[i, ..bits])); names(vec) <- bits
    if (any(is.na(vec))) return(NULL)
    p <- shap_one(vec)
    data.table(id_plot=W$id_plot[i], trait=names(p), phi=unname(p), period=period_lbl) }))
}
phi <- rbindlist(list(
  cbind(metric="Tmax", phi_for("Tmax_all","All period")),
  cbind(metric="Tmax", phi_for("Tmax_hot","10% hottest")),
  cbind(metric="VPD",  phi_for("VPD_all","All period")),
  cbind(metric="VPD",  phi_for("VPD_hot","10% hottest"))))
phi <- merge(phi, clu, by="id_plot")
phi[, period := factor(period, levels=c("All period","10% hottest"))]
saveRDS(phi, file.path(PIPE$OUT_DATA, "shapley_perplot.rds"))
fwrite(phi, file.path(PIPE$OUT_TAB, "tab_shapley_perplot.csv"))

cli_h2("Median φ per metric × period × trait")
print(phi[, .(median=round(median(phi,na.rm=TRUE),3)), by=.(metric,period,trait)][order(metric,period,median)])

# ---- additivity + regression asserts ----------------------------------------
W <- dcast(coal, id_plot ~ bit, value.var = "Tmax_all")[, .(id_plot, tot = get(PIPE$REF_BIT) - get(PIPE$BASE_BIT))]
sphi <- phi[metric=="Tmax" & period=="All period", .(s=sum(phi)), by=id_plot]
m <- merge(sphi, W, by="id_plot")
pipe_assert(all(pipe_near(m$s, m$tot, 1e-6)), "Additivity per plot: Σφ == ΔTmax(REF)-ΔTmax(Base)")  # always (math identity)
if (PIPE$strict) {
  pipe_assert(pipe_near(median(m$s), PIPE$ASSERT$sum_phi_tmax, 0.02),
              sprintf("Regression: median Σφ(ΔTmax,all) = %.3f (target %.3f)", median(m$s), PIPE$ASSERT$sum_phi_tmax))
} else cli_alert_info("[branch {PIPE$BRANCH}] median Σφ(ΔTmax,all) = {round(median(m$s),3)} (legacy target {PIPE$ASSERT$sum_phi_tmax}, not enforced)")
cli_alert_success("Stage 05 done.")
