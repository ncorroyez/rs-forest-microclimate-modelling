# ==============================================================================
# Fig 4 (v3.2.3) + annex v3.2.0-vs-v3.2.3 validation scatter, station forcing.
# From tab_hobo_station_validation_stationTool.csv (obs + 4 sim columns).
# Main = v3.2.3 station; annex = v3.2.0 vs v3.2.3 (ranks better but under-buffers).
#   Rscript c1_validation_v323_scatter.R
# ==============================================================================
suppressPackageStartupMessages({ library(data.table); library(ggplot2); library(patchwork); source("scripts/_article_style.R") })
D <- fread("out_files/Chapter1/tables/tab_hobo_station_validation_stationTool.csv")
obs <- D$obs
mets <- function(s){ ok<-is.finite(s)&is.finite(obs); c(
  r=cor(s[ok],obs[ok]), RMSE=sqrt(mean((s[ok]-obs[ok])^2)), bias=mean(s[ok]-obs[ok]),
  MAE=mean(abs(s[ok]-obs[ok])), amp=100*mean(s[ok&obs<0])/mean(obs[ok&obs<0])) }

lim <- range(c(obs, D[["v3.2.0 station"]], D[["v3.2.3 station"]]), na.rm=TRUE); lim <- lim+c(-1,1)*diff(lim)*0.06
mk <- function(col, sub){
  s <- D[[col]]; m <- mets(s)
  lab <- sprintf("italic(r)==%.2f~~RMSE==%.2f~~bias==%+.2f~~'amp'==%.0f*'%%'", m["r"],m["RMSE"],m["bias"],m["amp"])
  ggplot(data.table(obs,s), aes(obs,s)) +
    geom_abline(slope=1,intercept=0,linetype="dashed",colour="grey55") +
    geom_hline(yintercept=0,linetype="dotted",colour="grey80")+geom_vline(xintercept=0,linetype="dotted",colour="grey80") +
    geom_point(size=2.6,alpha=0.8,colour="#2C7FB8") +
    annotate("text",x=-Inf,y=Inf,hjust=-0.05,vjust=1.4,size=4,parse=TRUE,label=lab) +
    coord_equal(xlim=lim,ylim=lim) +
    labs(x=expression("Observed "*Delta*T[max]~"(HOBO, "*degree*"C)"),
         y=expression("Simulated "*Delta*T[max]~"("*degree*"C)"), subtitle=sub) +
    theme_article(12) }

# MAIN: v3.2.3 station
p_main <- mk("v3.2.3 station", "v3.2.3 (station forcing) — reproduces buffering sign & ~26% of amplitude")
ggsave_article("out_files/Chapter1/figures/Fig_validation_v323_station", p_main, 5.5, 5.8)

# ANNEX: v3.2.0 vs v3.2.3 side by side
p_ann <- mk("v3.2.0 station","(a) v3.2.0 — ranks well but under-buffers (~9%, wrong mean sign)") +
         mk("v3.2.3 station","(b) v3.2.3 — buffers realistically at lower cross-plot r") +
  plot_annotation(title="Version comparison (station forcing): ranking vs buffering fidelity",
                  theme=theme(plot.title=element_text(size=12,face="bold")))
ggsave_article("out_files/Chapter1/figures/Fig_validation_v320_vs_v323", p_ann, 10, 5.6)

cat("=== metrics (station forcing) ===\n")
for(c in c("v3.2.0 station","v3.2.3 station")){ m<-mets(D[[c]]); cat(sprintf("  %-15s r=%.3f bias=%+.2f RMSE=%.2f MAE=%.2f ampl-captée=%.0f%%\n",c,m["r"],m["bias"],m["RMSE"],m["MAE"],m["amp"])) }
cat("DONE -> Fig_validation_v323_station + Fig_validation_v320_vs_v323\n")
