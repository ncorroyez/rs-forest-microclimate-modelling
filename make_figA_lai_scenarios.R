# Annex figure: LAI fed to MuSICA per scenario, per HOBO plot, over the summer
# period (static => flat lines). One-sided LAI. Style from NC_Full skill.
suppressPackageStartupMessages({ library(data.table); library(ggplot2); library(rmusica); library(musica.tools) })
src <- list.files("R", pattern="\\.R$", full.names=TRUE); src <- src[!grepl("/(h1_|lovb_)", src)]; invisible(lapply(src, source)); source("Chapter3_config.R")
source("/home/corroyez/Documents/NC_Full/scripts/_article_style.R")
prep <- load_lai_prep(CFG_C3); df <- as.data.table(prep$df_plots); dopt <- CFG_C3$d_opt_m
# fills (same as scenario script)
df[is.na(LAI_ALS_DOPT), LAI_ALS_DOPT := LAI_ALS]
ropt <- mean(df$LAI_S2_DOPT,na.rm=TRUE)/mean(df$LAI_S2_ATBD,na.rm=TRUE)
df[is.na(LAI_S2_DOPT), LAI_S2_DOPT := LAI_S2_ATBD*ropt]
rA <- mean(df$LAI_ALS)/mean(df$LAI_S2_ATBD); bT <- mean(df$LAI_ALS_DOPT)/mean(df$LAI_S2_ATBD)
S <- df[, .(id_plot,
  `S2 ATBD (raw)`=LAI_S2_ATBD, `S2 opt (raw)`=LAI_S2_DOPT, `LiDAR full`=LAI_ALS,
  `LiDAR d_opt`=LAI_ALS_DOPT, `S2 ATBD x ratio->full`=LAI_S2_ATBD*rA,
  `S2 ATBD x ratio->d_opt`=LAI_S2_ATBD*bT,
  `fusion (S2/LiDAR by height)`=ifelse(Hmax<dopt, LAI_S2_DOPT, LAI_ALS))]
fwrite(S, "/home/corroyez/Documents/NC_Full/manuscripts/ch3/tables/Table4g_scenario_lai_perplot.csv")
m <- melt(S, id.vars="id_plot", variable.name="scenario", value.name="LAI")
ord <- m[, .(mu=mean(LAI)), by=scenario][order(-mu), scenario]
m[, scenario := factor(scenario, levels=ord)]
# static over summer: replicate flat over day-of-year 152..273
doy <- 152:273
mt <- m[, .(doy=doy, LAI=LAI), by=.(id_plot, scenario)]
mu <- m[, .(mu=mean(LAI)), by=scenario]
g <- ggplot(mt, aes(doy, LAI, group=id_plot)) +
  geom_line(alpha=0.18, colour="#1A9850", linewidth=0.3) +
  geom_hline(data=mu, aes(yintercept=mu), colour="#D7191C", linewidth=0.7) +
  geom_text(data=mu, aes(x=152, y=mu, label=sprintf("mean %.1f", mu)), inherit.aes=FALSE,
            hjust=0, vjust=-0.4, size=2.6, colour="#D7191C") +
  facet_wrap(~scenario, ncol=4) +
  labs(x="Day of year (summer 2021)", y="LAI (one-sided, m² m⁻²)") +
  theme_article(11)
outdir <- "/home/corroyez/Documents/NC_Full/outputs/figures_article_ch3"
ggsave_article(file.path(outdir,"FigA_scenario_lai_perplot"), g, 9, 4.6)
ggsave_article("/home/corroyez/Documents/NC_Full/manuscripts/ch3/figures/FigA_scenario_lai_perplot", g, 9, 4.6)
cat(sprintf("DONE. ratios rA=%.3f bT=%.3f ropt=%.3f | scenarios mean LAI:\n", rA, bT, ropt)); print(mu[order(-mu)])
