# ==============================================================================
# Réu 2026-06-26 #4 — Figure conférence MEB (Microclimate Ecology & Biogeography).
# Mesure du microclimat : (a) schéma du dispositif (HOBO 1 m sous capuchon blanc
# en sous-bois vs station CHS 41 1,5 m en prairie) ; (b) série temporelle
# illustrative micro (HOBO, 1 plot ouvert P1 + 1 dense P4) vs macro (station) sur
# une vague de chaleur — montre le tamponnement mesuré.
#   Rscript c1_meb_figure.R
# Out: FigStation_MEB_field_setup.png
# ==============================================================================
suppressPackageStartupMessages({
  library(data.table); library(lubridate); library(ggplot2); library(patchwork)
  source("R/cluster_relabel.R")
})
PAL <- c(P1="#D7191C", P2="#FDAE61", P3="#74C476", P4="#1A9850")

# ---- (a) schéma du dispositif -------------------------------------------------
th  <- seq(0, 2*pi, length.out=80)
can <- data.frame(x=2.3 + 1.8*cos(th), y=16 + 7.5*sin(th))     # canopée elliptique
gnd <- data.frame(xmin=0, xmax=10, ymin=-0.6, ymax=0)
trunk <- data.frame(xmin=2.05, xmax=2.55, ymin=0, ymax=9)
grass <- data.frame(x=c(5.8,6.3,6.8,8.2,8.7,9.2), xend=c(5.8,6.3,6.8,8.2,8.7,9.2), y=0, yend=0.9)
pA <- ggplot() +
  geom_rect(data=gnd, aes(xmin=xmin,xmax=xmax,ymin=ymin,ymax=ymax), fill="#8C6D46") +
  geom_polygon(data=can, aes(x,y), fill="#2E7D32", alpha=0.5) +
  geom_rect(data=trunk, aes(xmin=xmin,xmax=xmax,ymin=ymin,ymax=ymax), fill="#5D4037") +
  geom_segment(data=grass, aes(x=x,xend=xend,y=y,yend=yend), colour="#7CB342", linewidth=0.7) +
  geom_hline(yintercept=c(1,1.5), linetype=3, colour="grey60") +
  annotate("text", x=9.9, y=1,   label="1 m",   hjust=1, vjust=1.5,  size=2.9, colour="grey35") +
  annotate("text", x=9.9, y=1.5, label="1,5 m", hjust=1, vjust=-0.6, size=2.9, colour="grey35") +
  # HOBO (sous-bois, 1 m) — étiquette dans l'espace clair du sous-bois
  annotate("point", x=2.3, y=1, size=3.4, colour="#D7191C") +
  annotate("point", x=2.3, y=1, size=1.6, colour="white", shape=6) +
  annotate("segment", x=2.3, xend=0.45, y=1, yend=4.6, colour="#D7191C", linewidth=0.3) +
  annotate("text", x=0.4, y=5.4, label="HOBO 1 m\n(capuchon blanc\nanti-radiation)",
           hjust=0, size=3, colour="#D7191C", fontface="bold", lineheight=0.9) +
  # station (prairie, 1,5 m)
  annotate("segment", x=7.5, xend=7.5, y=0, yend=1.5, colour="grey30", linewidth=0.9) +
  annotate("point", x=7.5, y=1.5, size=3.4, colour="grey15") +
  annotate("text", x=7.5, y=4.4, label="Station CHS 41\n1,5 m (prairie)",
           size=3, colour="grey20", fontface="bold", lineheight=0.9) +
  annotate("text", x=2.3, y=24.8, label="canopée", size=3.1, colour="#1B5E20") +
  annotate("text", x=2.3, y=-1.4, label="sous-bois (placette)", size=3, colour="#5D4037") +
  annotate("text", x=7.5, y=-1.4, label="ouvert (référence macro)", size=3, colour="grey25") +
  coord_cartesian(xlim=c(0,10), ylim=c(-2, 27)) +
  labs(y="Hauteur (m)", title="(a) Dispositif de mesure du microclimat") +
  theme_void(base_size=11) +
  theme(axis.title.y=element_text(angle=90, size=9), axis.text.y=element_text(size=8),
        axis.line.y=element_line(colour="grey60"), axis.ticks.y=element_line(colour="grey60"),
        plot.title=element_text(size=10.5, face="bold"))

# ---- (b) série temporelle micro vs macro sur une vague de chaleur -------------
cl <- fread("out_files/Chapter1/tables/tab_hobo_perplot_cluster.csv")
p1id <- cl[P=="P1"][which.max(dTmax_obs), id_plot]  # P1 le plus ouvert/amplifiant (exemple contrasté)
p4id <- cl[P=="P4"][which.min(dTmax_obs), id_plot]  # P4 le plus tamponnant (exemple contrasté)
hobo <- fread("in_files/Blois_data_temperature.csv")
hobo <- hobo[position_sensor=="a" & id_plot %in% c(p1id,p4id)]
hobo[, time := as.POSIXct(datetime, format="%Y-%m-%d %H:%M:%S", tz="UTC")]
hobo <- hobo[, .(t=mean(t_hobo,na.rm=TRUE)), by=.(id_plot, time=floor_date(time,"hour"))]

raw <- fread("MetHor2021.txt", sep=";", header=TRUE, encoding="Latin-1"); setnames(raw,1,"code")
tcol <- grep("instantan", names(raw), value=TRUE)[1]
sta <- raw[code=="CHS 41", .(Date, H=`Heure (TU)`, t=as.numeric(get(tcol)))]
sta[, `:=`(time=as.POSIXct(Date,format="%d/%m/%Y",tz="UTC")+lubridate::hours(H%/%100), id_plot="Station 1,5 m")]
sta <- sta[is.finite(t), .(id_plot, time, t)]

# hot window: 5 jours autour du pic de la station
peak <- sta[which.max(t), as.Date(time)]
win <- seq(peak-2, peak+2, by="day")
D <- rbind(hobo, sta)[as.Date(time) %in% win]
D[, grp := fcase(id_plot==p1id, sprintf("P1 ouvert (%s)", p1id),
                 id_plot==p4id, sprintf("P4 dense (%s)", p4id),
                 default="Station 1,5 m (macro)")]
COLS <- c("#D7191C","#1A9850","grey25"); names(COLS) <- c(sprintf("P1 ouvert (%s)",p1id), sprintf("P4 dense (%s)",p4id), "Station 1,5 m (macro)")
pB <- ggplot(D, aes(time, t, colour=grp)) +
  geom_line(linewidth=0.6) +
  scale_colour_manual(values=COLS, name=NULL) +
  labs(x=NULL, y="Température de l'air (°C)",
       title="(b) Tamponnement mesuré pendant une vague de chaleur",
       subtitle=sprintf("%s au %s — exemples contrastés (placette la plus ouverte vs la plus dense) : le sous-bois dense écrête les pics diurnes",
                        format(min(win),"%d/%m"), format(max(win),"%d/%m/%Y"))) +
  theme_bw(base_size=11) + theme(legend.position="bottom",
       plot.title=element_text(size=10.5,face="bold"), plot.subtitle=element_text(size=8,colour="grey35"))

ggsave("out_files/Chapter1/figures/FigStation_MEB_field_setup.png",
       pA + pB + plot_layout(widths=c(0.8,1.2)), width=12, height=5, dpi=200, bg="white")
cat(sprintf("plots: ouvert=%s dense=%s ; fenêtre %s..%s\n", p1id, p4id, min(win), max(win)))
cat("DONE -> FigStation_MEB_field_setup.png\n")
