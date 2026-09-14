# ==============================================================================
# Réu 2026-06-26 #10 — Carte des 53 placettes HOBO par archétype structural P1–P4.
#   Rscript c1_archetype_map.R
# Out: FigStation_archetype_map.png
# ==============================================================================
suppressPackageStartupMessages({ library(sf); library(data.table); library(ggplot2); library(ggrepel) })
PAL <- c(P1="#D7191C", P2="#FDAE61", P3="#74C476", P4="#1A9850")  # thermal PAL_CLUSTER (article std)
LAB <- c(P1="P1 ouvert", P2="P2", P3="P3", P4="P4 dense")

g  <- st_read("in_files/data_Blois_utm31n.geojson", quiet=TRUE)
cl <- as.data.table(readRDS("outputs/figures_pipeline_z05/data/clusters.rds"))
d  <- merge(as.data.table(g)[, .(id_plot, x=coord_x_utm31n, y=coord_y_utm31n)], cl, by="id_plot")
d[, Cluster := factor(Cluster, levels=names(PAL))]
d[, num := sub("41_","",id_plot)]
cat(sprintf("placettes cartographiées : %d ; par archétype :\n", nrow(d))); print(d[, .N, by=Cluster][order(Cluster)])

lg <- sprintf("%s (n=%d)", LAB[levels(d$Cluster)], d[, .N, by=Cluster][order(Cluster)]$N)

p <- ggplot(d, aes(x, y, colour=Cluster)) +
  geom_point(size=3.2, alpha=0.9) +
  ggrepel::geom_text_repel(aes(label=num), size=2.4, colour="grey25", max.overlaps=20, seed=1) +
  scale_colour_manual(values=PAL, labels=lg, name="Archétype structural") +
  coord_sf(crs=st_crs(32631), datum=st_crs(32631)) +
  labs(x="Easting UTM 31N (m)", y="Northing UTM 31N (m)",
       title="Carte des 53 placettes HOBO par archétype structural (P1 ouvert → P4 dense)",
       subtitle="Clusters k-means gelés (LAI, Hmax, fCover) ; Blois, été 2021") +
  theme_bw(base_size=11) +
  theme(legend.position="right", plot.title=element_text(size=10.5,face="bold"),
        plot.subtitle=element_text(size=8,colour="grey35"), panel.grid.minor=element_blank())

ggsave("out_files/Chapter1/figures/FigStation_archetype_map.png", p, width=9, height=7, dpi=200, bg="white")
cat("DONE -> FigStation_archetype_map.png\n")
