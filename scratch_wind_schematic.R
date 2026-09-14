# Pedagogical schematic: why U(h+2m)/U(10m) < 1 over tall canopy.
# Two neutral log wind profiles sharing the same u*/k: open field (ERA5, small z0)
# vs forest (d=0.7h, z0=0.1h). ERA5 reads its wind at 10 m over the smooth field;
# MuSICA needs it at h+2 m over the rough forest — much lower on the forest profile.
suppressPackageStartupMessages({library(ggplot2);library(ggtext)})
h <- 25; d <- 0.7*h; z0f <- 0.1*h; z0e <- 0.1              # example canopy 25 m
Uf <- function(z) ifelse(z > d + z0f, log((z - d)/z0f), NA) # forest profile
Ue <- function(z) ifelse(z > z0e,      log(z/z0e),       NA) # open-field profile
zz <- seq(0.2, 35, 0.1)
df <- rbind(data.frame(z=zz, U=Uf(zz), prof="Forest (real canopy)"),
            data.frame(z=zz, U=Ue(zz), prof="Open field (ERA5)"))
u_h2 <- Uf(h+2); u_10 <- Ue(10)
PAL <- c("Forest (real canopy)"="#1A9850", "Open field (ERA5)"="#2c7fb8")
p <- ggplot(df, aes(U, z, colour=prof)) +
  annotate("rect", xmin=-Inf, xmax=Inf, ymin=0, ymax=h, fill="#8fbf6f", alpha=0.16) +   # canopy volume
  annotate("segment", x=-Inf, xend=Inf, y=0, yend=0, colour="#6b4a2a", linewidth=1.1) +   # ground
  geom_hline(yintercept=d,   linetype="dashed", colour="grey55", linewidth=0.3) +
  geom_hline(yintercept=h,   linetype="dashed", colour="#1A9850", linewidth=0.3) +
  geom_hline(yintercept=10,  linetype="dashed", colour="#2c7fb8", linewidth=0.3) +
  geom_line(linewidth=1.1, na.rm=TRUE) +
  # the two sampling points + their wind values
  annotate("point", x=u_h2, y=h+2, colour="#1A9850", size=3) +
  annotate("point", x=u_10, y=10,  colour="#2c7fb8", size=3) +
  annotate("segment", x=u_h2, xend=u_10, y=h+2, yend=h+2, linewidth=0.4, colour="grey30",
           arrow=arrow(length=unit(2,"mm"), ends="both")) +
  scale_colour_manual(values=PAL, name=NULL) +
  scale_y_continuous(breaks=c(0,10,d,h,h+2), labels=c("0","10","d = 17.5","h = 25","h+2")) +
  labs(x="wind speed  U  (arbitrary, ∝ log height)", y="height z (m)") +
  coord_cartesian(xlim=c(0, 5.2), ylim=c(0,35), clip="off") +
  theme_bw(13) + theme(panel.grid=element_blank(), legend.position=c(0.98,0.02),
        legend.justification=c(1,0), legend.background=element_rect(fill="white",colour="grey80")) +
  annotate("richtext", x=u_10, y=10.6, label="**ERA5 reads U here**<br>(10 m, smooth field)", hjust=1.08, vjust=0,
           size=3.2, fill=NA, label.color=NA, colour="#2c7fb8") +
  annotate("richtext", x=u_h2, y=h+2, label="**MuSICA needs U here**<br>(h+2 m, above forest)", hjust=1.06, vjust=0.3,
           size=3.2, fill=NA, label.color=NA, colour="#1A9850") +
  annotate("text", x=(u_h2+u_10)/2, y=h+4, label=sprintf("factor = %.2f", u_h2/u_10), size=3.6, colour="grey20") +
  annotate("richtext", x=0.15, y=d-1.4, label="*z*<sub>0</sub> = 0.1h just above d", hjust=0, vjust=1,
           size=2.9, fill=NA, label.color=NA, colour="grey45") +
  annotate("text", x=0.5, y=h-2, label="canopy", colour="#3f6b2a", size=3.6, fontface=2, hjust=0)
ggsave("/tmp/claude-1001/wind_profile_schematic.png", p, width=8, height=6, dpi=200, bg="white")
cat(sprintf("h=%d: U(h+2)=%.2f  U(10)=%.2f  ratio=%.2f\nDONE\n", h, u_h2, u_10, u_h2/u_10))
