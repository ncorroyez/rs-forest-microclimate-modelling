######################################@
## Define call function
# devtools::load_all("../../musica.tools/musica.tools/")
library(musica.tools)
library(rmusica)
library(ncdf4)
library(magrittr)
library(dplyr)
library(tidyr)
library(ggplot2)

out1 <- nc_open("out_files/musica_out_FR-Bil_2022_run1.nc")
out2 <- nc_open("out_files/musica_out_FR-Bil_2022_run2.nc")


df_Tair <- get_variable(out1, "Tair_z")
names(attributes(df_Tair))
attr(df_Tair, "relative_height")


ggplot_variable(filter(df_Tair, nair == 5), out.type = "standard")

filter(df_Tair, nair == 5) %>%
  select(-nair) %>%
ggplot_variable(out.type = "daily_heatmap")

ggplot_variable(df_Tair, out.type = "heatmap", layer.y = TRUE) +
  ggtitle("Mon Titre")
# +
#   scale_y_discrete("Mon Axe Y")


out.all <- list("run1" = out1,
                "run2" = out2)

shinymusica(out.all)
