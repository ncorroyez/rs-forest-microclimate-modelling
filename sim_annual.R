suppressPackageStartupMessages({library(terra);library(sf);library(ncdf4);library(lubridate);library(dplyr);library(tidyr);library(stringr);library(purrr);library(rmusica);library(musica.tools)})
src<-list.files("R",pattern="\\.R$",full.names=TRUE);src<-src[!grepl("/(h1_|lovb_)",src)];invisible(lapply(src,source));source("Chapter3_config.R")
prep<-load_lai_prep(CFG_C3)
sc<-make_all_scenarios_c3(prep$ts_by_plot,list_year=CFG_C3$list_year,d_opt=CFG_C3$d_opt_m,mode=CFG_C3$scenarios_mode)["DYN_S2_ANNUAL"]
stopifnot(!is.null(sc$DYN_S2_ANNUAL)); cat("simulating DYN_S2_ANNUAL\n")
ds<-CFG_C3$date_seq;dm<-extract_macro_daily(CFG_C3$forcing_file,ds);hd<-read_hobo_daily(CFG_C3$hobo_temp_csv,ds,dm,CFG_C3$ids_to_remove)
validate_scenarios_at_hobos(prep$df_plots,hd,sc,file.path(CFG_C3$out_dir,"nc"),dm,ds,CFG_C3$forcing_file,CFG_C3$musica_cmd,FALSE)
cat(sprintf("DYN_S2_ANNUAL nc: %d\nSIM_DONE\n",length(list.files(file.path(CFG_C3$out_dir,"nc","DYN_S2_ANNUAL"),pattern="\\.nc$"))))
