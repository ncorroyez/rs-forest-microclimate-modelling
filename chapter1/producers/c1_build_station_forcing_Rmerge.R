# ==============================================================================
# Station forcing via LOCAL R MERGE (robust; bypasses prep_site_forcing tool bugs).
# Base = musica_in_Blois.nc (ERA5 secondaries, already on MuSICA schema), subset
# to 2021; GRAFT the station CHS41 Tair + Qair(from RH) + Rainf(precip). Everything
# else (SWdown, LWdown, PSurf, wind, CO2) stays ERA5 — exactly what CHS41 lacks.
# Block-A coherence: the above-canopy T driving MuSICA = station T = ΔTmax ref.
# Nothing overwritten → out_files/MuSICA_in_CHS41-Blois_2021-station_Rmerge.nc
#   Rscript c1_build_station_forcing_Rmerge.R
# ==============================================================================
suppressPackageStartupMessages({ library(ncdf4); library(data.table); library(lubridate) })
BASE <- "in_files/musica_in_Blois.nc"
STA  <- "prep_site_forcing/MetHor2021_CHS41.csv"
DST  <- "out_files/MuSICA_in_CHS41-Blois_2021-station_Rmerge.nc"

# ---- base forcing time axis, subset to 2021 ----------------------------------
nc <- nc_open(BASE)
tu <- ncatt_get(nc,"time","units")$value; th <- ncvar_get(nc,"time")
t0 <- as.POSIXct(sub(".*since ","",tu),tz="UTC"); tv <- t0 + th*3600
keep <- which(format(tv,"%Y")=="2021")
cat(sprintf("base 2021 steps: %d (%s -> %s, %s-centred)\n", length(keep),
            format(tv[keep][1],"%Y-%m-%d %H:%M"), format(tv[keep][length(keep)],"%H:%M"),
            format(tv[keep][1],"%M")))
hr <- floor_date(tv[keep],"hour")                       # hour key for join

# ---- station: read + hourly key ----------------------------------------------
sta <- fread(STA, sep=";")
sta[, dt := as.POSIXct(datetime, format="%d.%m.%Y %H:%M", tz="UTC")]
sta[, hr := floor_date(dt,"hour")]
sta[Tair <= -900, Tair := NA][RH <= -900, RH := NA][precip <= -900, precip := NA]
setkey(sta, hr)
S <- sta[.(hr)]                                          # aligned to base 2021 steps
cat(sprintf("station matched: Tair %d, RH %d, precip %d / %d steps\n",
            sum(is.finite(S$Tair)), sum(is.finite(S$RH)), sum(is.finite(S$precip)), length(keep)))

# ---- read base vars (2021 subset), graft station ------------------------------
getv <- function(v){ a <- ncvar_get(nc,v); if (length(dim(a))>=1) { d<-dim(a); tdim<-which(d==length(th));
  if(length(tdim)==1){ idx<-lapply(seq_along(d),function(i) if(i==tdim) keep else seq_len(d[i])); do.call(`[`, c(list(a),idx,list(drop=FALSE))) } else a } else a }
Tair <- getv("Tair"); Qair <- getv("Qair"); PSurf <- getv("PSurf"); Rainf <- getv("Rainf")
dtemplate <- dim(Tair)                                   # e.g. [z? , y, x, time] or [time,...]
# flatten helper: base vars here are (time,z,y,x) or (time,y,x); station is 1 per time
esat <- function(Tc) 610.94*exp(17.625*Tc/(Tc+243.04))   # Pa (Alduchov-Eskridge)

Tk_sta <- S$Tair + 273.15                                # station Tair K
# Qair (kg/kg) from station RH + station T + base PSurf. PSurf per-time vector:
psv <- as.numeric(PSurf); if (length(psv)!=length(keep)) psv <- rep(median(psv,na.rm=TRUE), length(keep))
e   <- (S$RH/100) * esat(S$Tair)                         # Pa
q_sta <- 0.622*e/(psv - 0.378*e)                         # specific humidity kg/kg
rain_sta <- S$precip/3600                                # mm/h -> kg/m2/s (base Rainf units)

# replace along the time dimension (base vars are (time,z,y,x)/(time,y,x): time is dim 1 here? check)
# musica_in_Blois: Tair(time,z,y,x). After 2021 subset, time is last-collapsed dim we kept.
# Simplest robust path: operate on the raw arrays with time index = `keep`.
put_time <- function(arr, vec){                          # arr dims incl a time axis of length(keep); set each time slice to vec[t]
  d <- dim(arr); td <- which(d==length(keep))[1]
  idx <- slice.index(arr, td); out <- arr
  for(t in seq_len(length(keep))) out[idx==t] <- ifelse(is.finite(vec[t]), vec[t], arr[idx==t][1])
  out
}
Tair2  <- put_time(Tair,  Tk_sta)
Qair2  <- put_time(Qair,  q_sta)
Rainf2 <- put_time(Rainf, rain_sta)

# ---- write new nc = clone base(2021) with Tair/Qair/Rainf replaced ------------
newdims <- list()
tdimv <- ncdim_def("time", tu, th[keep], unlim=TRUE)
# rebuild dims from base, swapping time
dimmap <- lapply(nc$dim, function(dd) if(dd$name=="time") tdimv else ncdim_def(dd$name, dd$units, dd$vals, unlim=dd$unlim))
vars <- lapply(nc$var, function(v){
  dims <- lapply(v$dim, function(dd) dimmap[[dd$name]])
  ncvar_def(v$name, v$units, dims, missval=if(is.null(v$missval)) NA else v$missval,
            prec=if(v$prec=="double") "double" else "float")
})
out <- nc_create(DST, vars)
for (v in nc$var){
  val <- if (v$name=="Tair") Tair2 else if (v$name=="Qair") Qair2 else if (v$name=="Rainf") Rainf2 else getv(v$name)
  ncvar_put(out, v$name, val)
}
nc_close(out); nc_close(nc)
cat(sprintf("WROTE %s\n", DST))

# ---- sanity ------------------------------------------------------------------
c2 <- nc_open(DST)
Ta <- ncvar_get(c2,"Tair"); Qa <- ncvar_get(c2,"Qair"); cat(sprintf(
  "check: Tair K [%.1f,%.1f] mean %.1f | Qair [%.4f,%.4f] | n_time=%d\n",
  min(Ta,na.rm=T),max(Ta,na.rm=T),mean(Ta,na.rm=T),min(Qa,na.rm=T),max(Qa,na.rm=T), c2$dim$time$len))
nc_close(c2)
