---
title: "Correction du vent ERA5 → sommet de canopée"
subtitle: "Note pour Rémi — MuSICA Chapitre 1 (Blois, 2021)"
date: 2026-07-10
---

# L'idée en deux phrases

ERA5 donne le vent à **10 m au-dessus d'un sol ouvert**. MuSICA a besoin du vent **juste au-dessus de chaque
canopée** (h + 2 m). On corrige avec un profil de vent logarithmique : au-dessus d'une canopée de 15 à 35 m, le
vent réel n'est que ~0,4 fois la valeur ERA5 (donc le vent brut était ~2,3 fois trop fort).

# La formule

$$
U(h+2) = U(10) \cdot \frac{\ln(2 + h - d) - \ln(z_0)}{\ln(10) - \ln(z_{0,\text{ERA}})}
$$

- $d = 0{,}7\,h$ : hauteur de déplacement
- $z_0 = 0{,}1\,h$ : rugosité de la canopée
- $z_{0,\text{ERA}} = 0{,}44$ m : rugosité ERA5 de la maille de Blois

Le facteur ne dépend que de la hauteur $h$ du plot et vaut ~0,41–0,47 pour $h$ = 15 à 35 m.

![Vent réel au-dessus de la canopée rapporté au vent ERA5 à 10 m : avec la rugosité ERA5 de Blois (0,44 m), le vent au sommet du couvert vaut ~0,4 fois la valeur à 10 m entre 15 et 35 m, soit un forçage brut ~2,3 fois trop fort.](figures/article_v323/G1_wind_profile_correction.png){width=70%}

# Le code (R)

```r
library(ncdf4)
Z0_ERA_BLOIS <- 0.44   # rugosité ERA5 (fsr), maille de Blois

# facteur multiplicatif U(h+2m)/U(10m) pour une canopée de hauteur h (m)
wind_factor <- function(h, z0ERA = Z0_ERA_BLOIS) {
  d <- 0.7 * h; z0 <- 0.1 * h
  (log(2 + h - d) - log(z0)) / (log(10) - log(z0ERA))
}

# renvoie un forçage (chemin RELATIF au projet) dont le vent est multiplié par wind_factor(hmax).
# Cache par facteur arrondi -> ~40 fichiers couvrent tous les plots.
windcorr_forcing <- function(hmax, forc_base = "in_files/FR-Blo_2021_v2.nc",
                             cache_dir = "out_files/windcorr_forc", z0ERA = Z0_ERA_BLOIS) {
  dir.create(cache_dir, recursive = TRUE, showWarnings = FALSE)
  f  <- wind_factor(hmax, z0ERA)
  ff <- file.path(cache_dir, sprintf("forc_f%.2f.nc", round(f, 2)))
  if (!file.exists(ff) || file.size(ff) < 1e5) {
    file.copy(forc_base, ff, overwrite = TRUE)
    nc <- nc_open(ff, write = TRUE)
    for (v in c("Wind_E", "Wind_N")) ncvar_put(nc, v, ncvar_get(nc, v) * f)  # même direction, norme réduite
    nc_close(nc)
  }
  ff
}
```

Usage, avant chaque run (à la place du forçage brut) :

```r
source("R/wind_correction.R")
forcing <- windcorr_forcing(plot$Hmax)
```

**Deux points** : (1) on multiplie les deux composantes `Wind_E`/`Wind_N`, donc on réduit la norme du vent sans
changer sa direction ; (2) le chemin renvoyé doit être **relatif au projet** (`setup_musica` monte le forçage
via `ln -s ../<forcing>`, un `/tmp/...` absolu casserait le lien).

# Obtenir z0,ERA (une fois, par site)

C'est le champ ERA5 **`forecast_surface_roughness`**, moyenné sur la maille du site et la période d'étude.
Récupération via l'API Copernicus CDS (`cdsapi`) :

```python
import cdsapi, xarray as xr

cdsapi.Client().retrieve(
    "reanalysis-era5-single-levels",
    {
        "product_type": "reanalysis",
        "variable": "forecast_surface_roughness",
        "year": "2021",
        "month": ["06", "07", "08", "09"],
        "day":  [f"{d:02d}" for d in range(1, 32)],
        "time": [f"{h:02d}:00" for h in range(24)],
        "area": [47.75, 1.0, 47.25, 1.5],   # [Nord, Ouest, Sud, Est] autour de Blois
        "format": "netcdf",
    },
    "z0_blois_2021.nc",
)

ds  = xr.open_dataset("z0_blois_2021.nc")
z0  = float(ds["fsr"].sel(latitude=47.5, longitude=1.25, method="nearest").mean())
print(f"z0_ERA (Blois, juin-sept 2021) = {z0:.2f} m")   # -> 0.44
```

Pour **un autre site**, il suffit de changer `area` / `latitude` / `longitude` et la période : le reste du code
est identique.

# Effet sur les résultats (Chapitre 1)

Décision (courriel Jérôme, 8 juillet) : **intégration complète**, le forçage corrigé devient la baseline.
L'ordre d'attribution est **inchangé** (quantité de feuilles domine, profil vertical dernier). Sur la validation
HOBO, la correction dégrade légèrement la corrélation de ΔTmax (r 0,60 → 0,51) et augmente le biais chaud
(+0,56 → +0,65 °C) : elle est adoptée pour la **justesse physique**, pas parce qu'elle améliore l'ajustement.

![Effet de la correction sur les 53 loggers : le taux de tampon (buffering slope) est quasi inchangé tandis que l'amplitude de ΔTmax se dégrade un peu (biais chaud de +0,56 à +0,65 °C), l'effet croissant avec la hauteur de canopée.](figures/article_v323/windcorr_compare.png){width=85%}
