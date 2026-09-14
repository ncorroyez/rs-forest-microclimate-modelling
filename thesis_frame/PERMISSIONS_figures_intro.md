# Figures empruntées — crédits et formalités

L'introduction reproduit cinq figures qui ne sont pas dessinées par un script de `review/`.
Chacune est créditée en légende. Les autorisations formelles sont une formalité à régler
avant le dépôt, pas un préalable à l'écriture.

**Numérotation.** Les numéros ci-dessous sont ceux du manuscrit actuel, après l'insertion de
la Figure 1.4 (bilan d'énergie) le 11/09/2026. Zellweger et MuSICA ont donc glissé de 1.5 et
1.6 vers 1.6 et 1.7 ; l'ancien numéro est rappelé entre parenthèses. Chaque ligne nomme aussi
son fichier, qui lui ne bouge pas.

| Fig. | Source | Crédit en légende | Formalité |
|---|---|---|---|
| 1.1 | figure de communication ANR MaCCMic (`Fig_intro_maccmic.jpg`) | « Figure from the ANR MaCCMic project » | accord de Jérôme Ogée (coordinateur) ; demander aussi une version haute résolution, celle en anglais plafonne à 768 px |
| 1.2 | De Frenne et al. 2021, *Global Change Biology*, Fig. 4 (`Fig_intro_vertical_profiles.png`) | cité | Wiley, RightsLink, gratuit en thèse |
| 1.3 | Bouwen 2025, thèse, Fig. I.3 panneau A, d'après le Département de la santé des forêts (`Fig_intro_regeneration.png`) | cité | source primaire = DSF (service public) ; prévenir Klara par courtoisie |
| 1.6 (ex-1.5) | Zellweger et al. 2019, *Trends in Ecology & Evolution*, Fig. 1 (`Fig_intro_statistical_route.png`) | cité | Elsevier, RightsLink, gratuit en thèse |
| 1.7 (ex-1.6) | diagramme général de MuSICA (`Fig_intro_musica_scheme.jpg`) | « Figure by Jérôme Ogée, author of MuSICA » | accord de Jérôme Ogée |

RightsLink : depuis la page de l'article, « Get rights and content » → *reuse in a
thesis/dissertation*. La licence est délivrée immédiatement et sans frais.

**Écartée** : Bramer et al. 2018, *Advances in Ecological Research*, Fig. 1 (les facteurs
forçants et leur étendue verticale). Extraite et regardée : niveaux de gris,
photo-composite, style daté. Elle jurerait avec les onze figures dessinées, et la prose de
§1.2.1 tient sans elle. Le fichier reste dans `figures/_src/bramer-000.png` si tu veux la
récupérer.

## Resolution des figures empruntees (mesure du 11/09/2026)

Resolution effective a la taille imprimee (bloc de texte 17 cm, marges 2,4 cm) :

| Fig. | pixels | largeur imprimee | dpi effectif | verdict |
|---|---|---|---|---|
| 1.1 | 768 x 342 | 14,3 cm | **137** | sous le seuil de depot |
| 1.2 | 912 x 912 | 14,3 cm | **162** | sous le seuil de depot |
| 1.3 | 918 x 799 | 12,3 cm | **189** | sous le seuil de depot |
| 1.6 (ex-1.5, Zellweger) | 1634 x 2271 | 10,7 cm | 388 | acceptable |
| 1.7 (ex-1.6, MuSICA) | 2321 x 1615 | 15,2 cm | 387 | acceptable |

**La re-extraction est epuisee, ce n'est pas un defaut de decoupe.** Verifie avec
`pdfimages -list` sur les PDF sources :

- De Frenne et al. 2021 (Fig. 1.2) : le PDF editeur contient un raster **912 x 912 a
  145 ppi** en page 24, exactement le fichier deja utilise. Il n'existe rien de plus fin
  dans le PDF.
- Bouwen 2025 (Fig. 1.3) : la these contient un raster **1458 x 689 a 231 ppi** en page 17 ;
  le decoupage actuel vient d'un rendu de page a 300 dpi, donc deja au-dessus de la source.
- MaCCMic (Fig. 1.1) : version anglaise plafonnee a 768 px, pas de source vectorielle sous
  la main.

**Action, a joindre a la demande d'autorisation (le courrier est de toute facon a envoyer) :**
demander le fichier source (vectoriel ou raster haute resolution) aux trois interlocuteurs,
Jerome Ogee pour 1.1 et 1.7, Pieter De Frenne ou Wiley pour 1.2, Klara Bouwen pour 1.3.
A defaut, reduire la largeur imprimee : 1.2 a 70 % donne 204 dpi, 1.3 a 60 % en donne 239,
au prix de la lisibilite.
