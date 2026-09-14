envoyer code, réu, histoire, slides

Grosses diff S2/LiDAR sur distribution verticale uniforme ? capteur =/= niveau information
Tester MuSICA LAI LiDAR, Hauteur LiDAR, mais distribution uniforme (c hanger 1 variable par 1 variable, éventuellement voir si info LAI/Hauteur isolable), fCover 1 ou varie
Q1 S2 trop orientée S2 marche pas LiDAR marché
FORMS SPOT avec LiDAR HD ?

H1 S2 c'est trop un résultat -> mettre de côté ?
manque plan exp associés aux H
S2 mentionné que à la fin, focus sur la distribution verticale
redire Hyp que LiDAR valable tout l'été

H2 identification imp relative LAI Hmax fCover LAD (on a pas une idée si claire de la hiérarchie)
LAI fCover suffisant ? ou Hmax LAD nécessaire -> ça amène S2, comment on fait pour intégrer S2 sans les infos de hauteurs ? -> fusion, complémentarité ALS S2
ou alors ? LAI fCover ne suffisent pas, distrib joue, sur impact air...

hypothèses doivent être testables =/= interprétations de résultats
voir les outputs possibles de MuSICA en plus de Tair

H3, si je peux pas le montrer -> discussion, références
dire que y a un buffering + important, sans en expliquer les causes ?
végétation basse fait un isolant + épais ? circulation de l'air varie ?

effet jour dans le GAMM c'est quoi exactement ?
attention LAI 10m à 3 bandes vs 20m à 8 bandes

échantillonnage sur fCover Hauteur LAI et après seulement FPCA ?
légitime de faire FPCA avant ? représentativité biaisée

superposer time-series MuSICA S2, LiDAR, LiDAR uniforme, HOBOs
temp fortes disparaissent
peut-on regarder toutes les time-series ensemble ?

regarder par HOBO gamme température (sachant qu'à priori on a 45 buf 8 amp), eg 1 buf 1 amp

etre au point sur FPC, tracer fonctions (contribution des couches LAD aux FPC, un peu à la manière d'une vraie ACP où on dit métrique1 contribue pour X% à PC1, etc)
FPC1 vs FPC2, FPC1 vs FPC3, FPC2 vs FPC3
plotter dans espace composantes, pour chaque profil LAD on doit pouvoir dire que c'est X% FPC1, Y% FPC2, Z% FPC3

interprétation des FPCs est-ce que ça a vraiment un sens vis-à-vis de la végétation, et donc du DeltaTmax
se focaliser sur FPC1 ??? mais peu de varex
explorer les fonctions de plot du package fda
ou trouver un autre moyen que les FPCs pour expliquer profils

signification des composantes sur la végétation réelle compliqué à faire... est-ce grave ?

problème norm sur profils LAD ABCDEF

intercomparaison MuSICA LiDAR vrais profils vs MuSICA LiDAR profils uniforme

je me sur-complexifie la tâche -> il faut simplifier !!!!! pourquoi on regarde les profils ? typologie profils ? regarder MuSICA par type de profils/cluster

FPCA sur profils avec les hauteurs

points en suspens: avoir la forme du profil vertical, est-ce important ou pas, faut-il faire FPCA avant ou après échantillonnage ? profils type dans MuSICA et comparer

4 variables dérivées du LiDAR pour MuSICA: LAI, profil LAD, Hmax, fCover -> analyse de sensibilité / importance, analyse de manière + systématique, ce qui pourrait justifier les différences de résultat MuSICA-S2 / MuSICA-LiDAR, même sans les HOBOs, pour se focaliser sur MuSICA
varier l'information du LiDAR: introduire telle var, puis telle var...

data variable: jour, heure ? c'est super corrélé
points au hasard dans la série temporelle ?