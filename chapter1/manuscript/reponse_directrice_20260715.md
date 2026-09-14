**Objet : Re: Chapitre 1 — profil vertical, rescaling LAI, capteurs HOBO**

Bonjour [Prénom],

Merci beaucoup pour ces retours, ils tombent juste et m'aident à resserrer l'argumentaire. Réponses point par point.

---

**1. Rescaling du LAI (k, θ, angles de tir)**

Je te suis : je fais le **calcul exact des angles de tir** (la formule de l'article S2/LiDAR) et je **contrôle contre les profils `lidR`** pour vérifier la cohérence (baisse légère attendue le long des profils, ou égalité si les tirs sont quasi verticaux). Je pars sur le **k optimisé**, puisque les valeurs de l'étude S2/LiDAR correspondent à cette résolution et restent dans la gamme littérature.

Un point rassurant pour le microclimat : un rescaling k/clumping est **quasi-uniforme spatialement**, donc il déplace surtout la *baseline* du LAI, pas la *hiérarchie* d'attribution (qui repose sur des perturbations ±1 SD locales, autour de chaque canopée). Je l'ai **vérifié** par un test de sensibilité (rescale uniforme k 0.5→0.65, soit ×0.769, sur les 4 archétypes) : sur ΔTmax, le levier LAI rétrécit proportionnellement mais **reste dominant et croît toujours P1→P4** en canopée fermée (P3, P4), le profil restant ≤0,03 °C. Seule exception, le type mi-ouvert **P2**, où le levier LAI est déjà quasi-nul (parité avec le profil, les deux ≤0,04 °C) — mais c'est la zone de parité que je documente déjà, pas un renversement. Bref : **le classement tient là où l'histoire se lit.** Je ferai le calcul d'angle de tir exact (Éq. 2 de notre papier S2/LiDAR, `l = dz·⟨1/cos θᵢ⟩`) avec contrôle `lidR`, et je peux en faire une comparaison des deux k dans le papier comme tu le suggères.

**2. « Seul le LAI compte, le profil ne sert à rien / nuit »**

Tu as raison de tiquer, et en creusant je te rejoins sur le fond : **le profil vertical a bien un mécanisme réel.** Mon test contrôlé H2 (Annexe C) montre qu'à LAI, hauteur et couverture **fixés**, concentrer le feuillage en haut refroidit le sous-bois **jusqu'à ~0,58 °C** dans le P4 dense. Donc le levier n'est pas nul.

Ce que je montre, c'est que ce **potentiel n'est pas réalisé** : les vrais peuplements se concentrent près de la mi-hauteur et varient trop peu en équilibre vertical pour l'exprimer. Du coup la variation ±SD *réalisée* ne bouge ΔTmax que de ≤0,06 °C.

Sur la Fig. 3 et la C2 : les effets marqués que tu vois **entre types** sont portés par la **densité**, pas par la forme du profil. La C2 traçait les **4 archétypes réels** (tous traits à leur valeur de centroïde), donc elle mélangeait densité + hauteur + couverture + forme — d'où l'ambiguïté. **Je l'ai corrigée exactement comme tu le proposes** : pour chaque type je superpose maintenant le profil réel (trait plein) et un **profil homogène de même LAI/hauteur/couverture** (pointillé). Résultat : les deux se superposent presque (**0,01–0,12 °C d'écart au niveau sous-canopée**), face à ~2,8 °C d'écart entre types. C'est H1 rendu visible : la séparation, c'est la densité, pas la structure verticale.

Enfin, je **corrige le registre du titre et du texte** : c'était trop tranchant (« seul le LAI », « inutile de mesurer en 3D ») par rapport à ma propre Annexe C. Je passe sur « la quantité de feuilles **davantage que** le profil vertical gouverne… » et je précise partout que le mécanisme vertical est **réel mais non réalisé** dans les canopées observées (small-but-non-null).

**3. HOBO sur troncs, face nord, à l'ombre**

Point important, je le prends pleinement. Les HOBO échantillonnent un microhabitat **près du tronc, à l'ombre**, qui n'est pas la même grandeur que la température d'air de canopée simulée par MuSICA (1D, à une hauteur). Une **partie du biais chaud (+0,65 °C) relève probablement de ce décalage de mesure**, pas seulement d'une ET/humidité du sol manquante — et contrairement au modèle statistique d'Eva (calé directement sur les HOBO), MuSICA ne peut pas l'absorber. Je l'ajoute donc explicitement dans la Discussion **comme co-cause du biais**, à côté de l'ET. Note : ce décalage étant surtout un *offset fixe*, il pèse sur le biais absolu mais pas sur le classement (qui repose sur la différence haut-bas de chaque sensibilité) — ce qui reste cohérent avec le reste.

---

Bilan : figure C2 refaite, titre/texte adoucis, caveat HOBO intégré, et je lance le run de sensibilité LAI rescalé pour la robustesse. Je te tiens au courant.

Pour se voir : je note que tu es absente **mardi 21** — on cale un créneau **sur le reste de la semaine prochaine** ? Je serai à la MTD cette semaine et la suivante avant de repasser en télétravail pour l'été.

Bonne semaine, et bon repos si tu coupes un peu,

Nathan
