---
title: "Chapitre 1 : passe CSI 3A"
subtitle: "manuscript_chap1_EN_native20_csi3pass_2026-09-02.md"
date: "2 septembre 2026"
geometry: margin=2.2cm
mainfont: DejaVu Serif
fontsize: 10pt
---

# Ce qui a été fait

- Fichier source : `manuscript_chap1_EN_native20_cites_2026-08-27.md`, celui qui a produit
  `manuscript_chap1_EN_native20_AoFS_2026-08-30.docx` (vérifié par quatre sondes de texte).
- Sortie : `manuscript_chap1_EN_native20_csi3pass_2026-09-02.md`. Le source n'est pas modifié.
- **61 corrections**, toutes au niveau de la phrase. Aucun paragraphe réécrit depuis zéro.
- Scripts dans `csi3pass/`, chaque correction étiquetée par le défaut qu'elle traite.

**Pourquoi pas le mode REWRITE du skill.** Le chapitre est déjà passé par `_voicepass_` et
`_rewrite_` en août, sa médiane est à 20 mots par phrase, il porte 76 clés de citation auditées et
521 valeurs numériques. Reconstruire les paragraphes depuis la science, c'est là qu'un chiffre ou
une référence se perd. La passe corrige la liste des défauts et laisse le reste.

# Résultat mesuré

| | avant | après |
|---|--:|--:|
| phrases de plus de 35 mots | 36 | **8** |
| phrase la plus longue | 84 mots | 59 mots (fragment de légende) |
| prose maniérée | 32 | **3** |
| pivots de négation | 4 | **0** |
| énumérations compressées | 1 | **0** |
| médiane mots par phrase | 20 | 19 |

Intégrité vérifiée mécaniquement, avant contre après : 521 valeurs décimales, 54 clés `[@...]`,
22 clés `@` en ligne, 45 renvois `Section N`, 39 renvois `Appendix X`, 48 renvois de figure et
68 occurrences de `°C`, tous inchangés. Le fichier compile.

# La décision de cadrage, à valider ou à refuser

Le chapitre construit son message autour du mot **`lever`**, et de deux mots qui l'accompagnent,
`exercise` et `potent`. Ces deux-là sont de la prose maniérée : un levier ne s'exerce pas, il a une
valeur. J'ai gardé `lever` comme nom technique et supprimé les deux autres.

| avant | après |
|---|---|
| The vertical profile is a mechanistic lever real canopies barely exercise. | The vertical profile has a mechanistic effect, and real canopies vary too little in shape to realize it. |
| The vertical profile is potent only when fully exercised | The vertical profile acts only at its full range |
| the design gradient exercises only a small part of either range | the design gradient spans only a small part of either range |
| The mechanism itself is potent | The mechanism itself is measurable. Isolated as a real-versus-uniform contrast, the profile moves ΔT~max~ by up to 0.14 °C |
| the profile proves potent in principle but barely exercised | the profile is large as a contrast and small over the range real stands span |

L'alternative serait de remplacer `lever` par `effect` ou `sensitivity` partout, ce qui est plus
strict mais fait environ 25 modifications dans un texte que vos encadrants ont lu trois fois avec
`lever` dedans, et fait perdre sa forme au Key message. Dites-moi si vous préférez cette option.

# Les corrections par défaut

## Prose maniérée, 17 corrections

`carries` (5), `potent` / `potency` (6), `exercise` / `exercised` (5), `buys`, `firms up`,
`room to move`, `fall as bins`, `credited with primacy`, `keeps mattering`, `does not survive`,
`margin to clear` (2), `sits` (3).

Exemples :

| avant | après |
|---|---|
| The coupling buys a more complete canopy-atmosphere representation | The coupling gives a more complete canopy-atmosphere representation |
| The profile contrast firms up in the intermediate archetypes | The profile contrast grows in the intermediate archetypes |
| Leaf quantity also has the most room to move within a structural type | Leaf quantity also varies most within a structural type |
| The archetypes fall as bins along it | The archetypes are bins along that function |
| An optical estimate must resolve leaf quantity well inside that margin | An optical estimate must resolve leaf quantity to better than 23% |
| A lever is read at the operating point where its archetype sits | A lever is read at the operating point of its archetype |

Trois occurrences de `carries` sont conservées : « the reanalysis carries a diurnally structured
warm bias » et deux tournures où `sit` décrit une position physique réelle du feuillage.

## Pivots de négation, 4 corrections

Deux rhétoriques, supprimés jusqu'à l'affirmative :

| avant | après |
|---|---|
| That ordering is built in rather than recovered (Section 2.3). | That ordering is built in (Section 2.3). |
| P1 is therefore not a young stand but a heterogeneous open class mixing… | P1 is therefore a heterogeneous open class, mixing… |

Deux défensifs, restructurés en parenthèse plutôt que supprimés, l'alternative écartée étant une
lecture qu'un relecteur ferait vraiment :

| avant | après |
|---|---|
| a lever near a bound is understated rather than inflated | a lever near a bound is understated (never inflated) |
| The sign flip corroborates the mechanism rather than signaling an artifact. | The sign flip corroborates the mechanism (it is not a solver artifact). |
| Saturation and truncation compress the dense end rather than inflate it. | Saturation and truncation compress the dense end (they do not inflate it). |

## L'énumération compressée de 84 mots, §3.2

C'était la phrase la plus longue du chapitre et la seule à cumuler deux-points, deux
points-virgules, trois parenthèses et quatre faits. Découpée en cinq phrases, sans tableau ni
liste, pour ne pas changer la structure du chapitre.

> **Avant.** The vertical arrangement of foliage is a potent lever, but only when fully exercised.
> Isolated as the real-versus-uniform contrast at fixed leaf area and height, it acts unevenly along
> the density gradient: it is small in the open P1 (0.01 °C); across the intermediate P2 and P3 it
> moves ΔT~max~ by up to 0.14 °C (P3, 0.14 °C; P2, 0.09 °C), exceeding the leaf-area step in 79% of
> P3 plots but only 58% of P2 plots; and in the dense P4 it is dwarfed by leaf quantity (0.08 versus
> 0.26 °C), its sign flipping from cooling to warming.

> **Après.** The vertical arrangement of foliage is a large lever at its full range. Isolated as the
> real-versus-uniform contrast at fixed leaf area and height, it acts unevenly along the density
> gradient. It is small in the open P1 (0.01 °C). Across the intermediate archetypes it moves
> ΔT~max~ by 0.09 °C in P2 and 0.14 °C in P3, exceeding the leaf-area step in 58% of P2 plots and
> 79% of P3 plots. In the dense P4 leaf quantity is three times larger (0.26 against 0.08 °C), and
> the profile contrast changes sign from cooling to warming.

## Phrases longues, 39 corrections

De 84, 79, 64, 59, 53, 51, 48, 47 mots et ainsi de suite jusqu'à 38. Chacune coupée en deux ou
trois, sans rien retirer. Le chapitre gagne 55 phrases et 15 mots au total : ce n'est pas une coupe,
c'est un redécoupage.

Les huit phrases qui restent au-dessus de 35 mots font 36 ou 37 mots et sont en Méthodes, où le
contenu paramétrique est irréductible. Les couper davantage donnerait une prose télégraphique.

# Ce qui n'a pas été touché

- Les cinq renvois `(Corroyez et al., in preparation)` vers le chapitre 3.
- Les légendes de figures et les tableaux. La règle du skill les plafonne à trois phrases et 60
  mots ; la légende de la Figure 4 en fait 190. C'est une passe à part, à décider.
- Les annexes A à H, soit environ 6 500 mots. La passe porte sur le corps, du Key message à la
  Conclusion.

# Suite

1. Lire les cinq passages du cadrage `lever` ci-dessus et dire si vous les gardez.
2. Reconstruire le docx AoFS depuis le nouveau `.md`, comme le 30 août.
3. Si vous voulez, je passe les légendes de figures, puis les annexes.
