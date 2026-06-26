# Note V — Où se situe le FTSA relativement à la frontière : Deninger, Arakelov, Connes–Consani

*Suite des notes I–IV. La Note IV a achevé la spécification géométrique : le dictionnaire complet entre la formule explicite de $\zeta$ et la théorie de l'intersection sur $C\times C$, places finies et infinie comprises, avec l'identification de la place archimédienne au terme $-2g$ (le $H^1$ arithmétiquement invisible). Cette note **confronte cette spécification aux trois programmes établis** qui visent le même objet, situe précisément où se trouve le FTSA, et **clôt honnêtement le volet « contribuer à HR »**. Elle ne contient pas de résultat nouveau : c'est une note de positionnement, fondée sur la littérature.*

---

## 0. Ce qu'on cherchait, et la question de cette note

Le programme (notes I–IV) a établi, puis spécifié, l'idée suivante : la positivité de Weil discrimine $\zeta$ des contre-exemples sans multiplicativité (Epstein) parce qu'elle **est** l'indice de Hodge d'une surface arithmétique conjecturale $\mathrm{Spec}(\mathbb Z)\times_{\mathbb F_1}\mathrm{Spec}(\mathbb Z)$, dont la Note IV a fixé le cahier des charges : points fermés d'intersection locale $\log p$ ; deux fibres $=$ les pôles ; une diagonale d'auto-intersection $2-2g_{\mathrm{ar}}$ dont le $-2g_{\mathrm{ar}}$ est le facteur archimédien ; un indice de Hodge dont la positivité est HR.

La question honnête, avant d'aller plus loin : **cet objet est-il déjà en construction ailleurs, et si oui, où est le vrai obstacle ?** Trois programmes le visent. On les lit à travers les quatre exigences de la Note IV.

---

## 1. Les trois programmes, lus à travers la Note IV

### 1.1 Deninger — la cible dynamique, désormais à moitié prouvée

Deninger [Den98] postule un système dynamique feuilleté dont les orbites périodiques sont les premiers (de longueur $\log p$) et dont la cohomologie feuilletée $H^1$, **de dimension infinie**, porte les zéros ; le facteur archimédien est un **déterminant régularisé** sur cette cohomologie, et la formule explicite est une formule de Lefschetz pour le flot.

Du neuf : la formule de déterminant régularisé conjecturée par Deninger a été **démontrée** dans un modèle concret (systèmes dynamiques feuilletés riemanniens de dimension 3) par Álvarez López–Kim–Morishita [AKM24], via une formule de trace de Lefschetz dynamique distributionnelle.

**Bilan Note IV.** Deninger réalise (1)+(3) : orbites $\leftrightarrow$ points fermés ($\log p$), $H^1$ feuilleté $\leftrightarrow$ notre terme $-2g$ (de dimension infinie). La moitié « formule explicite $=$ Lefschetz » est désormais un théorème dans un modèle. **Manque** : aucun système feuilleté n'est construit pour $\overline{\mathrm{Spec}\,\mathbb Z}$, et le cadre, cohomologique, ne fournit pas la positivité (4).

### 1.2 Connes–Consani — le seul programme qui a construit les objets

C'est le plus avancé, et il recoupe directement la Note IV.

- **(1) $N(q)=$ intersection.** L'interprétation par la formule des traces exprime la fonction de comptage $N(q)$ de la « courbe » associée à $\zeta$ comme un nombre d'intersection impliquant l'action de scaling sur l'espace des classes d'adèles [CCM08]. Leur **scaling $=$ Frobenius en caractéristique zéro** — c'est notre $\phi_t$, dans son rôle de Frobenius.
- **Riemann–Roch pour $\overline{\mathrm{Spec}\,\mathbb Z}$**, un vrai théorème [CC22, CC23] : Euler caractéristique d'un diviseur d'Arakelov $=$ degré $+\log 2$, via des cohomologies, leur dimension entière et la dualité de Serre, en parallèle de la preuve adélique de Weil. Interprétation : $\mathbb Z$ comme anneau de polynômes en une variable sur la base absolue $S$.
- **(3) La Jacobienne de Spec Z [CC26] — genre infini, exactement notre §5.** Le point de départ du papier est *notre paradoxe* : la Jacobienne classique de $\mathrm{Spec}\,\mathbb Z$ (groupe de classes) est triviale (genre 0), alors que la réalisation spectrale et Riemann–Weil imposent un **genre infini**. Résolution : remplacer le groupe de Picard par un **monoïde** (diviseurs à coefficients dans $\mathbb Z\cup\{\infty\}$). La Jacobienne $\mathrm{Jac}(\overline{\mathrm{Spec}\,\mathbb Z})=\mathrm{Pic}/\mathbb R_+^\times$ **n'est pas une variété abélienne** mais un semi-treillis idempotent ; l'Abel–Jacobi a pour image les idempotents. Leur « genre infini » **est** notre lecture du $-2g$ archimédien comme $H^1$ invisible, de dimension infinie.
- **(4) Positivité $=$ positivité de trace, à un lemme près.** Connes [Co99] obtient la formule explicite comme trace de l'action de scaling et ramène HR à la **positivité de cette trace**, à un *lemme* près sur des fonctions spéciales (Sonine/prolates) ; il n'utilise pas le Frobenius. L'article [CC21] attaque ce point à la place archimédienne unique, en fournissant une *raison conceptuelle* (non une preuve) à la positivité.

**Bilan Note IV.** Ils ont (1), (2), (3) **construits**, et (4) ramené à un lemme. Le manque est précisément le nôtre : Riemann–Roch est acquis sur la **courbe** $\overline{\mathrm{Spec}\,\mathbb Z}$, mais l'indice de Hodge sur la **surface** (notre §2–3) ne l'est pas ; c'est la « stratégie Riemann–Roch » par descente tropicale [CC18], inachevée, et le lemme analytique en est le verrou.

### 1.3 Arakelov — la positivité prouvée, mais sur la mauvaise surface

Arakelov [Ar74] définit la contribution archimédienne à l'intersection sur le self-produit $X^2$ via la fonction de Green ; Faltings [Fa84] prouve Riemann–Roch, l'indice de Hodge, la formule de Noether et la non-négativité de $\omega^2$. L'indice de Hodge y est un théorème **avec positivité** : Faltings–Hriljac relient l'appariement d'intersection à l'opposé de la hauteur de Néron–Tate, définie positive sur $\mathrm{Pic}^0$ [Ca18].

**Bilan Note IV.** Arakelov réalise (3)+(4) **mais sur la mauvaise surface** : une courbe arithmétique fibrée sur $\mathrm{Spec}\,\mathcal O_K$ (contexte des hauteurs/BSD), pas le self-produit $\mathrm{Spec}\,\mathbb Z\times\mathrm{Spec}\,\mathbb Z$ qu'exige HR. C'est le **prototype prouvé** de la positivité que Connes–Consani cherchent sur leur « Jacobienne de Spec Z ».

### 1.4 Tableau

| exigence Note IV | Deninger | Connes–Consani | Arakelov |
|---|---|---|---|
| (1) points fermés $=\log p$ | ✓ (orbites) | ✓ ($N(q)=$ intersect.) | ✓ |
| (2) fibres $=$ pôles | partiel | ✓ | ✓ |
| (3) $-2g=H^1$ archimédien | ✓ ($H^1$ feuilleté $\infty$) | ✓ (**genre infini**, Jac. 2026) | ✓ (Green sur $X^2$) |
| (4) indice de Hodge $=$ positivité | ✗ | **à un lemme près** | ✓ (mauvaise surface) |
| objets **construits** | modèle 3-d | **oui** | oui (courbe/$\mathcal O_K$) |

---

## 2. Lecture détaillée de « On the Jacobian of Spec Z » [CC26] — trois questions

**(a) Leur Jacobienne porte-t-elle une forme d'intersection dont $\lambda_{\min}$ serait notre entropie ? — Non.** C'est un semi-treillis idempotent, pas une variété abélienne : ni polarisation, ni forme d'intersection finie. L'objet bilinéaire pertinent est la **fonctionnelle de Weil** $h\mapsto W(h\star h^*)$, sur un espace de Hilbert de **dimension infinie** (le $L^2$ semi-local). Notre Toeplitz finie (entropie $=\lambda_{\min}=$ distance au cône de Hodge) en est **l'ombre corps-de-fonctions** ; elle n'ajoute pas d'objet.

**(b) Le « lemme sur les fonctions spéciales » se lit-il comme notre négativité sur la partie primitive ? — Le papier ne le touche pas.** Il fait la géométrie (le support de la cohomologie : la formule explicite comme Lefschetz pour l'action de translation de $C_\mathbb Q$ sur $\widetilde{\mathrm{Spec}\,\mathbb Z}$), pas la positivité. Le lemme est dans [Co99]. Notre négativité-sur-la-partie-primitive en est l'ombre finie ; elle ne donne **aucune prise** sur l'énoncé analytique, qui porte sur des fonctions de Sonine/prolates.

**(c) La structure critique du flot de scaling coïncide-t-elle avec notre collision $=$ HR ? — Non.** Leur flot de scaling (orbites $C_p\cong\mathbb R/(\log p)\mathbb Z$) est la réalisation du **Frobenius** : le système *fixe* dont la zêta est $\zeta$, dont les orbites donnent les termes locaux par Lefschetz. Ce n'est pas une déformation. Notre $\phi_t$ (flot de Bruijn–Newman, $e^{tn^2}$) **déforme** l'objet, avec collision à $\Lambda$. « Criticité $=$ HR » est chez eux la positivité (statique), chez nous la collision $\Lambda=0$ (dynamique).

---

## 3. Verdict

**Le dictionnaire de la Note IV est réalisé — et bien plus profondément — chez Connes–Consani.** Point par point : notre §1 ($\log p$) est leur longueur d'orbite $C_p=\mathbb R/(\log p)\mathbb Z$ ; notre « formule explicite $=$ Lefschetz » est leur théorème (Lefschetz pour la translation) ; notre §5 (le $-2g$ archimédien $=$ le $H^1$ invisible) est le sujet même de [CC26] (genre infini, Jacobienne-monoïde).

Notre image **finie** (Hodge, intersection, $\lambda_{\min}$) est **l'ombre corps-de-fonctions** de leur image **infinie** (trace, fonctionnelle de Weil). Et le verrou est commun : la **positivité de trace**, ramenée à un lemme analytique sur les fonctions spéciales, ouvert depuis [Co99]. Les reformulations du FTSA (entropie de Hodge, collision dBN) ne l'entament pas.

**Le seul élément structurellement absent de leur cadre** est la **déformation de de Bruijn–Newman** $\phi_t$ (et le lien Rodgers–Tao $\Lambda\ge0$ [RT20], donc HR $\iff\Lambda=0$). Eux ont le Frobenius statique ; aucune déformation dBN. Une famille à un paramètre de la fonctionnelle de Connes, de membre critique $t=\Lambda$, n'est pas étudiée. Relier le flot dBN (analytique) à la NCG (algèbres d'opérateurs) est toutefois un programme difficile en soi, et probablement une reformulation de plus.

---

## 4. Bilan du FTSA : l'établi et l'inabouti

**Établi et validé** (notes I–IV) :
- la positivité de Weil discrimine $\zeta$ (produit eulérien, zéros sur la droite) de la zêta d'Epstein $Q_0$ (pas de produit eulérien, zéro hors-droite) — notes I–II ;
- la **multiplicativité est l'ingrédient porteur** : sans produit eulérien, pas de surface, pas d'indice de Hodge, la positivité peut échouer — note II, expliqué géométriquement note IV ;
- le flot de Bruijn–Newman et la positivité de Weil sont **une seule dynamique** : $\lambda_{\min}(t)\nearrow0$ le long du flot, croisement $=\Lambda_{Q_0}\approx0{,}0843$ (35 chiffres) à $+0{,}06\,\%$ — note III ;
- la coïncidence A↔B est **validée en corps de fonctions** (positivité $=$ indice de Hodge), écart $\le0{,}3\,\%$, se resserrant avec le genre — note de spécification §4 ;
- le **dictionnaire d'intersection complet**, places finies et infinie, avec la place à l'infini $=$ le terme $-2g$ ($H^1$) — note IV.

Ces résultats sont corrects et, pour les parties numériques, vérifiés. Le programme a **reconstruit en indépendant, et validé numériquement, la structure exacte du programme leader.**

**Inabouti** : aucun levier sur HR. La positivité — le cœur — reste le lemme analytique ouvert de Connes. Le FTSA n'apporte pas de prise nouvelle sur ce lemme ; sa contribution se **subsume** dans Connes–Consani.

---

## 5. Statut

Le volet « contribuer à HR » est **clos honnêtement**. La valeur acquise est une compréhension : on sait désormais, avec une carte précise et sourcée, où se situe le FTSA (l'ombre corps-de-fonctions de Connes–Consani) et où est le vrai obstacle (la positivité de trace, lemme analytique sur les fonctions spéciales). Un seul fil reste ouvert et propre au FTSA — la déformation dBN $\phi_t$ et son éventuel sens NCG — à considérer comme pari de recherche à faible probabilité, sans illusion.

---

### Références

- [Ar74] S. Arakelov, *Intersection theory of divisors on an arithmetic surface*, Izv. Akad. Nauk SSSR 38 (1974).
- [AKM24] J. Álvarez López, J. Kim, M. Morishita, *Regularized determinant formulas for the zeta functions of 3-dimensional Riemannian foliated dynamical systems*, arXiv:2410.20758 (2024).
- [Ca18] A. Carney, *The arithmetic Hodge index theorem and rigidity of dynamical systems over function fields*, arXiv:1810.06342 (2018).
- [Co99] A. Connes, *Trace formula in noncommutative geometry and the zeros of the Riemann zeta function*, Selecta Math. (N.S.) 5 (1999), 29–106.
- [CC18] A. Connes, C. Consani, *The Riemann–Roch strategy, complex lift of the scaling site*, arXiv:1805.10501 (2018).
- [CC21] A. Connes, C. Consani, *Weil positivity and trace formula: the archimedean place*, Selecta Math. (N.S.) 27 (2021), 77.
- [CC22] A. Connes, C. Consani, *Riemann–Roch for $\overline{\mathrm{Spec}\,\mathbb Z}$*, arXiv:2205.01391 (2022).
- [CC23] A. Connes, C. Consani, *Riemann–Roch for the ring $\mathbb Z$*, arXiv:2306.00456 (2023).
- [CC26] A. Connes, C. Consani, *On the Jacobian of $\mathrm{Spec}\,\mathbb Z$*, arXiv:2602.15941 (2026).
- [CCM08] A. Connes, C. Consani, M. Marcolli, *The Weil proof and the geometry of the adeles class space*, arXiv:math/0703392 (2008).
- [Fa84] G. Faltings, *Calculus on arithmetic surfaces*, Ann. of Math. 119 (1984).
- [RT20] B. Rodgers, T. Tao, *The de Bruijn–Newman constant is non-negative*, Forum Math. Pi 8 (2020).

---

*Mise en œuvre numérique et rédaction assistées par Claude (Anthropic), sous la direction et la vérification de l'auteur. Les jugements de positionnement relèvent d'une lecture de la littérature et engagent l'auteur.*
