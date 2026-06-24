# Note IV — Le dictionnaire d'intersection : de la formule explicite à $C\times C$, et ce que doit être la place à l'infini

*Suite des notes I–III et de la note de spécification géométrique. La phase de validation-par-l'ombre étant close (la coïncidence A↔B est établie côté nombre et côté corps de fonctions), cette note **assemble le dictionnaire complet** entre la formule explicite de $\zeta$ et la théorie de l'intersection sur $C\times C$, places finies **et** infinie comprises. Elle précise, et **spécifie entièrement**, le contenu de la contrainte D3.*

*Le contenu mathématique de cette note est, en corps de fonctions, un **théorème** (Lefschetz, Künneth, indice de Hodge de Weil 1948). Notre apport est de le rendre concret sur de vraies courbes et de l'assembler en un cahier des charges pour le cas arithmétique. Le **transport vers $\mathbb Z$** reste conjectural — c'est le programme $\mathbb F_1$.*

---

## 0. Position du problème

La thèse du programme (notes I–II) est que la **positivité de Weil** discrimine $\zeta$ (produit eulérien, zéros sur la droite) des contre-exemples sans multiplicativité (Epstein, Davenport–Heilbronn). La note de spécification a localisé l'origine de cette positivité : en corps de fonctions, elle **est** l'indice de Hodge sur la surface $C\times C$, et la multiplicativité (le produit eulérien sur les points fermés) **est** ce qui fait exister la surface et sa théorie d'intersection.

Reste à savoir *exactement* quels objets de la formule explicite correspondent à quels objets géométriques. C'est l'objet de cette note. On procède en deux temps : le côté **premier** (les facteurs locaux), puis le côté **non-premier** (pôles et facteur archimédien), qui est le plus délicat car c'est là que vit la place à l'infini.

On fixe une courbe $C$ lisse projective de genre $g$ sur $\mathbb F_q$, de fonction zêta
$$
Z_C(T)=\frac{P(T)}{(1-T)(1-qT)},\qquad P(T)=\prod_{i=1}^{2g}(1-\alpha_iT),\quad |\alpha_i|=\sqrt q,
$$
et la surface $S=C\times C$, avec le diagonal $\Delta$, les graphes de Frobenius $\Gamma_{F^n}=\{(x,F^n x)\}$, et les deux fibres $f_1=[C\times\mathrm{pt}]$, $f_2=[\mathrm{pt}\times C]$.

---

## 1. Côté premier : facteur local $=$ intersection locale $=\deg(x)$

Par la formule des traces de Lefschetz,
$$
(\Gamma_{F^n}\cdot\Delta)=\#\mathrm{Fix}(F^n)=\#C(\mathbb F_{q^n})
=\mathrm{tr}(F^n|H^0)-\mathrm{tr}(F^n|H^1)+\mathrm{tr}(F^n|H^2)
=q^n+1-\sum_i\alpha_i^n .
$$
Les points fixes de $F^n$ se regroupent en **points fermés** (orbites de Galois) : un point fermé $x$ de degré $d\mid n$ apporte ses $d$ points géométriques. D'où la **décomposition locale** du nombre d'intersection,
$$
(\Gamma_{F^n}\cdot\Delta)=\#C(\mathbb F_{q^n})=\sum_{d\mid n}d\,a_d
=\sum_{x:\ \deg x\mid n}\deg(x),\qquad a_d=\#\{\text{points fermés de degré }d\},
$$
chaque point fermé contribuant **$\deg(x)$**. Or $\deg(x)$ est précisément le poids de la place $x$ dans la formule explicite (le von Mangoldt $\Lambda(x^m)=\deg x$). On a donc l'identité
$$
\boxed{\ \text{facteur de Weil local en }x\ =\ (\Gamma_{F^n}\cdot\Delta)_x\ =\ \deg(x)\ }
$$
— la formule explicite **est** une identité de nombres d'intersection, terme local par terme local.

**Vérification numérique** (`fork2_local_intersection.py`, courbes elliptiques réelles $q=5,7,11$, comptage du corps de base seul) : pour tout $n$, $\ \#C(\mathbb F_{q^n})=\sum_{d\mid n}d\,a_d\ $ exactement, avec $a_d\ge0$ entiers (points fermés authentiques) et $|\alpha_i|=\sqrt q$ (HR de Weil pour la courbe).

**Dictionnaire (côté premier).** Point fermé $x\leftrightarrow$ premier $p$ ; $\deg(x)\leftrightarrow\log p$ ; $(\Gamma_{F^n}\cdot\Delta)_x\leftrightarrow$ facteur local de la formule explicite en $p$.

---

## 2. Côté non-premier : pôles $\leftrightarrow$ fibres, facteur $\Gamma$ $\leftrightarrow$ genre

La formule explicite a, outre les termes premiers, des **termes de pôles** (les pôles de $\zeta$ en $s=0,1$) et un **terme archimédien** (le facteur $\Gamma$). On les lit sur l'auto-intersection de la diagonale.

### 2.1 Nombres d'intersection sur $C\times C$ (via Künneth)

Décomposition de Künneth de $H^2(C\times C)=H^0\!\otimes H^2\ \oplus\ H^1\!\otimes H^1\ \oplus\ H^2\!\otimes H^0$ :
$$
[\Gamma_{F^n}]=q^n\,f_1+f_2+\gamma_n,\qquad \gamma_n\in H^1\!\otimes H^1,
$$
où $\gamma_n$ est le graphe de $F^n$ vu comme endomorphisme de $H^1$. Les appariements :
$$
f_1^2=f_2^2=0,\quad f_1\cdot f_2=1,\qquad
\gamma_\phi\cdot\gamma_\psi=-\,\mathrm{Tr}\!\big(\phi\,\psi'\big),
$$
où $\psi'$ est l'**adjoint symplectique** (involution de Rosati) pour la forme d'intersection (alternée) sur $H^1$, vérifiant pour Frobenius $\ F'=qF^{-1}\ $ (car $F$ multiplie la forme symplectique par $q$). On en déduit, sans ambiguïté de signe (calibré sur la diagonale) :
$$
\Gamma_{F^n}\cdot f_1=q^n,\qquad \Gamma_{F^n}\cdot f_2=1,
$$
$$
\Gamma_{F^n}\cdot\Gamma_{F^m}=q^{\min(n,m)}\,N_{|n-m|}\ \ (n\ne m),\qquad
\Gamma_{F^n}^2=q^n(2-2g).
$$
En particulier, pour $n=0$ ($\Gamma_{F^0}=\Delta$) :
$$
\boxed{\ \Delta^2=2-2g\ }=\underbrace{\;2\;}_{2\,(f_1\cdot f_2)}+\underbrace{(-2g)}_{\gamma_0^2=-\mathrm{Tr}(\mathrm{id}\,|H^1)} .
$$

### 2.2 Lecture

La décomposition sépare exactement les deux blocs « triviaux » de la formule des traces :

- **$+2$ (les fibres)** $=2\,(f_1\cdot f_2)$ provient de $H^0$ et $H^2$, c'est-à-dire des **deux pôles** de $Z_C$ en $T=1$ et $T=1/q$ — les pôles de $\zeta$ en $s=0$ et $s=1$.
- **$-2g$ (le genre)** $=\gamma_0^2=-\mathrm{Tr}(\mathrm{id}|H^1)$ provient de la partie $H^1\!\otimes H^1$ — c'est le **terme archimédien** (le facteur $\Gamma$).

**Vérification numérique** (`fork2_archimedean_dictionary.py`) :

| courbe | $\Delta^2=2-2g$ | décomposition | signature de la forme d'intersection |
|---|---|---|---|
| $g=1$ (elliptique $q=5$) | $0$ | $+2-2$ | $(1,\,3^-,\,2^0)$ |
| $g=2$ (Weil $q=5$, $a_1=-4,a_2=10$) | $-2$ | $+2-4$ | $(1,\,5^-,\,2^0)$ |

---

## 3. L'indice de Hodge $=$ la positivité de Weil

La forme d'intersection sur $\mathrm{NS}(C\times C)$ a, par l'**indice de Hodge** (Hodge–Riemann pour les surfaces / Castelnuovo), une signature $(1,\rho-1)$ : **une seule** direction positive (la classe ample), toutes les autres négatives. La vérification ci-dessus le confirme : sur l'espace engendré par $\{f_1,f_2,\Gamma_{F^0},\dots,\Gamma_{F^k}\}$, on trouve **exactement une** valeur propre positive (les valeurs propres nulles ne sont que la redondance des puissances de Frobenius au-delà du degré $2g$).

La direction positive est portée par le plan hyperbolique des fibres (les pôles) ; la **partie primitive** (orthogonale aux fibres) est **définie négative**, et cette négativité **est** la positivité de Weil. Son cœur est la forme de Weil sur $H^1$ — la matrice de Toeplitz des moments normalisés $s_k=\sum_i(\alpha_i/\sqrt q)^k$, dont $\lambda_{\min}\approx0$ (le bord critique) pour une vraie courbe. C'est l'objet déjà rencontré (note de spécification, §4) : entropie $=\lambda_{\min}=$ distance au cône de Hodge.

---

## 4. Le dictionnaire complet

| objet analytique ($\zeta$) | corps de fonctions ($C/\mathbb F_q$) | géométrie sur $C\times C$ |
|---|---|---|
| zéros non triviaux $\rho$ | valeurs propres de $F$ sur $H^1$ | partie $H^1\!\otimes H^1$ de $\Delta$ |
| facteur local $\log p$ | point fermé $x$, $\deg(x)$ | $(\Gamma_{F^n}\cdot\Delta)_x=\deg x$ |
| pôles en $s=0,1$ | pôles de $Z_C$ en $T=1,1/q$ ($H^0,H^2$) | fibres $f_1,f_2$ ; $+2$ dans $\Delta^2$ |
| facteur $\Gamma$ (place $\infty$) | le genre (pas de place archimédienne) | $\gamma_0^2$ ($H^1$) ; $-2g$ dans $\Delta^2$ |
| équation fonctionnelle $s\leftrightarrow1-s$ | $\alpha_i\leftrightarrow q/\alpha_i$ | involution de Rosati $F'=qF^{-1}$ |
| positivité de Weil | indice de Hodge | signature $(1,\text{reste})$ sur $\mathrm{NS}(C\times C)$ |
| $\;$ | $\Delta^2=2-2g$ | $+2$ (pôles) $-2g$ (archim.) |

---

## 5. Ce qu'est la place à l'infini

Le côté non-premier dit quelque chose de net sur la question la plus difficile du programme $\mathbb F_1$. Dans le miroir corps de fonctions, **la place archimédienne de $\mathbb Q$ correspond au terme $-2g$ de l'auto-intersection de la diagonale**, c'est-à-dire au $H^1$ lui-même — aux **zéros**.

C'est cohérent avec la structure : en corps de fonctions, $H^1$ est de la cohomologie étale honnête (les $2g$ zéros sont *là*, géométriquement) et il n'y a *pas* de facteur $\Gamma$ séparé. En arithmétique, ce $H^1$ est **invisible** (on n'a pas la cohomologie), et le facteur $\Gamma$ archimédien en est exactement le **substitut** : il porte le $-2g$, le « genre arithmétique ». Autrement dit,

> la place à l'infini de $\mathbb Q$ est le lieu qui doit fournir le $H^1$ manquant.

C'est là que le programme $\mathbb F_1$ bute : construire $\mathrm{Spec}(\mathbb Z)\times_{\mathbb F_1}\mathrm{Spec}(\mathbb Z)$, c'est essentiellement **fabriquer ce $H^1$ archimédien**.

---

## 6. Ce que ça fixe pour D3

La contrainte D3 (« pourquoi la multiplicativité force la positivité ») est désormais **entièrement spécifiée**. La surface arithmétique $\mathrm{Spec}(\mathbb Z)\times_{\mathbb F_1}\mathrm{Spec}(\mathbb Z)$ doit posséder :

1. des **points fermés** $p$ d'intersection locale $\log p$ (côté premier, §1) ;
2. deux classes de **fibres** $=$ les pôles de $\zeta$ en $s=0,1$ (contribution $+2$, §2) ;
3. une **diagonale** d'auto-intersection $2-2g_{\mathrm{ar}}$, dont le $-2g_{\mathrm{ar}}$ **est** le facteur $\Gamma$ archimédien (le $H^1$ arithmétique, §5) ;
4. un **indice de Hodge** dont la signature $(1,\text{reste})$ est la positivité de Weil $=$ HR (§3).

On sait maintenant *exactement* à quoi doit ressembler l'objet ; il reste à le construire — et c'est là que le numérique s'arrête.

---

## 7. Établi / conjectural

| énoncé | statut |
|---|---|
| facteur local $=(\Gamma_{F^n}\cdot\Delta)_x=\deg x$ | **théorème** (Lefschetz + orbites de Galois) ; **vérifié** §1 |
| $\Delta^2=2-2g=+2-2g$ (pôles / genre) ; nombres d'intersection | **théorème** (Künneth, Rosati) ; **vérifié** §2 |
| positivité de Weil $=$ indice de Hodge, signature $(1,\text{reste})$ | **théorème** (Weil 1948) ; **vérifié** §3 |
| dictionnaire complet, place à l'infini $=$ terme $-2g$ ($H^1$) | **assemblage** (correct, non formalisé pour $\mathbb Z$) |
| transport vers $\mathbb Z$ : existence de la surface arithmétique | **conjectural** (programme $\mathbb F_1$), maintenant **entièrement spécifié** |
| le tout $\Rightarrow$ HR | **non** ; programme, pas preuve |

---

### Scripts

- `fork2_local_intersection.py` — vérifie §1 sur de vraies courbes (elliptiques, comptage du corps de base) : $|\alpha_i|=\sqrt q$, $a_d\ge0$ entiers, $\#C(\mathbb F_{q^n})=\sum_{d\mid n}d\,a_d$.
- `fork2_archimedean_dictionary.py` — vérifie §2–§3 (g=1 elliptique réel ; g=2 polynôme de Weil valide) : $\Delta^2=2-2g$ et sa décomposition, signature $(1,\text{reste})$ de la forme d'intersection.

---

*Mise en œuvre numérique et rédaction assistées par Claude (Anthropic), sous la direction et la vérification de l'auteur.*
