# Spécification géométrique : choix de l'espace et stratégie du flot (note de synthèse)

*Note de travail (conceptuelle) — programme FTSA, sous-problème 1.*
*Suite des notes I–III. Cette note ne contient aucun résultat numérique nouveau : elle consigne des **décisions de programme** et la stratégie de test. La séparation établi / spéculatif des notes précédentes reste en vigueur.*

---

## 0. Objet

On cherche un triplet $(\mathcal C,\ \phi_t,\ \pi)$ — espace de configurations, flot, projection — réalisant géométriquement la dynamique dont on a établi les deux ombres :

| | desideratum | statut |
|---|---|---|
| **D1** | $\pi\circ\phi_t$ = flot de de Bruijn–Newman sur les zéros | ombre calculée (Note III) |
| **D2** | forme de Weil = fonctionnelle monotone le long de $\phi_t$ | ombre calculée (Note III) |
| **D1↔D2** | (A) flot dBN et (B) descente de Weil = *même* dynamique | **établi numériquement** (Note III) |
| **D3** | produit eulérien *structurel* dans $\mathcal C$ | contrainte dure (Note II) |
| **D4** | points critiques de $\phi_t$ = configs Weil-positives = HR | objectif |

La Note III a transformé D1↔D2 d'hypothèse en fait : tout candidat doit produire **un seul** flot se projetant à la fois sur dBN et sur la descente de $\lambda_{\min}$, avec le même temps critique $\Lambda_{Q_0}$.

---

## 1. Fork 1 — l'espace $\mathcal C$ : tranché

Critères décisifs : **K1** multiplicativité native (D3), **K2** non-circularité (ne pas présupposer HR). Critères requis : **K3** zéros $=\pi(\mathcal C)$, **K4** Weil = fonctionnelle naturelle, **K5** flot $\phi_t$ naturel, **K6** réduction traçable au cas corps de fonctions.

* **Triplet spectral synthétique** $(\mathcal A,\mathcal H,D)$ — l'idée Hilbert–Pólya. Échoue à **K1** (il faut *fabriquer* le produit eulérien dans le spectre — or la Note II dit que c'est lui le moteur) et à **K2** (on construit $D$ pour avoir les zéros voulus : circulaire). Bon pour K5 seulement.
* **Espace des classes d'adèles** $\mathbb A_{\mathbb Q}/\mathbb Q^\times$ et site des fréquences (Connes 1999 ; Connes–Consani). Les zéros y sont le **spectre d'absorption** de l'action des classes d'idèles, la formule explicite de Weil y est une **formule des traces**, et HR $\iff$ **positivité** de cette trace. **K1 ✓✓** (les adèles se factorisent en $\prod_p\mathbb Q_p\times\mathbb R$ — le produit eulérien est natif, et les termes de places sont exactement notre décomposition place-par-place de la Note II), **K2 ✓✓**, **K4 ✓✓** (Weil = la trace), **K3 ✓**, **K6 ✓**. Seul **K5 ✗** : pas de flot dans le tableau standard.

**Décision.** $\boxed{\mathcal C=\mathbb A_{\mathbb Q}/\mathbb Q^\times}$ (espace des classes d'adèles / site des fréquences), **dont on dérive** le triplet spectral — et non l'inverse. Le triplet spectral n'est pas un concurrent mais un produit dérivé de la construction adélique. Tout le gap est concentré en **K5 : la construction de $\phi_t$**.

---

## 2. Fork 2 — le flot et son entropie : deux faces d'un même objet

$\phi_t$ (le flot) et $\mathcal S$ (la fonctionnelle de Lyapunov / entropie à la Perelman, monotone le long de $\phi_t$, stationnaire ssi HR) ne sont pas deux options exclusives : on a besoin des deux.

**Piste $\phi$ — un piège de signe à éviter.** Le réflexe « $s$ dual au scaling, dBN = chaleur, donc $\phi_t=e^{\tau\partial_u^2}$ (laplacien de scaling) » est **faux**. dBN agit côté noyau par *multiplication* $\Phi(u)\mapsto e^{tu^2}\Phi(u)$ — une amplification **anti-diffusive** des grandes longueurs $u$ (les $\log p^k$). Le flot de la chaleur de scaling est au contraire une *convolution* gaussienne (lissage), qui amortit le côté spectral en $e^{-\tau r^2}$ : l'opération **inverse**. C'est précisément cette anti-diffusion qui rend $\Lambda\ge 0$ non trivial. Le bon $\phi_t$ n'est donc pas le laplacien de scaling naïf ; sa réalisation géométrique sur $\mathbb A_{\mathbb Q}/\mathbb Q^\times$ est une **repondération gaussienne** des termes de la trace, dont l'opérateur générateur canonique reste à identifier.

**Piste entropie — on a déjà une longueur d'avance.** L'entropie candidate est la **fonctionnelle de Weil elle-même** : $\lambda_{\min}(t)$, dont la Note III a montré la montée monotone vers $0$. Sous-problème 2 et D2 sont le même objet. Ce qui reste ouvert est la version *propre* : un **scalaire unique** à formule de monotonie de type Bochner $\dfrac{d\mathcal S}{dt}=\int|\nabla(\cdot)|^2\ge 0$, plutôt qu'une valeur propre. Ancrages connus : le **critère de Li** ($\lambda_n\ge0\ \forall n\iff$ HR), candidat naturel à une entropie indexée ; et, côté corps de fonctions, la positivité de Weil *est* une inégalité d'indice (Hodge sur $C\times C$).

---

## 3. Stratégie de test : descendre dans l'analogue corps de fonctions

Le risque d'auto-illusion est ici maximal (construire un flot qui « marche » parce qu'on y a glissé la réponse). La parade : **le cas des corps de fonctions**, où Weil a *prouvé* HR par une positivité (indice de Hodge / Castelnuovo–Severi sur $C\times C$). C'est le seul terrain où l'on peut **tester** la construction au lieu de la décréter, car :

* $\mathcal C \leadsto C\times C$ ;
* positivité de Weil $\leadsto$ inégalité de Hodge (théorème) ;
* les « zéros » $\leadsto$ les $2g$ valeurs propres de Frobenius $\alpha_i$, HR $\iff |\alpha_i|=\sqrt q$ (angles réels) ;
* tout est **fini-dimensionnel et exact** (polynômes, petites matrices) — pas de précision à gérer.

La question test : *quel est l'analogue du flot $\phi_t$, et la forme d'intersection est-elle l'entropie monotone ?* Si la coïncidence A↔B établie sur l'Epstein a un analogue géométrique, il doit y être **visible et vérifiable**.

**Plan de l'expérience corps de fonctions** (miroir exact de l'étape 2) :
1. **Ancrage (non-circularité).** Config tous-angles-réels (HR vrai) $\Rightarrow$ vérifier que la forme de Weil reconstruite est $\succeq 0$ — réduction à la positivité prouvée de Weil.
2. **Test.** Config avec un **quadruplet hors-cercle** $\pm\beta\pm i\delta$ (analogue corps-de-fonctions de Davenport–Heilbronn / de l'Epstein) $\Rightarrow$ $\lambda_{\min}<0$ ; puis flot dBN-analogue $\Rightarrow$ suivre $\delta(t)\to 0$ (collision $\Lambda_{\mathrm{ff}}$) **et** $\lambda_{\min}(t)\nearrow 0$, et vérifier qu'ils coïncident.
3. **Identification de l'entropie.** Montrer que la forme de Weil y est l'avatar de la forme d'intersection (positivité de Bochner / Hodge), donc que l'« entropie » a un sens géométrique — ce qui manquait dans le cas Epstein.

---

## 4. Établi vs conjectural (état à l'entrée de la phase de test)

| brique | statut |
|---|---|
| $\mathcal C=\mathbb A_{\mathbb Q}/\mathbb Q^\times$ ; zéros = spectre d'absorption ; Weil = trace | **établi** (Connes 1999) |
| multiplicativité structurelle dans $\mathcal C$ (D3) | **établi** (adélique) + appuyé (Note II) |
| $\lambda_{\min}(t)\nearrow0$ le long de dBN (D1↔D2) | **établi numériquement** (Note III) |
| $\phi_t$ canonique sur $\mathcal C$ se projetant sur dBN | **ouvert** ; candidat naïf (laplacien scaling) **écarté** |
| entropie scalaire à formule de Bochner | **ouvert** ; candidate = fonctionnelle de Weil |
| réduction corps de fonctions (Hodge sur $C\times C$) | **théorème** (Weil) — pierre de touche du test |
| le tout $\Rightarrow$ HR | **non** ; programme, pas preuve |

---

*Rédaction assistée par Claude (Anthropic), sous la direction et la vérification de l'auteur.*
