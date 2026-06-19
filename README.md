# Positivité de Weil comme discriminant — expériences numériques sur ζ et les zêta d'Epstein

Code Python reproductible et deux notes de recherche explorant numériquement la **positivité de Weil** comme critère séparant la fonction ζ de Riemann (à produit eulérien) des fonctions zêta d'**Epstein** de formes quadratiques binaires (sans produit eulérien). Le cas étudié est le formulaire principal $Q_0=m^2+mn+6n^2$, de discriminant $-23$ et de nombre de classes $3$.

---

## ⚠️ Ce que ce dépôt est — et n'est pas

**Ce n'est PAS une preuve de l'hypothèse de Riemann**, ni une voie revendiquée vers une preuve.

Les faits mathématiques sous-jacents sont **classiques** :
- les fonctions zêta d'Epstein de nombre de classes $>1$ ont des zéros **hors** de la droite critique (Davenport–Heilbronn, 1936) ;
- l'hypothèse de Riemann équivaut à la **positivité de Weil** (critère de Weil) ;
- le **produit eulérien** (la multiplicativité) est l'ingrédient structurel qui distingue les deux situations.

Ce dépôt en propose une **illustration numérique** soignée et reproductible, ainsi qu'une **synthèse** qui délimite précisément ce que *toute* approche de HR doit accomplir, et exhibe un contre-exemple (la zêta d'Epstein) qui échoue exactement là où la théorie le prédit. C'est un **discriminant**, pas un théorème.

---

## Résultats (résumé)

**Note I — la positivité de Weil discrimine ζ de l'Epstein.**
- Localisation et raffinement d'un **zéro hors-droite** de $Z_{Q_0}$ :
  $\rho_0 = 0{,}9532604747946\ldots + 16{,}2902157203903\ldots\, i$, avec $|Z_{Q_0}(\rho_0)|=10^{-81}$, reproduisant le tableau de Davenport–Heilbronn (zéros sur la droite *et* hors-droite).
- Validation de la **positivité de Weil pour ζ** ($\mathcal W(\zeta)=+0{,}037>0$, formule explicite vérifiée à $10^{-32}$).
- Décomposition en fonctions $L$ : $Z_{Q_0}=\tfrac23\zeta_K+\tfrac43 L(\cdot,f)$, et contraste avec l'objet à produit eulérien $\zeta_K$.

**Note II — la positivité de Weil place par place ; la multiplicativité comme ingrédient porteur.**
- Un **certificat explicite** $\mathcal W^{Q_0}(h)<0$ pour une fonction test admissible.
- **Reconstruction place par place** de la forme quadratique de Weil depuis le côté *arithmétique* (coefficients de $-Z'/Z$ + terme archimédien), **sans utiliser les zéros**, validée contre le côté zéros à $10^{-10}$ : $\lambda_{\min}(\zeta)\approx 0$ ($\ge 0$) contre $\lambda_{\min}(\text{Epstein})=-0{,}70$ ($<0$).
- **Isolation de la multiplicativité** : à pôle et facteur archimédien *identiques* (Epstein vs $\zeta_K$, même corps, même $\Gamma_{\mathbb C}$), la positivité bascule uniquement avec la nature des coefficients ; la négativité vit *entièrement* dans la masse de $-Z'/Z$ aux entiers composés (le défaut de multiplicativité).

---

## Contenu

| Fichier | Rôle |
|---|---|
| `note_weil_epstein.md` | Note de recherche I (zéro hors-droite, positivité de Weil pour ζ, décomposition $L$) |
| `note_weil_multiplicativite.md` | Note de recherche II (certificat explicite, reconstruction place par place, multiplicativité) |
| `ftsa_weil.py` | Banc d'essai principal : zêta d'Epstein (forme à gamma incomplète), recherche/raffinement de zéros, positivité de Weil pour ζ, spectre de la forme de Weil |
| `weil_local.py` | Reconstruction place par place de la forme de Weil depuis l'arithmétique + validation contre le côté zéros |
| `weil_decomp.py` | Isolation de la multiplicativité (Epstein vs $\zeta_K$) et décomposition de $\lambda_{\min}$ |

---

## Installation et usage

```bash
pip install mpmath numpy
```

```bash
# zero hors-droite, positivite de Weil pour zeta, decomposition L
python ftsa_weil.py polish 0.9533 16.290 90      # rho_0 a 10^-81, quadruplet
python ftsa_weil.py weil-zeta gauss 0.02         # W(zeta) = +0.037
python ftsa_weil.py weil-epstein                 # certificat W^Q0(h) < 0
python ftsa_weil.py entropy-flow                 # spectre (cote zeros)

# reconstruction place par place + validation
python weil_local.py

# isolation de la multiplicativite
python weil_decomp.py
```

`weil_local.py` et `weil_decomp.py` requièrent `ftsa_weil.py` dans le même répertoire. Le recensement des zéros sur la droite est mis en cache dans `onz.json`.

---

## Établi vs spéculatif

- **Établi et reproductible** : tous les calculs numériques ci-dessus (zéro hors-droite, positivité de Weil, reconstruction place par place, isolation de la multiplicativité). Ils illustrent des faits classiques.
- **Spéculatif** : le cadre « FTSA » (modèle de fonction de partition d'un superfluide adélique) et le programme géométrique Connes–Consani $\mathbb F_1$ / flot d'entropie à la Perelman — clairement signalés comme tels dans les notes. Aucune de ces directions n'est ici démontrée.

---

## Références

- H. Davenport, H. Heilbronn. *On the zeros of certain Dirichlet series*. J. London Math. Soc. **11** (1936).
- A. Weil. *Sur les « formules explicites » de la théorie des nombres premiers* (1952) — critère de positivité.
- A. Connes. *Trace formula in noncommutative geometry and the zeros of the Riemann zeta function*. Selecta Math. **5** (1999).
- B. Rodgers, T. Tao. *The de Bruijn–Newman constant is non-negative*. Forum Math. Pi **8** (2020).
- P. Sarnak, et al. — sur le rôle du produit eulérien et l'analogie des corps de fonctions (contexte).

---

## Licence

- **Code** (`*.py`) : licence **MIT** — voir `LICENSE`.
- **Notes de recherche** (`*.md`) : **Creative Commons Attribution 4.0 (CC BY 4.0)** — voir `LICENSE-CC-BY-4.0.md`.

En cas de réutilisation des notes, merci de citer ce dépôt.

---

## Remerciements et transparence sur l'assistance par IA

La **conception, la direction scientifique et l'ensemble des décisions de recherche** de ce projet — le programme, les questions posées, le choix du contre-exemple (zêta d'Epstein), le cadrage de la positivité de Weil comme discriminant, et l'interprétation des résultats — sont le fait de l'auteur.

La **mise en œuvre** a été réalisée avec l'assistance de **Claude** (un modèle de langage d'Anthropic) : écriture et exécution du code Python, calculs numériques, dérivation des expressions de la formule explicite, et rédaction des notes — sous la direction de l'auteur et d'après ses décisions.

Tous les résultats numériques ont été **reproduits et vérifiés indépendamment par l'auteur** (chaque script a été exécuté). La responsabilité du contenu, y compris d'éventuelles erreurs, incombe à l'auteur.
