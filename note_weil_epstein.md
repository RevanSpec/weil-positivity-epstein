# La positivité de Weil discrimine ζ de l'Epstein : deux expériences numériques

*Note de recherche — programme FTSA / approche topologique de l'hypothèse de Riemann.*
*16 juin 2026. Toutes les valeurs sont calculées en précision arbitraire (mpmath) ; le banc d'essai `ftsa_weil.py` reproduit chaque chiffre.*

---

## 1. La question

Le fil conducteur de nos échanges a isolé, étape après étape, **la** maladie structurelle des constructions de type « gaz de primons » / FTSA : elles s'appuient sur une positivité **gratuite** (positivité d'Osterwalder–Schrader, équivalente à l'auto-adjonction du hamiltonien libre, ou plus prosaïquement à la positivité des coefficients de Dirichlet $r(k)\ge 0$), alors que la condition réellement ouverte est la **positivité de Weil** :
$$
\text{HR} \;\Longleftrightarrow\; \mathcal W(g)=\sum_{\rho} h(\gamma_\rho)\;\ge\;0
\quad\text{pour toute } g \text{ admissible},
$$
la somme portant sur les zéros non triviaux $\rho=\tfrac12+i\gamma_\rho$, avec $h(r)=\hat g(r)\ge 0$ sur $\mathbb R$ (fonction de type positif, $g=\varphi\star\tilde\varphi$).

Notre hypothèse de travail : **la positivité de Weil discrimine** un objet à produit eulérien (zéros sur la droite) d'un objet qui n'en a pas (zéros hors-droite), et **le produit eulérien — la multiplicativité — est l'ingrédient opératoire**. Les deux expériences ci-dessous la confirment de façon décisive.

---

## 2. Expérience 1 — la positivité de Weil pour ζ (validée à 32 chiffres)

On teste la formule explicite de Riemann–Weil avec la fonction test gaussienne $h(r)=e^{-a r^2}$ ($a=0{,}02$), de type positif, dont la transformée $g(u)=\tfrac{1}{2\sqrt{\pi a}}e^{-u^2/4a}\ge 0$.

**Côté zéros** (somme sur les zéros, $\gamma$ via les premiers zéros de ζ) :
$$
\mathcal W(\zeta)=\sum_\rho h(\gamma_\rho)=2\sum_{\gamma>0}e^{-a\gamma^2}= 0{,}0370825834626.
$$

**Côté arithmétique** (formule explicite) :
$$
\mathcal W(\zeta)=\underbrace{2h(\tfrac i2)}_{\text{pôle}}\;-\;\underbrace{g(0)\log\pi}_{\text{archi. constant}}\;+\;\underbrace{\frac{1}{2\pi}\int_{\mathbb R} h(r)\,\operatorname{Re}\psi\!\Big(\tfrac14+\tfrac{ir}{2}\Big)\,dr}_{\text{archimédien}}\;-\;\underbrace{2\sum_{n\ge 2}\frac{\Lambda(n)}{\sqrt n}\,g(\log n)}_{\text{premiers}}.
$$

Décomposition numérique :

| terme | valeur |
|---|---|
| pôle $2h(i/2)$ | $+2{,}010025$ |
| $-g(0)\log\pi$ | $-2{,}283406$ |
| archimédien $\frac1{2\pi}\!\int h\,\operatorname{Re}\psi(\tfrac14+\tfrac{ir}2)$ | $+0{,}315283$ |
| premiers | $-0{,}004820$ |
| **total** | $\mathbf{+0{,}0370825834626}$ |

**Écart entre les deux côtés : $4\times 10^{-32}$.** La formule explicite est donc exacte à la précision de travail, et
$$
\boxed{\;\mathcal W(\zeta) = +0{,}037 > 0\;}\qquad\text{(compatible avec HR).}
$$

Pour ζ, les quatre contributions s'équilibrent en un nombre positif ; en particulier le terme **archimédien** $\operatorname{Re}\psi(\tfrac14+\tfrac{ir}{2})\sim \log|r|$ — *absent de toute construction de type gaz de primons* — pèse $+0{,}315$, du même ordre que les autres.

---

## 3. Expérience 2 — un zéro hors-droite de la zêta d'Epstein

### 3.1 L'objet

Soit $Q_0(m,n)=m^2+mn+6n^2$, le **formulaire principal** de discriminant $-23$ (nombre de classes $h=3$). Sa zêta d'Epstein
$$
Z_{Q_0}(s)=\sum_{(m,n)\neq(0,0)} Q_0(m,n)^{-s}=\sum_{k\ge 1} r(k)\,k^{-s},\qquad r(k)=\#\{(m,n):Q_0=k\}\ge 0,
$$
a des **coefficients de Dirichlet positifs** (c'est, mot pour mot, la « positivité physique » de type OS sur laquelle s'appuie le FTSA). Elle possède un prolongement méromorphe, un pôle simple en $s=1$ de résidu $2\pi/\sqrt{23}=1{,}3101347$ (vérifié numériquement : $1{,}3101354$), et une équation fonctionnelle $s\leftrightarrow 1-s$ via la transformation thêta (Poisson). **Mais elle n'a aucun produit eulérien** : $Q_0$ représente une seule classe d'idéaux, pas le groupe de classes entier.

L'évaluation utilise la représentation à gamma incomplète (exacte, dérivée de la transformation thêta) :
$$
\pi^{-s}\Gamma(s)Z_{Q_0}(s)=\sum_k r(k)(\pi k)^{-s}\Gamma(s,\pi k)+\frac{1}{\sqrt{\det A}}\sum_k r(k)\,\beta_k^{\,s-1}\Gamma(1-s,\beta_k)+\frac{1}{\sqrt{\det A}\,(s-1)}-\frac1s,
$$
avec $A=\begin{psmallmatrix}1&1/2\\1/2&6\end{psmallmatrix}$, $\det A=23/4$, $\beta_k=\tfrac{4\pi}{23}k$. (Vérifiée contre la somme directe pour $\operatorname{Re}s>1$.)

### 3.2 Le zéro hors-droite

Une recherche fiable (précision adaptée $\approx 1{,}3\,t$ chiffres, $\mathrm{KMAX}\gtrsim 0{,}42\,\mathrm{dps}$) puis un raffinement de Newton à dps = 90, KMAX = 500 donnent un zéro **non trivial hors de la droite critique** :
$$
\boxed{\;\rho_0 = 0{,}9532604747946606862505\ldots \;+\; 16{,}2902157203903907929\ldots\, i\;}
$$
$$
\operatorname{Re}(\rho_0)-\tfrac12 = 0{,}4532604\ldots \neq 0,\qquad |Z_{Q_0}(\rho_0)| = 10^{-81}.
$$
Coefficients réels et équation fonctionnelle obligent : le **quadruplet** complet est constitué de zéros (chacun $|Z|\approx 10^{-81}$) :
$$
\{\rho_0,\ \bar\rho_0,\ 1-\rho_0,\ 1-\bar\rho_0\}=\{\,0{,}9533\pm 16{,}290\,i,\ \ 0{,}0467\pm 16{,}290\,i\,\}.
$$

### 3.3 Le contexte (tableau de Davenport–Heilbronn)

$Z_{Q_0}$ a aussi de **nombreux zéros sur la droite** $\sigma=\tfrac12$ (calculés ici, $\operatorname{Re}=\tfrac12$ exactement) :
$$
t \approx 4{,}942,\ \ 7{,}015,\ \ 10{,}494,\ \ 11{,}804,\ \ 18{,}828,\ \ 19{,}304,\ \dots
$$
Le quadruplet hors-droite ($t\approx\pm 16{,}29$) loge **précisément dans l'intervalle** entre les zéros sur-droite $t=11{,}804$ et $t=18{,}828$. C'est exactement le phénomène de **Davenport–Heilbronn (1936)** : une zêta d'Epstein de nombre de classes $>1$ a une proportion positive de zéros sur la droite *et* une infinité de zéros hors-droite (dont une infinité dans $\sigma>1$ ; celui que nous localisons est dans la bande $\tfrac12<\sigma<1$).

### 3.4 Le mécanisme : produit eulérien vs combinaison linéaire

Par les caractères du groupe de classes $\mathrm{Cl}=\mathbb Z/3\mathbb Z$ (avec $w=2$, $h=3$) :
$$
Z_{Q_0}(s)=\frac{w}{h}\sum_{\psi\in\widehat{\mathrm{Cl}}}L(s,\psi)=\frac23\,\zeta_K(s)+\frac43\,L(s,f),
$$
où $\zeta_K(s)=\zeta(s)\,L(s,\chi_{-23})$ est la zêta de Dedekind de $\mathbb Q(\sqrt{-23})$ (caractère trivial) et $L(s,f)$, $f=\eta(z)\eta(23z)$ la forme primitive de poids 1, niveau 23 (les deux caractères cubiques conjugués).

*Identité vérifiée numériquement* à $s=2$ : $Z_{Q_0}(2)=2{,}4385325$ contre $\tfrac23\zeta_K(2)+\tfrac43 L(2,f)=2{,}4385485$ (écart $1{,}6\times 10^{-5}$, dû à la troncature de la série de $f$).

- $\zeta_K$ et $L(\cdot,f)$ ont **chacune un produit eulérien** $\Rightarrow$ leurs zéros sont sur la droite (HR vérifiée pour les premiers).
- $Z_{Q_0}$ est leur **combinaison linéaire** $\tfrac23\zeta_K+\tfrac43 L(\cdot,f)$ : *une somme de produits eulériens n'est pas un produit eulérien*, et c'est de là que naissent les zéros hors-droite.

**Contraste direct** : à $\rho_0$, l'objet à produit eulérien construit sur le *même* corps est régulier,
$$
|\zeta_K(\rho_0)| = 1{,}456 \neq 0,
$$
et $|\zeta_K(\tfrac12+it)|$ plonge (zéros) près de $t=14$ et $t=19$ — *sur la droite* — sans jamais plonger à $\sigma=0{,}953$. Le retrait du produit eulérien (passage de $\zeta_K$ au seul morceau $Z_{Q_0}$) **déplace des zéros hors de la droite**.

---

## 4. Conclusion et conséquence pour le FTSA

Le critère de Weil tranche sans ambiguïté :

1. **ζ (produit eulérien)** : $\mathcal W(\zeta)=+0{,}037>0$, formule explicite exacte à $10^{-32}$ — zéros sur la droite.
2. **$Z_{Q_0}$ (pas de produit eulérien, mais $r(k)\ge 0$)** : zéro hors-droite $\rho_0$ à $10^{-81}$. Par le **théorème de Weil** ($\mathcal W(g)\ge 0\ \forall g \Leftrightarrow$ tous les zéros sur la droite), l'existence de $\rho_0$ équivaut à $\mathcal W^{Q_0}(g)<0$ pour une certaine $g$ admissible : **la positivité de Weil échoue**.

> **L'ironie est la leçon.** $Z_{Q_0}$ possède exactement la positivité dont se prévaut le FTSA — coefficients $r(k)\ge0$, fonction de partition « physique », auto-adjonction du gaz sous-jacent — et **viole pourtant HR**. La positivité OS / des coefficients est **nécessaire mais notoirement insuffisante**.

**Ce qu'un objet FTSA devrait reproduire, ce n'est pas la positivité OS, mais la positivité de Weil.** Concrètement, cela impose deux ingrédients qu'aucune construction de type gaz de primons ne contient :

- **le produit eulérien** (multiplicativité) — opératoirement responsable de la localisation des zéros, comme le montre le passage $Z_{Q_0}\to\zeta_K$ ;
- **le terme archimédien** $\operatorname{Re}\psi(\tfrac14+\tfrac{ir}{2})$ — dominant dans la formule explicite ($+0{,}315$ ci-dessus), absent du gaz de primons, et qui encode la place à l'infini.

C'est exactement la cible du programme géométrique (Connes–Consani $\mathbb F_1$, flot d'entropie à la Perelman) évoqué en ligne de mire : une fonctionnelle dont le point critique coïnciderait avec la **positivité de Weil–Connes** — et non avec une positivité spectrale gratuite.

---

## Annexe — le banc d'essai

`ftsa_weil.py` (autonome, `pip install mpmath`) encode toutes les formules ci-dessus et se fait évoluer au fil des résultats :

```
python ftsa_weil.py selftest                  # representation integrale vs somme directe
python ftsa_weil.py online 2 20               # zeros de Z_Q0 sur sigma=1/2
python ftsa_weil.py hunt 0.46 1.10 14 18 0.03 0.2   # chasse + raffinement des zeros
python ftsa_weil.py verify 0.9533 16.290      # confirmation d'un zero (precision croissante)
python ftsa_weil.py weil-zeta gauss 0.02      # positivite de Weil pour zeta
```

**Pistes immédiates pour la suite** (à coder dans le même fichier) :
- exhiber une valeur **négative explicite** $\mathcal W^{Q_0}(g)<0$ (construction de $g$ de type positif concentrée près de $\gamma_0$, à la Bombieri — délicate car un simple pic à $\pm 16{,}29$ n'est pas de type positif) ;
- principe de l'argument / comptage des zéros hors-droite de $Z_{Q_0}$ par hauteur ;
- formule explicite d'Epstein (facteur $\Gamma(s)$, côté « coefficients » $-Z'/Z$ sans von Mangoldt) pour retrouver $\mathcal W^{Q_0}<0$ par le côté arithmétique ;
- transposition au programme $\mathbb F_1$ / flot d'entropie.


---

*Mise en œuvre numérique et rédaction assistées par Claude (Anthropic), sous la direction et la vérification de l'auteur.*
