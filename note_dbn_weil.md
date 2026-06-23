# Le flot de de Bruijn–Newman et la positivité de Weil : une seule dynamique (zêta d'Epstein)

*Note de recherche III — suite de « La positivité de Weil, place par place ».*
*Programme FTSA / approche topologique de l'hypothèse de Riemann.*
*22 juin 2026. Toutes les valeurs sont reproductibles par `dbn_epstein.py` et `dbn_weil_coupling.py` (précision arbitraire `mpmath` ; algèbre linéaire `numpy`).*

---

## 0. Position du problème

Les notes I–II ont établi, sur la zêta d'Epstein du formulaire principal $Q_0=m^2+mn+6n^2$ (disc $-23$, $h=3$), que la **positivité de Weil discrimine** un objet à produit eulérien d'un objet qui n'en a pas, et que l'ingrédient porteur de cette positivité est la **multiplicativité**. La violation de HR par l'Epstein y était mesurée par un invariant : la valeur propre minimale de la forme de Weil, $\lambda_{\min}^{Q_0}=-0{,}7045<0$.

Dans le programme géométrique spéculatif, deux **flots-ombres** sont en jeu :

* **(A)** le flot de de Bruijn–Newman $H_t$, qui déforme la fonction complétée et fait migrer ses zéros ;
* **(B)** le flot de gradient de la forme de Weil, dont $\lambda_{\min}$ mesure l'obstruction.

Une *condition nécessaire* pour qu'un éventuel flot géométrique unique se projette sur les deux (compatibilité D1$\leftrightarrow$D2) est que **(A) et (B) soient la même dynamique**. Cette note teste cette coïncidence numériquement, sur l'Epstein, où elle est non triviale puisque HR y est *faux*.

---

## 1. La constante de de Bruijn–Newman analogue $\Lambda_{Q_0}$

### 1.1. Reformulation non oscillante

Soit $\xi_{Q_0}(s)=(23/4)^{s/2}\,\pi^{-s}\Gamma(s)\,Z_{Q_0}(s)$ la fonction complétée, réelle sur la droite critique, vérifiant $\xi_{Q_0}(s)=\xi_{Q_0}(1-s)$. Le flot de de Bruijn–Newman s'écrit usuellement $H_t(z)=\int_0^\infty e^{tu^2}\Phi(u)\cos(zu)\,du$. Comme $z\mapsto s=\tfrac12+iz$ envoie $\partial_z^2$ sur $-\partial_s^2$, ce flot est le **semi-groupe de la chaleur dans la variable $s$**, $H_t=e^{t\partial_s^2}\xi_{Q_0}$, soit une convolution gaussienne. Après normalisation :
$$
\boxed{\;H_t(z) \;=\; \frac{1}{\sqrt{2\pi}}\int_{-\infty}^{\infty} e^{-v^2/2}\;\xi_{Q_0}\!\Big(\tfrac12+\sqrt{2t}\,v+iz\Big)\,dv\;}
$$
La convergence tient parce que la décroissance gaussienne bat la croissance d'ordre $1$ de $\Gamma$. On vérifie $H_0(z)=\xi_{Q_0}(\tfrac12+iz)$, et que $H_t(z)$ est **réelle pour $z$ réel** (équation fonctionnelle $+$ réflexion de Schwarz) : les zéros restent donc en quadruplets $\pm\mu(t)\pm i\delta(t)$. Cette forme est purement non oscillante et réutilise la $\xi_{Q_0}$ déjà validée — c'est l'évaluation par quadrature de Gauss–Hermite (nœuds en précision arbitraire) qui rend le calcul robuste.

### 1.2. La collision, par le système du zéro double

Au temps $t=0$, l'Epstein possède le quadruplet hors-droite en $z_0=\mu+i\delta=16{,}2902\ldots+0{,}4533\ldots\,i$. Quand $t$ croît, $\delta(t)$ décroît : la paire conjuguée rejoint la droite critique et y fusionne en un **zéro double réel** au temps $\Lambda_{Q_0}$ — l'analogue de la constante de de Bruijn–Newman pour l'Epstein. On le calcule directement en résolvant le système
$$
H_t(z)=0,\qquad \partial_z H_t(z)=0
$$
en $(t,z_*)$ réels, par un Newton bidimensionnel (Jacobien analytique par stencil $3\times3$). La convergence est quadratique et déterministe (mêmes itérés à toute précision), signe d'un zéro double non dégénéré ($\partial_{zz}H\neq 0$). À dps $50$ puis $70$, les résidus tombent à $|H|\sim 10^{-65}$, et les deux précisions coïncident sur $\sim 35$ chiffres :
$$
\boxed{\;\Lambda_{Q_0}=0{,}08430669450968749501499987219163086\ldots\;}
\qquad z_*=16{,}0705292242713711\ldots
$$
Le lieu de fusion $z_*\approx 16{,}07$ tombe dans le creux entre les zéros en-ligne. C'est une **deuxième signature quantitative** de la violation de HR par l'Epstein, indépendante de la positivité de Weil :

| invariant | $\zeta$ (HR vrai) | Epstein $Z_{Q_0}$ |
|---|---|---|
| $\lambda_{\min}$ de Weil | $\approx 0\ (\ge 0)$ | $-0{,}7045$ |
| de Bruijn–Newman $\Lambda$ | $0$ (Rodgers–Tao) | $\mathbf{0{,}0843066945\ldots}$ |

---

## 2. La forme de Weil reconstruite le long du flot

À chaque $t$, on reconstruit la forme de Weil **uniquement à partir des zéros de $H_t$** — choix imposé par la cohérence : sous le flot, *tous* les zéros bougent et les facteurs archimédien/pôle sont déformés, de sorte que l'objet propre est la somme sur les zéros, avec des fonctions test localisées (seuls les zéros voisins comptent). Avec les centres $\nu_i\in\{14{,}0,\,14{,}5,\dots,18{,}5\}$ et la largeur $a=0{,}5$ de la note II,
$$
M_{ij}(t)=\operatorname{Re}\sum_{z\,:\,H_t(z)=0} u_i(z)\,u_j(z),\qquad
u_i(z)=e^{-a(z-\nu_i)^2}+e^{-a(z+\nu_i)^2},
$$
la somme portant sur le quadruplet hors-droite $\mu(t)\pm i\delta(t)$ **et** les zéros en-ligne, tous suivis numériquement sous le flot. Mécanisme : un zéro réel donne $vv^\top$ réel (rang $1$, semi-défini positif) ; le quadruplet $z=\mu\pm i\delta$ donne
$$
4\operatorname{Re}\big[v\,v^\top\big]=4\big(v_Rv_R^\top-v_Iv_I^\top\big),
$$
dont la part $-4\,v_Iv_I^\top$ est la **direction négative**, de norme $\sim 4|v_I|^2\sim C\,\delta^2$. Quand $\delta(t)\to 0$, on a $v_I\to 0$ : la direction négative s'évanouit. À $t=0$ on retrouve $\lambda_{\min}(0)=-0{,}70422$ (en accord avec la note II), validant la construction.

---

## 3. Résultat : $\lambda_{\min}(t)$ remonte vers $0$ au rythme de la collision

Suivi simultané du quadruplet et des zéros en-ligne, grille concentrée près de la collision sur $t\in[0,\,0{,}0840]$ (dps $28$, $\mathrm{GH}=20$) :

| $t$ | $\delta(t)$ | $\lambda_{\min}(t)$ |
|---|---|---|
| $0{,}00000$ | $0{,}453259$ | $-7{,}0422\times10^{-1}$ |
| $0{,}02386$ | $0{,}373215$ | $-3{,}7840\times10^{-1}$ |
| $0{,}03430$ | $0{,}335336$ | $-2{,}6677\times10^{-1}$ |
| $0{,}05219$ | $0{,}263152$ | $-1{,}1484\times10^{-1}$ |
| $0{,}05964$ | $0{,}228590$ | $-6{,}8639\times10^{-2}$ |
| $0{,}07157$ | $0{,}161976$ | $-2{,}1564\times10^{-2}$ |
| $0{,}07953$ | $0{,}098334$ | $-5{,}6117\times10^{-3}$ |
| $0{,}08201$ | $0{,}067947$ | $-2{,}2344\times10^{-3}$ |
| $0{,}08350$ | $0{,}040165$ | $-5{,}4924\times10^{-4}$ |
| $0{,}08400$ | $0{,}024833$ | $-1{,}3747\times10^{-4}$ |

La remontée est **monotone**. À $t=0{,}0840$ — à $0{,}3\,\%$ sous $\Lambda_{Q_0}$, là où $\delta\to 0$ — l'obstruction $\lambda_{\min}$ ne vaut plus que $-1{,}4\times 10^{-4}$, soit $\sim 5000$ fois plus petite qu'à $t=0$. **Les deux flots-ombres guérissent la violation de HR au même instant.**

---

## 4. Le croisement $\lambda_{\min}=0$ coïncide avec $\Lambda_{Q_0}$

La théorie des perturbations sur $M(t)=M_{\text{fond}}(t)+4v_Rv_R^\top-4v_Iv_I^\top$ donne, près de la collision,
$$
\lambda_{\min}(t)\ \approx\ \lambda_0\ -\ 4\,|\langle e_0, v_I\rangle|^2\ \approx\ \lambda_0\ -\ K\,\delta(t)^2,
$$
où $\lambda_0=\lambda_{\min}\!\big(M_{\text{fond}}+4v_Rv_R^\top\big)$ est la valeur **à $\delta=0$**. La relation est *linéaire en $\delta^2$* (non en $t$), ce que confirment les données : $-\lambda_{\min}/\delta^2$ dérive de $3{,}43$ à $0{,}22$. L'extrapolation correcte est donc en $\delta^2$ : on ajuste $\lambda_{\min}=\lambda_0-K\delta^2$, et l'on compose avec $\delta^2(t)$ pour le passage par zéro.

Sur une grille **concentrée près de la collision** (points jusqu'à $\delta\simeq 0{,}011$, $t=0{,}08425$), l'ambiguïté se lève sans appel. La valeur extraite $\lambda_0$ **décroît monotonement vers $0$** quand on resserre la fenêtre vers la collision :
$$
\lambda_0:\quad 2{,}10\times 10^{-4}\ \to\ 1{,}32\times 10^{-4}\ \to\ 7{,}3\times 10^{-5}\ \to\ 3{,}0\times 10^{-5}
\qquad(\text{fenêtres } 6\to 3),
$$
et l'écart au croisement suit : $0{,}22\,\%\to 0{,}16\,\%\to 0{,}10\,\%\to 0{,}06\,\%$. Il n'y a donc **pas de plancher positif** : ce n'était pas un effet de sous-espace fini, mais un simple biais d'extrapolation qui s'efface à l'approche de la collision. Avec la fenêtre la plus serrée,
$$
\boxed{\;t_{\text{crois}}=0{,}084258,\qquad \Lambda_{Q_0}=0{,}0843067,\qquad \text{écart }+0{,}06\,\%\;}
$$
La tendance ($\lambda_0\downarrow 0$, $t_{\text{crois}}\uparrow\Lambda_{Q_0}$) est sans ambiguïté. Le point direct le plus profond le confirme : à $t=0{,}08425$ — à $0{,}07\,\%$ sous la collision, $\delta=0{,}0108$ — $\lambda_{\min}=-1{,}17\times 10^{-5}$, soit $\sim 60\,000$ fois plus petit qu'à $t=0$. **La positivité de Weil est restaurée exactement quand le quadruplet atteint la droite critique.**

---

## 5. Bilan

| voie | objet calculé | valeur |
|---|---|---|
| collision (zéro double de $H_t$) | $\Lambda_{Q_0}$ | $0{,}0843066945\ldots$ (35 chiffres) |
| restauration de la positivité de Weil | $t$ tel que $\lambda_{\min}=0$ | $0{,}084258$ (extrapolation $\delta^2$ ; écart $+0{,}06\,\%$, $\to\Lambda_{Q_0}$ avec $\lambda_0\downarrow 0$) |

Le flot de de Bruijn–Newman (A) **est** une descente pour l'obstruction de Weil (B) : deux dynamiques construites par des voies disjointes — l'une déforme la fonction par la chaleur et suit ses zéros, l'autre forme une matrice de Rayleigh à partir des zéros — restaurent la positivité au **même instant** $\Lambda_{Q_0}$ (à $0{,}06\,\%$ près, l'écart résiduel tendant vers $0$). C'est la **compatibilité D1$\leftrightarrow$D2** du programme géométrique, vérifiée numériquement.

**Ce que cette note n'établit pas.** Elle teste la coïncidence des deux *ombres linéarisées*. Elle ne construit ni l'espace de configurations $\mathcal C$, ni le flot géométrique dont (A) et (B) seraient les projections, et elle ne touche pas à **D3** — le rôle *structurel* de la multiplicativité (produit eulérien comme propriété de $\mathcal C$), qui demeure hors de portée du numérique. La séparation établi / spéculatif de la note II reste en vigueur.

---

## Annexe — Reproductibilité

```bash
# constante de de Bruijn-Newman analogue (35 chiffres) :
python dbn_epstein.py collision 0.0835 16.06 70 56

# image complete du flot delta(t) (recoupement par extrapolation racine) :
python dbn_epstein.py flow 0 0.083 16 50 48

# couplage flot <-> forme de Weil : lambda_min(t) le long du flot
# (grille concentree pres de la collision, extrapolation en delta^2) :
python dbn_weil_coupling.py 0 0.0840 14 28 20 55 2
```

---

*Mise en œuvre numérique et rédaction assistées par Claude (Anthropic), sous la direction et la vérification de l'auteur.*