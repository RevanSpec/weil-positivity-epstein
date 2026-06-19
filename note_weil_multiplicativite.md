# La positivité de Weil, place par place : la multiplicativité comme ingrédient porteur

*Note de recherche II — suite de « La positivité de Weil discrimine ζ de l'Epstein ».*
*Programme FTSA / approche topologique de l'hypothèse de Riemann.*
*17 juin 2026. Toutes les valeurs sont reproductibles par `ftsa_weil.py`, `weil_local.py`, `weil_decomp.py` (précision arbitraire `mpmath` ; algèbre linéaire `numpy`).*

---

## 0. Rappel et question

La note I a établi, par deux expériences numériques, que la **positivité de Weil discrimine** un objet à produit eulérien d'un objet qui n'en a pas :
$$
\text{HR} \iff \mathcal W(g)=\sum_\rho h(\gamma_\rho)\ge 0 \ \text{ pour toute } g \text{ admissible } (h=\hat g\ge 0 \text{ sur }\mathbb R),
$$
avec $\mathcal W(\zeta)=+0{,}037>0$ d'un côté, et de l'autre un **zéro hors-droite** de la zêta d'Epstein du formulaire principal $Q_0=m^2+mn+6n^2$ (disc $-23$, $h=3$),
$$
\rho_0 = 0{,}953260474794660686250509\ldots + 16{,}290215720390390792963\ldots\, i,\qquad |Z_{Q_0}(\rho_0)|=10^{-81}.
$$

Mais deux insatisfactions subsistaient. **(i)** L'échec de la positivité de Weil pour l'Epstein n'était établi qu'*indirectement*, par le théorème de Weil ($\rho_0$ hors-droite $\Rightarrow \exists\, g$ avec $\mathcal W^{Q_0}(g)<0$) ; on n'avait pas exhibé de certificat explicite. **(ii)** La positivité était lue *sur les zéros* — ce qui, pour tester HR, est circulaire. Cette note II répond aux deux : on exhibe le certificat, on **reconstruit la fonctionnelle de Weil place par place** (sans jamais regarder les zéros), et on **localise l'origine** de l'échec. Cette origine est la **multiplicativité**.

---

## 1. Un certificat explicite $\mathcal W^{Q_0}(h)<0$

Par l'équation fonctionnelle et la conjugaison, $\rho_0$ vient en quadruplet $\{\rho_0,\bar\rho_0,1-\rho_0,1-\bar\rho_0\}$, soit en coordonnées $\gamma=(\rho-\tfrac12)/i$ les quatre points
$$
\pm\mu \pm i\delta,\qquad \mu=\operatorname{Im}\rho_0=16{,}290,\quad \delta=\operatorname{Re}\rho_0-\tfrac12=0{,}4533 .
$$
On choisit une fonction test **admissible** $h=k^2\ge 0$ sur $\mathbb R$, avec $k$ réelle paire possédant un **nœud en $\pm\mu$** et de largeur $w$ :
$$
k(r)=(r-\mu)\,e^{-(r-\mu)^2/2w^2}-(r+\mu)\,e^{-(r+\mu)^2/2w^2}.
$$
Le quadruplet hors-droite contribue alors $4\operatorname{Re}\!\big[h(\mu+i\delta)\big]$, négatif (car $k(\mu+i\delta)\approx i\delta\,e^{\delta^2/2w^2}$ est quasi imaginaire pur, donc $k^2<0$), tandis que les zéros **sur la droite** — tous à distance $\gtrsim 1{,}6$ de $\mu$ — sont filtrés. La somme $\mathcal W^{Q_0}(h)=\sum_\rho h(\gamma_\rho)$ donne, pour toute largeur testée :

| $w$ | sur-droite $(\ge 0)$ | hors-droite | $\mathcal W^{Q_0}(h)$ |
|---|---|---|---|
| 0,4 | $4\times 10^{-17}$ | $-2{,}968$ | $\mathbf{-2{,}968}$ |
| 0,6 | $2\times 10^{-7}$ | $-1{,}454$ | $\mathbf{-1{,}454}$ |
| 0,8 | $5{,}6\times 10^{-4}$ | $-1{,}132$ | $\mathbf{-1{,}132}$ |
| 1,0 | $0{,}0226$ | $-1{,}009$ | $\mathbf{-0{,}987}$ |

Plus le filtre est étroit, plus la part sur-droite s'évanouit et plus la contribution négative du quadruplet domine seule. C'est le **certificat explicite** que la positivité de Weil échoue pour l'Epstein — rigoureux à la qualité près de notre recensement des zéros près de $t=16{,}29$, qui est complète (aucun zéro sur la droite à moins de $1{,}6$ de $\mu$).

---

## 2. Le spectre de la forme quadratique de Weil

La fonctionnelle de Weil, restreinte à $g=\varphi\star\varphi^\ast$ (de sorte que $h=|\hat\varphi|^2\ge 0$), est une **forme quadratique hermitienne**. Sur une base réelle paire $\{\hat\varphi_i\}$,
$$
\mathcal W(\varphi\star\varphi^\ast)=\sum_{i,j}c_i\bar c_j\,M_{ij},\qquad M_{ij}=\sum_\rho \hat\varphi_i(\gamma_\rho)\,\hat\varphi_j(\gamma_\rho),
$$
et **HR $\iff M\succeq 0$**. Pour $\zeta$ (zéros réels) $M$ est une matrice de Gram, donc $\lambda_{\min}\ge 0$ ; pour l'Epstein, le quadruplet hors-droite (points complexes $\mu\pm i\delta$) rend $M$ indéfinie. Une base de 10 gaussiennes centrées près de $\mu$ donne
$$
\lambda_{\min}(\zeta)\approx 0\ (\ge 0)\qquad\text{contre}\qquad \lambda_{\min}(\text{Epstein})<0 .
$$
Le **flot d'entropie** — descente de gradient sur le quotient de Rayleigh $R(c)=c^\top M c/\|c\|^2$, monotone décroissante — converge vers $\lambda_{\min}$ (le vecteur propre minimal *est* le certificat optimal du §1). Ici $M$ est encore bâtie *à partir des zéros* : ce paragraphe ne fait que rendre manifeste l'équivalence « zéro hors-droite $\Leftrightarrow$ valeur propre négative ». Le test réellement non circulaire est le suivant.

---

## 3. Reconstruction place par place (sans les zéros)

On reconstruit $\mathcal W$ depuis le **côté arithmétique** de la formule explicite, où les zéros n'apparaissent jamais. Pour une fonction $L$ de fonction complétée $\Lambda(s)=Q^{s/2}\gamma_\infty(s)L(s)$ :
$$
\mathcal W(g)=\underbrace{2h(\tfrac i2)}_{\text{pôle}}+\underbrace{\frac{1}{2\pi}\int_{\mathbb R} h(r)\,2\operatorname{Re}\frac{G'}{G}\!\big(\tfrac12+ir\big)\,dr}_{\text{place archimédienne}}-\underbrace{2\sum_{k\ge 2}\frac{c(k)}{\sqrt k}\,g(\log k)}_{\text{places finies}},\quad G=Q^{s/2}\gamma_\infty .
$$

- **$\zeta$** : $c(n)=\Lambda(n)$ (von Mangoldt, supportée sur les puissances de premiers), archimédien $-\log\pi+\operatorname{Re}\psi(\tfrac14+\tfrac{ir}{2})$.
- **Epstein** : complétée $\Lambda(s)=\pi^{-s}\Gamma(s)Z_{Q_0}(s)$, conducteur $N=\det A=\tfrac{23}{4}$, d'où l'archimédien $\log\tfrac{23}{4}-2\log\pi+2\operatorname{Re}\psi(\tfrac12+ir)$ ; et coefficients de $-Z'/Z$,
$$
c(k)=\sum_{d\mid k} r(d)\,\log d\ \cdot\ b(k/d),\qquad b=\text{inverse de Dirichlet de } r .
$$
Ces $c(k)$ sont **non multiplicatifs** : $c(6)=3{,}58$, $c(12)=4{,}97$ (non nuls aux composés — la signature de l'absence de produit eulérien), l'identité exacte $\sum_{d\mid k}r(d)c(k/d)=r(k)\log k$ étant vérifiée à $1{,}4\times10^{-14}$.

**Validation.** Après recensement complet des zéros sur la droite (12 jusqu'à $T=22$, par changements de signe de $\xi(\tfrac12+it)$ réelle, dont un très bas à $t=1{,}309$ que les balayages grossiers rataient), le côté arithmétique reproduit le côté zéros à
$$
\sim 10^{-10}\ (\text{Epstein})\qquad \text{et}\qquad \sim 10^{-16}\ (\zeta) .
$$
La formule explicite de l'Epstein est donc exacte.

**Résultat.** La matrice de Weil reconstruite *uniquement* à partir de $c(k)$ et du terme archimédien donne
$$
\lambda_{\min}(\zeta)=-3{,}3\times10^{-13}\approx 0\ (\ge 0,\ \text{Weil OK}),\qquad \boxed{\lambda_{\min}(\text{Epstein})=-0{,}7045<0\ (\text{Weil ECHOUE})},
$$
en accord avec le côté zéros à **6 chiffres** (écart entrywise $5\times10^{-7}$), et le flot d'entropie y converge ($R:\ +2{,}30\to-0{,}7045$). **La positivité de Weil, calculée place par place, sépare $\zeta$ de l'Epstein.** (Le recensement complet corrige au passage la valeur : $\lambda_{\min}=-0{,}70$, et non $-1{,}35$ — ce dernier était un artefact de zéros manquants, qui ajoutaient des contributions positives.)

---

## 4. La multiplicativité isolée — l'expérience décisive

Reste à savoir *quel* ingrédient, place par place, injecte la négativité. On le tranche par une expérience parfaitement contrôlée. L'Epstein $Z_{Q_0}$ et la zêta de Dedekind $\zeta_K=\zeta\,L(\cdot,\chi_{-23})$ sont **tous deux de degré 2 sur le même corps** $\mathbb Q(\sqrt{-23})$. Ils partagent donc le **même pôle** et le **même facteur archimédien** $\Gamma_{\mathbb C}(s)$, par la coïncidence
$$
\log\tfrac{23}{4}-2\log\pi \;=\; \log 23 - 2\log 2\pi .
$$
La **seule** variable est la nature des coefficients de $-L'/L$ :
$$
\text{Epstein}:\ c(k)\ \text{(non multiplicatif)}\qquad\text{vs}\qquad \zeta_K:\ \Lambda(n)\big(1+\chi_{-23}(n)\big)\ \text{(Euler, puissances de premiers)} .
$$
À pôle et archimédien identiques, on obtient :

| forme de Weil (mêmes pôle et archimédien) | $\lambda_{\min}$ |
|---|---|
| Epstein, $c(k)$ **non multiplicatif** | $\mathbf{-0{,}7045}$  (Weil échoue) |
| $\zeta_K$, $\Lambda(n)(1+\chi)$, **Euler** | $-1{,}45\times10^{-7}\approx 0$  (Weil OK) |
| Epstein, **termes composés retirés** | $+2{,}78\times10^{-8}\ge 0$  (positivité **restaurée**) |

Et la **décomposition place par place** de $\lambda_{\min}(\text{Epstein})$, au vecteur propre optimal, est éloquente :
$$
\underbrace{0}_{\text{pôle}}\ \underbrace{+\,1{,}764}_{\text{archimédien}}\ \underbrace{-\,0{,}701}_{-\text{coeff. } p^m}\ \underbrace{-\,1{,}767}_{-\text{coeff. COMPOSÉS}}\ =\ -0{,}7045 .
$$

Trois faits décisifs en ressortent :

1. **La multiplicativité bascule le signe** : tout étant gelé par ailleurs, passer de coefficients Euler ($\zeta_K$) à des coefficients non multiplicatifs (Epstein) fait basculer $\lambda_{\min}$ de $\approx 0$ à $-0{,}70$. La positivité de Weil **est exactement ce que la multiplicativité fournit**.
2. **La négativité vit dans les termes composés** : retirer la masse de $-Z'/Z$ aux entiers composés (le *défaut de multiplicativité*, nul pour un produit eulérien) restaure $\lambda_{\min}\ge 0$. C'est le terme $-\!\sum_{k\ \text{composé}}c(k)k^{-1/2}(\cdots)=-1{,}77$ qui détruit, à lui seul, la positivité.
3. **L'archimédien est du côté positif** : le pôle est nul (le vecteur optimal a un nœud) et la place à l'infini contribue $+1{,}76$ — elle « pousse » vers la positivité, loin d'être coupable.

---

## 5. Bilan

On a démontré numériquement, au niveau de la forme quadratique de Weil et **place par place**, ce qui n'était au départ qu'un diagnostic structurel :

> **Le produit eulérien — la multiplicativité — est l'ingrédient porteur de la positivité de Weil.** Son absence s'incarne précisément dans la masse de $-L'/L$ aux entiers composés, et c'est elle *seule* qui détruit la positivité ; le pôle et la place archimédienne sont, eux, du côté de la positivité.

Le tableau complet, du plus faible au plus fort :

- **positivité OS / coefficients $r(k)\ge 0$** : nécessaire mais *insuffisante* (l'Epstein la possède et viole HR) ;
- **positivité de Weil** : équivalente à HR, et *discriminante* ;
- **multiplicativité** (produit eulérien) : l'ingrédient *opératoire*, isolé ci-dessus ;
- **place archimédienne** $\operatorname{Re}\psi(\tfrac12+ir)$ : présente, correcte, et *positive*.

**Conséquence pour toute approche de type FTSA.** La cible n'est pas de reproduire une positivité spectrale gratuite (OS), mais la **positivité de Weil** — c'est-à-dire d'**encoder structurellement le produit eulérien** (la multiplicativité), accompagné du terme archimédien. Le numérique a donné la cible exacte ; il ne peut pas, par lui-même, bâtir la géométrie qui l'atteindrait.

Il faut être clair sur le statut : tout ceci est un **discriminant**, pas une preuve. On a délimité, de façon tranchée et reproductible, ce que *toute* démonstration de HR doit faire — et exhibé un contre-exemple (l'Epstein) qui échoue exactement là où l'on prédit.

---

## Annexe A. Le programme géométrique (spéculatif)

La cible étant fixée, la suite est conceptuelle : construire l'objet géométrique (Connes–Consani $\mathbb F_1$ fusionné à un flot à la Perelman) dont une **fonctionnelle d'entropie monotone** aurait pour point critique cette positivité de Weil, la multiplicativité y étant encodée par structure. Trois sous-problèmes restent ouverts, et difficiles :

1. **L'espace et le flot** : définir précisément l'objet sur lequel couler (triplets spectraux, métriques sur l'espace des classes d'adèles, données $\mathbb F_1$) et le « flot de Ricci » associé. Le flot de de Bruijn–Newman $(H_t)$ — concret — en serait l'*ombre linéarisée*.
2. **L'entropie monotone** : l'analogue de la $\mathcal W$-entropie de Perelman, à monotonie démontrable, dont les points critiques soient les configurations Weil-positives.
3. **Point critique $\Leftrightarrow$ positivité de Weil**, avec la **multiplicativité encodée structurellement** (et non imposée) — sachant Rodgers–Tao $\Lambda_{\text{dBN}}\ge 0$ et HR $\iff \Lambda_{\text{dBN}}=0$.

Le risque dominant y est l'auto-illusion : bâtir un flot qui « marche » parce qu'on y aurait glissé la réponse. Toute avancée devra distinguer scrupuleusement l'établi du spéculatif.

## Annexe B. Reproductibilité

```
# le zero hors-droite, la positivite de Weil pour zeta, la decomposition L
python ftsa_weil.py polish 0.9533 16.290 90      # rho_0 a 10^-81, quadruplet
python ftsa_weil.py weil-zeta gauss 0.02         # W(zeta)=+0.037
python ftsa_weil.py weil-epstein                 # certificat W^Q0(h) < 0  (sec. 1)
python ftsa_weil.py entropy-flow                 # spectre cote zeros       (sec. 2)

# reconstruction place par place + validation + spectre arithmetique
python weil_local.py                             # sec. 3 (validation 10^-10, lambda_min)

# isolation de la multiplicativite
python weil_decomp.py                            # sec. 4 (Epstein vs zeta_K, decomposition)
```


---

*Mise en œuvre numérique et rédaction assistées par Claude (Anthropic), sous la direction et la vérification de l'auteur.*
