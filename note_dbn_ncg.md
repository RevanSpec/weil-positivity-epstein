# Note VI — Le flot de de Bruijn–Newman dans le cadre métaplectique de Connes–Consani–Moscovici : localisation et obstruction

*Suite de la Note V. Celle-ci avait situé le FTSA relativement à la frontière et identifié l'obstacle contraignant — le lemme analytique de positivité de trace, ouvert depuis Connes 1999 — en notant qu'un seul élément du programme restait structurellement absent du cadre de Connes–Consani : la déformation de de Bruijn–Newman $\phi_t$. Cette note consigne l'exploration de ce thread et son verdict honnête. Elle ne prétend à aucun levier sur HR ; elle localise précisément où et pourquoi le flot dBN rejoint le cadre NCG, et où le calcul fidèle bute.*

---

## 0. La question

Le flot dBN déforme $\xi\to\xi_t$ ; ses zéros redeviennent réels pour $t\ge\Lambda$, et HR $\iff\Lambda=0$ (Rodgers–Tao : $\Lambda\ge0$). C'est le seul objet du FTSA absent du cadre de Connes–Consani–Moscovici (CCM), qui repose sur le scaling fixe ($=$ Frobenius en caractéristique zéro) et l'opérateur prolate. Question : le flot dBN a-t-il un sens dans ce cadre, et sa monotonie (Note III : $\lambda_{\min}$ de la forme de Weil $\nearrow0$ à la collision) éclaire-t-elle la positivité ?

---

## 1. Le cadre métaplectique de CCM

CCM réalisent $\mathfrak{sl}(2,\mathbb R)$ sur $L^2(\mathbb R)$ (base de Hermite) :
$$
\varpi(h)=x\partial_x+\tfrac12\ (=i\,\mathcal S,\ \text{scaling}),\quad
\varpi(e_+)=i\pi x^2,\quad
\varpi(e_-)=\tfrac{i}{4\pi}\partial_x^2,\quad
\varpi(k)=\pi x^2-\tfrac1{4\pi}\partial_x^2\ (\text{Hermite}).
$$
En forme « paire cyclique » $(D,\xi)$ — $D$ multiplication par $s$, mesure $\propto|\Gamma(\tfrac14+\tfrac{i}{2}s)|^2$, matrice de Jacobi $a_n=\tfrac12\sqrt{(2n+1)(2n+2)}$, nombre $N$ — l'**opérateur prolate** est
$$
\boxed{\ \omega(D,\xi,\lambda)=-D^2+\lambda^2 N\ }
$$
(« le carré du scaling plus la graduation »). Le **radical de la forme quadratique de Weil** est le **spectre négatif** de $\omega$ — l'**espace de Sonin** — qui est aussi la **range** de la carte $\mathcal E(f)(u)=u^{1/2}\sum_{n\ge1}f(nu)$. La partie positive code les zéros bas.

---

## 2. Identification précise du flot dBN

Sur la fonction spectrale $\Xi(z)$ ($z$ $=$ hauteur du zéro $=$ spectre de $\mathcal S$), le flot dBN est $e^{-t\partial_z^2}$ (chaleur rétrograde, $\partial_t H_t=-\partial_z^2 H_t$). Or $\mathcal S^2=z^2$ et $\partial_z^2$ engendrent, avec $z\partial_z+\tfrac12$, un $\mathfrak{sl}(2)$ (l'algèbre de l'oscillateur côté spectral). Comme $\varpi(e_-)\propto i\,\partial^2$ est anti-hermitien, le sous-groupe $e^{\tau e_-}$ à $\tau$ **réel** est **unitaire** (évolution libre). Donc :

$$
\boxed{\ \text{flot dBN}\ =\ e^{-t\partial_z^2}\ =\ e^{\tau e_-}\ \text{avec}\ \tau\ \text{IMAGINAIRE}\ }
$$

— le semi-groupe de la chaleur est la **continuation analytique en temps imaginaire** de l'évolution libre unitaire métaplectique. C'est sa nature non unitaire qui le distingue des flots du groupe.

---

## 3. L'obstruction est invariante par conjugaison

Le flot dBN et le carré du scaling $\mathcal S^2$ sont les deux générateurs paraboliques **opposés** ($e_\mp$) du même $\mathfrak{sl}(2)$. Tentation : faire « tourner » l'opérateur prolate par le flot dBN pour fermer le spectre négatif (Sonin). Mais :

> **La conjugaison par tout opérateur inversible préserve le spectre.** Donc aucune orbite de $SL(2,\mathbb R)$ par conjugaison — unitaire ou non, y compris le flot dBN — ne vide l'espace de Sonin.

Vérifié numériquement (`metaplectic_dbn_explorer.py`) : les flots unitaires préservent le spectre à $10^{-12}$ ; le flot dBN comme conjugaison aussi (et sa norme d'opérateur explose, $1\to10^{14}$ pour $t=0\to0{,}05$, confirmant la chaleur non unitaire). **L'obstruction de Sonin est $SL(2)$-invariante.** Le seul rôle concevable du flot dBN est donc comme déformation de l'**état** $\xi$ / de la **forme** de Weil — pas comme conjugaison de l'opérateur. C'est le verrou conceptuel, que la théorie des groupes seule ne franchit pas.

---

## 4. Le mur du calcul fidèle

Construit fidèlement (`prolate_weil_dbn.py`, $\omega=-D^2+\lambda^2 N$, carte $\mathcal E$), le modèle à **une seule place archimédienne** :
- ne reproduit **pas** les zéros de $\zeta$ (les valeurs propres positives ne donnent ni $14{,}13$ ni $21{,}02\dots$) ;
- a un spectre négatif (Sonin) **croissant avec la troncature** (il est de dimension infinie — l'UV).

Le vrai accord aux zéros et la vraie forme de Weil exigent le cadre **semilocal complet** (toutes les places, carte $\mathcal E$ adélique, domaines d'auto-adjonction précis) — l'appareil technique de CCM. Un test fidèle de la monotonie de la Note III en NCG demanderait de réimplémenter ce cadre : entreprise lourde, exposée aux erreurs, au **gain attendu d'un simple test de cohérence** (il confirmerait que la collision dBN coïncide avec le bord de positivité, sans prouver $\Lambda=0$).

---

## 5. Verdict

Le thread a livré une **localisation structurelle propre**, pas un levier :
- le flot dBN est le membre **en temps imaginaire** du $\mathfrak{sl}(2)$ métaplectique de CCM ;
- l'espace de Sonin (radical de Weil) est **invariant par conjugaison** — le flot dBN ne peut l'éliminer par une orbite de groupe ;
- son seul rôle possible (déformation d'état) est précisément le verrou, et le calcul fidèle qui le testerait bute sur le cadre semilocal de CCM.

Le flot dBN reste une **reformulation** de HR ($\Lambda=0$), non un outil de preuve. L'obstacle contraignant demeure celui de la Note V : le lemme analytique de positivité de trace. **Le thread dBN↔NCG est clos.**

---

## 6. Établi / conjectural

| énoncé | statut |
|---|---|
| opérateur prolate de CCM $=-D^2+\lambda^2 N$ ; Sonin $=$ spectre négatif $=$ radical de Weil | **théorème** (CCM) ; opérateur construit fidèlement |
| flot dBN $=$ $e_-$ métaplectique en temps imaginaire | **établi** (identification ; vérifié numériquement) |
| Sonin invariant par conjugaison $SL(2,\mathbb R)$ | **établi** (élémentaire, mais clarifiant : exclut la voie « conjugaison ») |
| rôle du flot dBN comme déformation d'état fermant la positivité | **ouvert** (le verrou ; non franchi) |
| test fidèle en NCG (Note III semilocale) | **hors portée tractable** (exige le cadre semilocal de CCM) |
| levier sur HR | **non** ; reformulation, pas preuve |

---

### Scripts

- `metaplectic_dbn_explorer.py` — générateurs métaplectiques (relations $\mathfrak{sl}(2)$ à $10^{-13}$, Casimir $-3/4$), opérateur prolate, flots de groupe, flot dBN ; démontre l'invariance par conjugaison et la nature non unitaire (temps imaginaire) du flot dBN.
- `prolate_weil_dbn.py` — opérateur prolate exact $\omega=-D^2+\lambda^2 N$ (forme paire cyclique), carte $\mathcal E(f)(u)=u^{1/2}\sum f(nu)$, flot dBN gaussien ; montre crûment la limite « une seule place » (zéros non reproduits, Sonin de dimension infinie).

---

*Mise en œuvre numérique et rédaction assistées par Claude (Anthropic), sous la direction et la vérification de l'auteur. Les jugements de positionnement relèvent d'une lecture de la littérature et engagent l'auteur.*
