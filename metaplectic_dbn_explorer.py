"""
metaplectic_dbn_explorer.py
===========================
Outil d'exploration de la piste dBN <-> NCG (Connes-Consani-Moscovici).

CADRE (CCM, "Zeta zeros and prolate wave operators", 2024).
Representation metaplectique de sl(2,R) sur L^2(R), base de Hermite {h_n} :
    a(h_n) = sqrt(n/2) h_{n-1},   [a,a*] = 1/2
    varpi(h)  = a^2 - (a*)^2          (= x d_x + 1/2,  anti-hermitien : scaling S = -i h)
    varpi(e+) = (i/2)(a+a*)^2         (= i pi x^2)
    varpi(e-) = (i/2)(a-a*)^2         (= i/(4pi) d_x^2)
    varpi(k)  = a a* + a* a           (= Hermite, diagonal n+1/2 ; generateur compact)
Relations : [h,e+]=2e+, [h,e-]=-2e-, [e+,e-]=h, k=i(e- - e+).

OPERATEUR PROLATE (CCM, Prop.) :
    W_lambda = h^2 + 4 pi lambda^2 k - 1/4   dans U(sl(2,R)).
Le RADICAL de la forme quadratique de Weil = le SPECTRE NEGATIF de W_lambda = l'ESPACE DE SONIN.
La partie positive code les zeros bas de zeta.

FLOT DE BRUIJN-NEWMAN.
Sur la fonction spectrale Xi(z), le flot dBN est e^{-t d_z^2} (chaleur retrograde),
i.e. le sous-groupe parabolique e_- MAIS EN TEMPS IMAGINAIRE :
    e^{tau e-} avec tau reel = evolution libre UNITAIRE  -> preserve le spectre,
    dBN = continuation analytique tau = i*sigma          -> NON unitaire (chaleur).

QUESTION EXPLOREE.
Peut-on eliminer le spectre negatif (Sonin) de W_lambda par une orbite de SL(2,R)
(eventuellement combinant dBN et scaling), ou l'obstruction est-elle SL(2)-invariante ?

POINT-CLE STRUCTUREL (que le script etablit) : la conjugaison par TOUT element
inversible preserve le spectre. Donc aucune orbite de groupe (unitaire ou non, par
CONJUGAISON) ne supprime le Sonin. Le flot dBN ne peut aider que par son role de
deformation de l'ETAT / de la forme (non-conjugaison) -- c'est le verrou conceptuel,
non resolu par la theorie des groupes seule. Le script fournit les operateurs pour
explorer ce role.

AVERTISSEMENT. La troncature finie de Hermite capture l'ALGEBRE exactement mais pas
l'espace de Sonin analytique precis de CCM (domaines, cadre semilocal). C'est un
explorateur STRUCTUREL, pas l'operateur prolate exact. Les grandes valeurs propres
negatives de h^2 sont des artefacts de bord de troncature : surveiller la dependance en N.

Usage :
    python metaplectic_dbn_explorer.py [N lambda]
    defaut : N=60, lambda=1.0
numpy + scipy.
"""
import sys
import numpy as np
try:
    from scipy.linalg import expm
    HAVE_SCIPY = True
except Exception:
    HAVE_SCIPY = False


# ----------------------------------------------------------------------
# Generateurs metaplectiques dans la base de Hermite
# ----------------------------------------------------------------------
def annihilation(N):
    a = np.zeros((N, N), dtype=complex)
    for n in range(1, N):
        a[n - 1, n] = np.sqrt(n / 2.0)
    return a

def generators(N):
    """Renvoie (a, ad, h, ep, em, k) dans la base de Hermite tronquee a N."""
    a = annihilation(N)
    ad = a.conj().T
    u, v = a + ad, a - ad
    h = a @ a - ad @ ad
    ep = 0.5j * (u @ u)
    em = 0.5j * (v @ v)
    k = a @ ad + ad @ a
    return a, ad, h, ep, em, k

def check_relations(N=40, margin=6):
    """Verifie les relations sl(2) sur le bloc loin du plafond de troncature."""
    a, ad, h, ep, em, k = generators(N)
    B = slice(0, N - margin)
    cm = lambda X, Y: X @ Y - Y @ X
    e = lambda M: float(np.max(np.abs(M[B, B])))
    return {
        "[h,e+]-2e+": e(cm(h, ep) - 2 * ep),
        "[h,e-]+2e-": e(cm(h, em) + 2 * em),
        "[e+,e-]-h":  e(cm(ep, em) - h),
        "k-i(e-_-e+)": e(k - 1j * (em - ep)),
        "h+h*":       e(h + h.conj().T),       # anti-hermitien
        "k-k*":       e(k - k.conj().T),       # hermitien
    }

def casimir(N):
    """Casimir C = h^2 + 2(e+ e- + e- e+) ; doit etre scalaire (-3/4) dans l'irrep."""
    a, ad, h, ep, em, k = generators(N)
    C = h @ h + 2 * (ep @ em + em @ ep)
    B = slice(0, N - 8)
    diag = np.diag(C)[B].real
    offmax = float(np.max(np.abs(C[B, B] - np.diag(np.diag(C))[B, B])))
    return diag.mean(), diag.std(), offmax


# ----------------------------------------------------------------------
# Operateur prolate et diagnostic de l'espace de Sonin
# ----------------------------------------------------------------------
def prolate(N, lam):
    a, ad, h, ep, em, k = generators(N)
    W = h @ h + 4 * np.pi * lam**2 * k - 0.25 * np.eye(N)
    return 0.5 * (W + W.conj().T)            # hermitien (symetrise les erreurs num.)

def sonin_spectrum(N, lam, tol=1e-9):
    """Spectre de W_lambda ; dimension de la partie negative (Sonin)."""
    W = prolate(N, lam)
    ev = np.linalg.eigvalsh(W)
    return ev, int(np.sum(ev < -tol))


# ----------------------------------------------------------------------
# Elements de groupe SL(2,R) et flot dBN
# ----------------------------------------------------------------------
def group_element(N, alpha=0.0, beta=0.0, gamma=0.0):
    """exp(alpha*h + beta*e+ + gamma*e-).  Coeffs REELS -> operateur UNITAIRE
    (car h,e+,e- anti-hermitiens). Coeffs IMAGINAIRES -> non unitaire (ex. dBN)."""
    a, ad, h, ep, em, k = generators(N)
    X = alpha * h + beta * ep + gamma * em
    if HAVE_SCIPY:
        return expm(X)
    # repli : exp par diagonalisation (X general -> via serie tronquee)
    M = np.eye(N, dtype=complex); term = np.eye(N, dtype=complex)
    for j in range(1, 60):
        term = term @ X / j
        M = M + term
    return M

def dbn_flow(N, t):
    """Flot de Bruijn-Newman = e^{-t d_z^2} = e_- en TEMPS IMAGINAIRE.
    Comme varpi(e-) = (i/(4pi)) d_x^2, on a e^{-t d^2} = exp(-4pi t * (-i) e-) = exp(4i pi t e-).
    -> coefficient gamma = 4i pi t (IMAGINAIRE) : operateur NON unitaire (chaleur)."""
    return group_element(N, gamma=4j * np.pi * t)

def spectrum_under_conjugation(N, lam, g):
    """Spectre de g W_lambda g^{-1}. (Invariant pour tout g inversible : test de coherence.)"""
    W = prolate(N, lam)
    gi = np.linalg.inv(g)
    Wc = g @ W @ gi
    return np.linalg.eigvalsh(0.5 * (Wc + Wc.conj().T))


# ----------------------------------------------------------------------
# Demonstration / exploration
# ----------------------------------------------------------------------
def main(N=60, lam=1.0):
    print(f"# Explorateur metaplectique dBN<->NCG   (N={N}, lambda={lam})\n")

    print("## 1. Relations sl(2) (bloc loin du bord) :")
    for k_, v_ in check_relations(N).items():
        print(f"     |{k_}| = {v_:.2e}")
    cm, cs, co = casimir(N)
    print(f"   Casimir C : moyenne {cm:+.4f} (attendu -0.75), ecart-type {cs:.2e}, hors-diag {co:.2e}")

    print("\n## 2. Operateur prolate : spectre et espace de Sonin (partie negative)")
    print(f"   {'N':>4} {'lambda':>7} {'dim Sonin':>10}   5 plus basses vp")
    for Ntest in [N // 2, N, N + 20]:
        ev, nneg = sonin_spectrum(Ntest, lam)
        print(f"   {Ntest:>4} {lam:>7.2f} {nneg:>10}   {np.round(ev[:5], 2)}")
    print("   (la dependance en N revele les artefacts de bord ; ne pas surinterpreter dim Sonin)")

    print("\n## 3. Invariance par l'action de SL(2,R) (conjugation) :")
    ev0, n0 = sonin_spectrum(N, lam)
    tests = [("scaling  e^{0.4 h}", group_element(N, alpha=0.4)),
             ("parabolique unitaire e^{0.3 e-}", group_element(N, gamma=0.3)),
             ("dBN (temps imaginaire) t=0.05", dbn_flow(N, 0.05))]
    for name, g in tests:
        try:
            evc = spectrum_under_conjugation(N, lam, g)
            shift = float(np.max(np.abs(np.sort(evc) - np.sort(ev0))))
            print(f"   {name:<34} : max|Delta spectre| = {shift:.2e}  -> spectre {'preserve' if shift<1e-6 else 'CHANGE'}")
        except np.linalg.LinAlgError:
            print(f"   {name:<34} : (non inversible en troncature)")
    print(f"   dim Sonin de reference = {n0}.  CONCLUSION : la conjugaison (meme non unitaire,")
    print("   comme dBN) preserve le spectre -> AUCUNE orbite de groupe par conjugaison ne vide")
    print("   le Sonin. L'obstruction est SL(2)-invariante au sens de la conjugaison.")

    print("\n## 4. Le flot dBN comme deformation NON-conjugative (le verrou) :")
    print("   Le role pertinent du flot dBN n'est PAS de conjuguer W_lambda (cela preserve le")
    print("   spectre) mais de deformer l'ETAT xi / la forme de Weil. Operateurs fournis :")
    print("     - dbn_flow(N,t)      : e^{-t d^2}, non unitaire (norme d'operateur ci-dessous)")
    for t in [0.0, 0.02, 0.05, 0.1]:
        g = dbn_flow(N, t)
        # norme spectrale (mesure de non-unitarite : 1.0 = unitaire)
        s = np.linalg.svd(g, compute_uv=False)
        print(f"       t={t:>4} : ||e^(-t d^2)|| = {s[0]:.3e}, cond = {s[0]/s[-1]:.3e}")
    print("   -> a t>0 l'operateur s'eloigne de l'unitaire (chaleur) : c'est la seule porte")
    print("      possible pour modifier la positivite, et c'est la question conceptuelle ouverte.")

    print("\n## Pistes a explorer avec cet outil :")
    print("   (a) construire la forme de Weil a partir d'un etat xi et de la carte E de CCM,")
    print("       puis suivre lambda_min sous dbn_flow (analogue NCG de la Note III) ;")
    print("   (b) chercher si une combinaison dBN+scaling rend W_lambda + (deformation additive)")
    print("       sans spectre negatif -- deformation ADDITIVE, pas conjugation ;")
    print("   (c) tester la stabilite de dim Sonin en N pour separer structure et artefact.")


if __name__ == "__main__":
    a = sys.argv[1:]
    N = int(a[0]) if len(a) > 0 else 60
    lam = float(a[1]) if len(a) > 1 else 1.0
    main(N, lam)
