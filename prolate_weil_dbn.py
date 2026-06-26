"""
prolate_weil_dbn.py
===================
Piste (a) : la forme de Weil de CCM via la carte E, et le flot dBN -- analogue NCG
de la Note III. Outil FIDELE au cadre de Connes-Consani-Moscovici (cas archimedien).

CONSTRUCTION FIDELE (CCM, "Zeta zeros and prolate wave operators", section 2).
Paire cyclique (D, xi) : D auto-adjoint = multiplication par s sur L^2(R, dmu),
mesure archimedienne  dmu propto |Gamma(1/4 + i s/2)|^2 ds. Dans la base des polynomes
orthogonaux (xi_n), D est la matrice de JACOBI tridiagonale de coefficients
    a_n = (1/2) sqrt((2n+1)(2n+2)),
et N le nombre (N xi_n = n xi_n). L'OPERATEUR PROLATE est
    omega(D, xi, lambda) = -D^2 + lambda^2 N
(Definition 'prolop'). Restreint aux secteurs pair/impair (grading gamma=exp(i pi N)),
ce sont les matrices de Jacobi J^+ / J^- (Prop. 'jacobprol') :
    J^+_{n,n}   = -a_{2n}^2 - a_{2n-1}^2 + 2n   lambda^2,   J^+_{n,n+1} = -a_{2n}  a_{2n+1}
    J^-_{n,n}   = -a_{2n}^2 - a_{2n+1}^2 + (2n+1)lambda^2,   J^-_{n,n+1} = -a_{2n+1}a_{2n+2}
Le RADICAL de la forme quadratique de Weil = la PARTIE NEGATIVE du spectre de omega
= l'ESPACE DE SONIN = la RANGE de la carte
    E(f)(u) := u^{1/2} sum_{n>=1} f(n u).

FLOT dBN. Sur la variable log-echelle s = log u, le flot dBN est la multiplication
par la gaussienne e^{t s^2} (chaleur retrograde dans la variable spectrale).

CE QUE CET OUTIL FAIT / NE FAIT PAS (honnete).
- Il construit l'operateur prolate EXACT de CCM (forme paire cyclique) et la carte E.
- Il N'EST QU'A UNE PLACE (archimedienne) : il ne reproduit pas les zeros de zeta
  (cela exige le cadre SEMILOCAL complet de CCM, toutes les places + domaines precis),
  et le spectre negatif (Sonin) est de dimension infinie + sensible a la troncature.
- C'est donc un outil pour explorer la STRUCTURE (operateur prolate, carte E, action du
  flot dBN sur le radical), pas pour calculer la positivite de Weil de zeta elle-meme.

Usage : python prolate_weil_dbn.py [M lambda]   (defaut M=120, lambda=2.0)
numpy seul.
"""
import sys
import numpy as np


# ----------------------------------------------------------------------
# Paire cyclique : matrice de Jacobi D, nombre N, operateur prolate omega
# ----------------------------------------------------------------------
def a_coef(n):
    return 0.5 * np.sqrt((2 * n + 1) * (2 * n + 2)) if n >= 0 else 0.0

def jacobi_D(M):
    """D = multiplication par s (matrice de Jacobi tridiagonale, base des OP)."""
    D = np.zeros((M, M))
    for n in range(M - 1):
        D[n, n + 1] = D[n + 1, n] = a_coef(n)
    return D

def number_N(M):
    return np.diag(np.arange(M, dtype=float))

def prolate_full(M, lam):
    """omega = -D^2 + lambda^2 N (forme directe, base des OP)."""
    D = jacobi_D(M)
    return -D @ D + lam**2 * number_N(M)

def prolate_sectors(M, lam):
    """Matrices de Jacobi J^+ (pair) et J^- (impair), Prop. jacobprol."""
    Jp = np.zeros((M, M)); Jm = np.zeros((M, M))
    for n in range(M):
        Jp[n, n] = -a_coef(2*n)**2 - a_coef(2*n-1)**2 + 2*n*lam**2
        Jm[n, n] = -a_coef(2*n)**2 - a_coef(2*n+1)**2 + (2*n+1)*lam**2
        if n + 1 < M:
            Jp[n, n+1] = Jp[n+1, n] = -a_coef(2*n) * a_coef(2*n+1)
            Jm[n, n+1] = Jm[n+1, n] = -a_coef(2*n+1) * a_coef(2*n+2)
    return Jp, Jm


# ----------------------------------------------------------------------
# Carte E (range = Sonin = radical de la forme de Weil)
# ----------------------------------------------------------------------
def E_map(f_vals, s_grid):
    """E(f)(u) = u^{1/2} sum_{n>=1} f(n u),  u = e^s.  f donne sur la grille s.
    Renvoie E(f) sur la meme grille u=e^s (interpolation pour f(n u))."""
    u = np.exp(s_grid)
    out = np.zeros_like(u)
    fmax_s = s_grid[-1]
    n = 1
    while True:
        s_shift = s_grid + np.log(n)          # log(n u) = s + log n
        if s_shift[0] > fmax_s:               # f essentiellement nul au-dela de la grille
            break
        fn = np.interp(s_shift, s_grid, f_vals, left=0.0, right=0.0)
        out += fn
        n += 1
        if n > 10000:
            break
    return np.sqrt(u) * out


# ----------------------------------------------------------------------
# Flot dBN : multiplication par la gaussienne e^{t s^2} dans la variable s
# ----------------------------------------------------------------------
def dbn_weight(s_grid, t):
    return np.exp(t * s_grid**2)


# ----------------------------------------------------------------------
# Exploration
# ----------------------------------------------------------------------
def main(M=120, lam=2.0):
    print(f"# Outil fidele piste (a) : prolate de CCM + carte E + flot dBN   (M={M}, lambda={lam})\n")

    print("## 1. Operateur prolate omega = -D^2 + lambda^2 N (forme paire cyclique CCM)")
    Jp, Jm = prolate_sectors(M, lam)
    evp = np.linalg.eigvalsh(Jp); evm = np.linalg.eigvalsh(Jm)
    ev = np.sort(np.concatenate([evp, evm]))
    nneg = int(np.sum(ev < -1e-9))
    print(f"   dim(spectre<0) = {nneg}  [candidat Sonin ; infini en realite + effet troncature]")
    print(f"   6 plus basses vp > 0 : {np.round(ev[ev>0][:6],3)}")
    # dependance en M : signaler l'instabilite (Sonin infini-dim)
    for Mt in [M//2, M, M+40]:
        Jp2, Jm2 = prolate_sectors(Mt, lam)
        n2 = int(np.sum(np.concatenate([np.linalg.eigvalsh(Jp2), np.linalg.eigvalsh(Jm2)]) < -1e-9))
        print(f"     M={Mt:>4} -> dim(spectre<0)={n2}")
    print("   => dim Sonin croit avec M (UV de dimension infinie) : ne pas surinterpreter.")

    print("\n## 2. Carte E : E(f)(u) = u^{1/2} sum_{n>=1} f(n u)  (range = Sonin = radical de Weil)")
    s = np.linspace(-6, 6, 1200)
    # un etat test even (gaussienne centree) et son image par E
    f = np.exp(-0.5 * s**2)
    Ef = E_map(f, s)
    print(f"   ||f|| = {np.linalg.norm(f):.3f}, ||E(f)|| = {np.linalg.norm(Ef):.3f}")
    print(f"   E(f) concentre la 'somme theta' ; sa range engendre l'espace de Sonin.")

    print("\n## 3. Flot dBN e^{t s^2} et son action sur E (le verrou : deforme-t-il le radical ?)")
    print(f"   {'t':>6} {'||e^{ts^2} f||':>14} {'<E(f),E(e^{ts^2}f)>/||.||':>26}")
    for t in [0.0, 0.02, 0.05, 0.1]:
        w = dbn_weight(s, t)
        ft = w * f
        Eft = E_map(ft, s)
        cos = float(np.dot(Ef, Eft) / (np.linalg.norm(Ef) * np.linalg.norm(Eft) + 1e-30))
        print(f"   {t:>6.2f} {np.linalg.norm(ft):>14.3e} {cos:>26.4f}")
    print("   (cos -> 1 : dBN preserverait la direction du radical ; cos qui derive : il la deforme)")

    print("\n## Lecture honnete :")
    print("   - L'operateur prolate construit est EXACT (forme paire cyclique de CCM).")
    print("   - Mais a UNE place : ne reproduit pas les zeros de zeta ni la vraie forme de Weil")
    print("     (cela exige le cadre SEMILOCAL complet : toutes les places, domaines, carte E adelique).")
    print("   - Pour un vrai test de la Note III en NCG, il faudrait reimplementer ce cadre semilocal")
    print("     -- entreprise technique lourde (l'appareil de CCM), au gain attendu = test de coherence.")


if __name__ == "__main__":
    a = sys.argv[1:]
    M = int(a[0]) if len(a) > 0 else 120
    lam = float(a[1]) if len(a) > 1 else 2.0
    main(M, lam)
