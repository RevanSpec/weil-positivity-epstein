"""
fork2_function_field.py -- exploration du Fork 2 dans l'analogue CORPS DE FONCTIONS.

Cadre. Pour une courbe C/F_q de genre g, zeta_C(T)=P(T)/((1-T)(1-qT)),
P(T)=prod_{i=1}^{2g}(1-alpha_i T), alpha_i = sqrt(q) e^{i theta_i}.
HR (Weil, THEOREME) : |alpha_i|=sqrt(q), i.e. tous les theta_i reels.
Sur le cercle critique on pose x=cos(phi) ; la fonction completee reelle est
    Xi(x) = prod_j (x - cos theta_j)     (polynome reel en x ; zeros aux angles).
Un zero HORS-CERCLE = quadruplet d'angles +-beta +- i delta  ->  en x, une paire
complexe conjuguee cos(beta +- i delta) (analogue corps-de-fonctions de l'Epstein
/ Davenport-Heilbronn).

Les deux pistes du Fork 2, testees ici :

  PISTE phi (le flot). Analogue de de Bruijn-Newman : H_t multiplie le coefficient
  de frequence n par e^{t n^2}. En corps de fonctions la frequence duale a l'angle
  phi est l'entier n, et cos(n phi)=T_n(cos phi)=T_n(x). Donc :
      Xi_t(x) = somme_n ( b_n e^{t n^2} ) T_n(x),   b_n = coeffs de Tchebychev de Xi.
  C'est un flot EXACT, fini, sur les coefficients. On teste s'il pousse le
  quadruplet hors-cercle vers le cercle (delta(t) -> 0, collision Lambda_ff).
  [Question geometrique ouverte, signalee : ce flot-ombre se releve-t-il en une
   deformation de la correspondance de Frobenius sur C x C ? Le script teste l'ombre.]

  PISTE entropie. La positivite de Weil en corps de fonctions = positivite de la
  MESURE SPECTRALE sum delta_{theta_i} sur le cercle (Bochner), i.e. de la matrice
  de Toeplitz des moments  s_k = sum_i (alpha_i/sqrt q)^k = 2 sum_j T_k(x_j) :
      T_{nm} = s_{|n-m|},   lambda_min(T) >= 0  <=>  angles reels  <=>  HR.
  Ces s_k sont les donnees de comptage de points (Lefschetz : s_k ~ q^k+1-#C(F_{q^k})).
  lambda_min(T) est l'avatar de la FORME D'INTERSECTION sur C x C (indice de Hodge),
  dont la definitude est le THEOREME de Weil. On teste sa monotonie le long du flot.

  COINCIDENCE A<->B attendue : collision (delta->0) et croisement (lambda_min=0)
  au MEME temps Lambda_ff -- la meme dynamique, dans le cadre ou la positivite est
  un theoreme.

Tout est en numpy float64 (exact a la precision machine ; objets finis).

Usage :
  python fork2_function_field.py [delta0 beta NT t_max N]
  defaut : 0.40 1.45 26 0.20 8
"""
import sys
import numpy as np
from numpy.polynomial import chebyshev as C

# angles "sur le cercle" (HR vrai) du fond, en (0, pi) :
THETA_ON = [0.45, 1.05, 1.85, 2.55]

def x_zeros(theta_on, beta, delta0):
    """zeros en x = cos(phi). delta0>0 : ajoute un quadruplet hors-cercle."""
    xr = [np.cos(t) for t in theta_on]
    if delta0 == 0.0:
        return np.array(xr + [np.cos(beta), np.cos(beta)], dtype=complex)
    X = np.cos(beta + 1j * delta0)            # complexe
    return np.array(xr + [X, np.conj(X)], dtype=complex)

def cheb_of_xi(xz):
    """coeffs de Tchebychev b_0..b_D de Xi(x)=prod (x - x_j)."""
    pc = np.poly(xz)                          # base puissances, degre decroissant
    pc = np.real(pc[::-1])                     # croissant, reel (racines conj-fermees)
    return C.poly2cheb(pc)

def flowed_roots(b, t):
    """racines en x de Xi_t(x) = sum_n (b_n e^{t n^2}) T_n(x)."""
    bt = b * np.exp(t * np.arange(len(b)) ** 2)
    return C.chebroots(bt)

def smoments(xz, N):
    """s_k = 2 sum_j T_k(x_j), k=0..N (reels si xz conj-ferme)."""
    s = np.empty(N + 1)
    for k in range(N + 1):
        ck = np.zeros(k + 1); ck[k] = 1.0
        s[k] = 2.0 * np.real(np.sum(C.chebval(xz, ck)))
    return s

def lam_min_toeplitz(xz, N):
    """plus petite valeur propre de la Toeplitz des moments (forme de Weil/Bochner GLOBALE)."""
    s = smoments(xz, N)
    T = np.array([[s[abs(n - m)] for m in range(N + 1)] for n in range(N + 1)])
    return float(np.linalg.eigvalsh(T)[0])

# fonctions test LOCALISEES sur le quadruplet (miroir exact de l'etape 2 / Epstein)
NU_OFFSETS = [-0.6, -0.4, -0.2, 0.0, 0.2, 0.4, 0.6]
A_LOC = 3.0

def lam_min_local(xz, beta, a=A_LOC):
    """forme de Weil LOCALISEE : test-fcts gaussiennes en phi centrees sur beta.
    Mêmes ingredients que weil_local/etape 2 : M_ij = Re sum_{angles} u_i u_j."""
    phis = np.arccos(xz)                       # angles (complexes si hors-cercle)
    ang = np.concatenate([phis, -phis])        # +- phi (les 2g valeurs propres)
    nus = [beta + o for o in NU_OFFSETS]
    n = len(nus)
    M = np.zeros((n, n), dtype=complex)
    for ph in ang:
        v = np.array([np.exp(-a * (ph - nu) ** 2) + np.exp(-a * (ph + nu) ** 2)
                      for nu in nus])
        M += np.outer(v, v)
    return float(np.linalg.eigvalsh(M.real)[0])

def off_circle_root(roots, prev=None):
    """racine en x hors de l'axe reel (plus grand |Im|), ou la plus proche de prev."""
    if prev is None:
        return roots[np.argmax(np.abs(roots.imag))]
    return roots[np.argmin(np.abs(roots - prev))]

def main(delta0=0.40, beta=1.45, NT=26, t_max=0.20, N=8):
    print(f"# FORK 2 -- analogue corps de fonctions")
    print(f"#   angles fond (HR vrai) : {THETA_ON}")
    print(f"#   quadruplet hors-cercle : beta={beta}, delta0={delta0}; Toeplitz N={N}, a_loc={A_LOC}")

    # --- ANCRAGE : config tous-reels (HR vrai) -> formes de Weil >= 0 ---
    xz_anchor = x_zeros(THETA_ON, beta, 0.0)
    lam_T0 = lam_min_toeplitz(xz_anchor, N)
    lam_L0 = lam_min_local(xz_anchor, beta)
    print(f"\n  ANCRAGE (HR vrai) : lambda_min Toeplitz = {lam_T0:+.6e} ; "
          f"localisee = {lam_L0:+.6e}   ({'OK >=0' if min(lam_T0,lam_L0) > -1e-9 else '?'})")

    # --- TEST : quadruplet hors-cercle, puis flot dBN-analogue ---
    b = cheb_of_xi(x_zeros(THETA_ON, beta, delta0))
    print(f"\n  Flot Xi_t(x) = sum_n b_n e^(t n^2) T_n(x) :")
    print(f"  {'t':>9}  {'|Im x_off|':>12}  {'delta_phi':>11}  {'lam_min[Toeplitz]':>18}  {'lam_min[localise]':>18}")
    ts = np.linspace(0.0, t_max, NT)
    prev = None; rows = []
    for t in ts:
        r = flowed_roots(b, t)
        Xoff = off_circle_root(r, prev); prev = Xoff
        imx = abs(Xoff.imag)
        dphi = abs(np.arccos(Xoff).imag)
        lamT = lam_min_toeplitz(r, N)
        lamL = lam_min_local(r, beta)
        print(f"  {t:9.5f}  {imx:12.8f}  {dphi:11.8f}  {lamT:18.10e}  {lamL:18.10e}")
        rows.append((t, imx, dphi, lamT, lamL))

    rows = np.array(rows)
    msk = rows[:, 1] > 1e-6
    tt = rows[msk, 0]; d2 = rows[msk, 1] ** 2
    print("\n  Collision (|Im x_off|->0) vs croisements (lambda_min->0) :")
    # collision
    for k in (3, 4, 5):
        if len(tt) >= k:
            p, q = np.polyfit(tt[-k:], d2[-k:], 1); Lam_ff = -q / p
            out = f"    {k} pts : Lambda_ff(collision) = {Lam_ff:.6f}"
            for col, name in ((3, "Toeplitz"), (4, "localise")):
                ll = rows[msk, col]
                bb, aa = np.polyfit(d2[-k:], ll[-k:], 1)
                d2star = -aa / bb; tcr = (d2star - q) / p
                out += f" | {name}: lam0={aa:+.2e}, croise a t={tcr:.5f}"
            print(out)
    print(f"\n  Lecture : la forme LOCALISEE (lam0->0, croise ~ Lambda_ff) reproduit A<->B ;")
    print(f"            la Toeplitz GLOBALE (lam0>0, croise avant) montre que la coincidence")
    print(f"            croisement=collision est propre a la forme localisee / complete.")

if __name__ == "__main__":
    a = sys.argv[1:]
    d0 = float(a[0]) if len(a) > 0 else 0.40
    bt = float(a[1]) if len(a) > 1 else 1.45
    NT = int(a[2]) if len(a) > 2 else 26
    tm = float(a[3]) if len(a) > 3 else 0.20
    N = int(a[4]) if len(a) > 4 else 8
    main(d0, bt, NT, tm, N)
