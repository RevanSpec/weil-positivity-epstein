"""
dbn_weil_coupling.py -- ETAPE 2 du test-pont : couple le flot de de Bruijn-Newman
a la forme quadratique de Weil.

A chaque t, on reconstruit la forme de Weil COTE-ZEROS a partir des zeros de H_t^{Q0}
(quadruplet hors-droite mu(t) +- i*delta(t)  ET  zeros en-ligne, TOUS flottes sous le
flot de la chaleur), avec EXACTEMENT les memes fonctions test localisees que weil_local
(centres nus = 14.0..18.5, largeur a = 0.5) :

    M_ij(t) = Re sum_{z zero de H_t}  u_i(z) u_j(z),
    u_i(z)  = exp(-a (z-nu_i)^2) + exp(-a (z+nu_i)^2),

puis lambda_min(t) = plus petite valeur propre de M(t).

MECANISME : un zero reel donne v v^T reel (rang 1, >=0). Le quadruplet z = mu +- i delta
donne 4 Re[v v^T] = 4 (v_R v_R^T - v_I v_I^T) : la part -4 v_I v_I^T est la direction
NEGATIVE, de norme ~ 4 |v_I|^2 ~ C delta^2. Quand delta(t) -> 0 elle s'evanouit.

PREDICTION (compatibilite D1<->D2) : lambda_min(t) -> 0 quand t -> Lambda_Q0,
LINEAIREMENT (lambda_min ~ -C' (Lambda_Q0 - t)), avec passage par zero en
t = Lambda_Q0 = 0.0843066945...  Ce serait la preuve que le flot de de Bruijn-Newman (A)
EST une descente pour l'obstruction de Weil (B) -> meme dynamique, vue des deux cotes.

Cross-check : a t=0 on doit retrouver lambda_min ~ -0.7045 (valeur de weil_local).

Depend de dbn_epstein.py (Ht, gh_nodes_mp) + numpy + mpmath.

Usage :
  python dbn_weil_coupling.py [T0 T1 NT dps NGH KMAX pw]
  recommande : python dbn_weil_coupling.py 0 0.0840 14 28 20 55 2
     T1 pres de Lambda_Q0=0.08431 ; grille concentree pres de la collision (pw=2).
     ~15-25 min CPU. Pour resserrer encore : T1=0.0842, NT=16, pw=2.5.
"""
import sys
import cmath
import numpy as np
import mpmath as mp
from dbn_epstein import Ht, gh_nodes_mp, _kmax, MU0, DELTA0

A = 0.5
NUS = [14.0 + 0.5 * i for i in range(10)]          # 14.0, 14.5, ..., 18.5
# zeros EN-LIGNE (coord gamma = z) a t=0, recenses et valides ; on garde la fenetre
# utile [~9.5, 23] (hors fenetre : contribution < 1e-12 aux fonctions test localisees).
ONZ0 = [10.494, 11.804, 12.972, 14.63, 18.828, 19.304, 20.348, 21.885]

def u_c(z, nu):
    return cmath.exp(-A * (z - nu) ** 2) + cmath.exp(-A * (z + nu) ** 2)

def lam_min(gcoords):
    """plus petite valeur propre de M_ij = Re sum_z u_i(z) u_j(z) (numpy float64)."""
    n = len(NUS)
    M = np.zeros((n, n), dtype=complex)
    for z in gcoords:
        v = np.array([u_c(z, nu) for nu in NUS], dtype=complex)
        M += np.outer(v, v)
    return float(np.linalg.eigvalsh(M.real)[0])

def _seeds_real(x):
    e = mp.mpf('1e-4'); return (x, x + e, x - e)
def _seeds_cplx(z):
    e = mp.mpf('1e-4'); return (z, z + e, z + 1j * e)

def main(T0=0.0, T1=0.0840, NT=14, dps=28, NGH=20, KMAX=55, pw=2.0):
    mp.mp.dps = dps
    nodes = gh_nodes_mp(NGH, dps)
    tol = mp.mpf(10) ** (-(dps - 12))
    print(f"# ETAPE 2 : couplage flot de de Bruijn-Newman <-> forme de Weil")
    print(f"#   t in [{T0},{T1}], {NT} pas, dps={dps}, KMAX={KMAX}, GH={NGH}")
    print(f"#   nus = 14.0..18.5 (10), a = {A} ; {len(ONZ0)} zeros en-ligne + quadruplet, TOUS flottes")
    print(f"#   Lambda_Q0 (collision, ref) = 0.0843066945")
    print(f"#   grille concentree pres de la collision (puissance {pw}); Lambda_Q0 ref = 0.0843066945")
    print(f"#   {'t':>9}  {'delta(t)':>13}  {'lambda_min(t)':>18}  {'-lam/delta^2':>13}")
    T0m, T1m = mp.mpf(T0), mp.mpf(T1)
    ts = [T1m - (T1m - T0m) * (mp.mpf(NT - 1 - i) / (NT - 1)) ** pw for i in range(NT)]
    zq = mp.mpc(MU0, DELTA0)                       # quadruplet : on suit mu + i*delta
    zon = [mp.mpf(str(g)) for g in ONZ0]           # zeros en-ligne (reels)
    rows = []
    for t in ts:
        tt = mp.mpf(t)
        f = lambda z: Ht(z, tt, dps, KMAX, nodes)
        # --- quadruplet : 1 zero suivi -> 4 par symetrie (H_t pair et reel) ---
        zq = mp.findroot(f, _seeds_cplx(zq), solver='muller', tol=tol, maxsteps=100)
        mu_t, del_t = mp.re(zq), mp.im(zq)
        quad = [complex(mu_t, del_t), complex(mu_t, -del_t),
                complex(-mu_t, del_t), complex(-mu_t, -del_t)]
        # --- zeros en-ligne : suivis, restent reels -> +-G ---
        gon = []
        for i in range(len(zon)):
            zr = mp.re(mp.findroot(f, _seeds_real(zon[i]), solver='muller',
                                   tol=tol, maxsteps=100))
            zon[i] = zr
            gon += [complex(zr, 0.0), complex(-zr, 0.0)]
        lam = lam_min(gon + quad)
        d = float(del_t)
        rij = (-lam / (d * d)) if d > 1e-9 else float('nan')
        print(f"  {float(tt):9.5f}  {d:13.8f}  {lam:18.10e}  {rij:13.6f}")
        rows.append((float(tt), d, lam))
    # --- estimation du croisement lambda_min(t)=0 : deux methodes ---
    rows = np.array(rows)
    LAMQ0 = 0.08430669450968749
    print("\n  Croisement lambda_min(t) = 0 :  [lin-t] = ajustement lineaire en t ;")
    print("                                  [lin-d2] = lambda ~ lambda0 - K*delta^2 (recommande)")
    print(f"    {'fen.':>4}  {'[lin-t] t0':>12}  {'lambda0':>12}  {'[lin-d2] t_cross':>17}  {'ecart/Lam':>10}")
    for k in (3, 4, 5, 6):
        if len(rows) >= k:
            tt = rows[-k:, 0]; d2 = rows[-k:, 1] ** 2; ll = rows[-k:, 2]
            s, b = np.polyfit(tt, ll, 1); cross_t = -b / s
            beta, alpha = np.polyfit(d2, ll, 1)          # lam = beta*delta^2 + alpha
            d2star = -alpha / beta                        # delta^2 ou lam=0
            p, q = np.polyfit(tt, d2, 1)                  # delta^2 = p*t + q
            t_cross = (d2star - q) / p
            ec = 100 * (LAMQ0 - t_cross) / LAMQ0
            print(f"    {k:>4}  {cross_t:>12.6f}  {alpha:>+12.3e}  {t_cross:>17.6f}  {ec:>+8.2f}%")
    print(f"    (Lambda_Q0 = {LAMQ0:.12f} ; lambda0>0 => croisement legerement SOUS Lambda_Q0)")

if __name__ == "__main__":
    a = sys.argv[1:]
    if a:
        T0, T1, NT = float(a[0]), float(a[1]), int(a[2])
        dps = int(a[3]) if len(a) > 3 else 28
        NGH = int(a[4]) if len(a) > 4 else 20
        KMAX = int(a[5]) if len(a) > 5 else 55
        pw = float(a[6]) if len(a) > 6 else 2.0
        main(T0, T1, NT, dps, NGH, KMAX, pw)
    else:
        main()
