"""
dbn_epstein.py -- flot de de Bruijn-Newman pour la zeta d'Epstein Z_{Q0}.

But : estimer la constante de de Bruijn-Newman ANALOGUE  Lambda_{Q0} > 0,
le temps t auquel le zero hors-droite de Z_{Q0} (Q0 = m^2+mn+6n^2, disc -23)
rejoint la droite critique sous le flot de la chaleur H_t.

Reformulation NON-OSCILLANTE (semi-groupe de la chaleur dans s). Comme
z -> s = 1/2 + i z envoie d_z^2 sur -d_s^2, on a H_t = e^{t d_s^2} xi, soit
une convolution gaussienne ; apres normalisation de la gaussienne :

    H_t(z) = (1/sqrt(2pi)) * Int_R  e^{-v^2/2} xi(1/2 + sqrt(2t) v + i z) dv
           = E[ xi(1/2 + sqrt(2t) V + i z) ],   V ~ N(0,1).

xi_{Q0}(s) = (23/4)^{s/2} Lambda(s),  Lambda(s) = pi^{-s} Gamma(s) Z_{Q0}(s),
Lambda par la representation a gamma incomplete (exacte, validee dans ftsa_weil).

PRECISION : les noeuds de Gauss-Hermite sont calcules EN PRECISION ARBITRAIRE
(raffinement Newton de He_n a partir des graines numpy), sinon la precision
plafonne a ~15 chiffres (noeuds float64). C'est ce qui permet d'augmenter dps.

Depend de ftsa_weil.py (rep_numbers) + mpmath + numpy.

--------------------------------------------------------------------------
RUNS RECOMMANDES (haute precision) :

  # 1) validation rapide
  python dbn_epstein.py validate

  # 2) flot dense, dps 50, 48 noeuds : trace delta(t) + extrapolation racine
  python dbn_epstein.py flow 0 0.083 16 50 48

  # 3) valeur PRECISE de Lambda_Q0 par le zero double (recommande) :
  python dbn_epstein.py collision 0.0835 16.06 50 48
  #   -> resout {H_t(z)=0, d_z H_t(z)=0} en (t, z*).  Pour plus de chiffres :
  python dbn_epstein.py collision 0.0835 16.06 70 56
--------------------------------------------------------------------------
"""
import sys
import mpmath as mp
import numpy as np
from ftsa_weil import rep_numbers

# --- zero hors-droite connu rho0 ; coordonnee z = (rho-1/2)/i = mu + i*delta ---
RHO0_RE = mp.mpf('0.9532604747946606862505')
RHO0_IM = mp.mpf('16.290215720390390792963')
MU0    = RHO0_IM                  # ~ 16.2902
DELTA0 = RHO0_RE - mp.mpf('0.5')  # ~ 0.45326
Z0     = mp.mpc(MU0, DELTA0)

def _kmax(dps):
    return int(4.4 * dps) + 15

# ----- Hermite "probabiliste" He_n et noeuds de Gauss-Hermite en precision arbitraire
def _he_pair(n, x):
    """(He_n(x), He_{n-1}(x)) par recurrence He_{k+1} = x He_k - k He_{k-1}."""
    if n == 0:
        return mp.mpf(1), mp.mpf(0)
    hm1, h = mp.mpf(1), mp.mpf(x)     # He_0, He_1
    for k in range(1, n):
        hm1, h = h, x * h - k * hm1
    return h, hm1

def gh_nodes_mp(n, dps):
    """Noeuds/poids de Gauss-Hermite pour Int e^{-x^2/2} f dx, en precision dps."""
    work = dps + 10
    mp.mp.dps = work
    x0, _ = np.polynomial.hermite_e.hermegauss(n)        # graines float64
    fact = mp.factorial(n)
    sq2pi = mp.sqrt(2 * mp.pi)
    xs, ws = [], []
    for xg in x0:
        x = mp.mpf(float(xg))
        for _ in range(80):                              # Newton sur He_n
            Hn, Hnm1 = _he_pair(n, x)
            dx = Hn / (n * Hnm1)                          # He_n' = n He_{n-1}
            x -= dx
            if abs(dx) < mp.mpf(10) ** (-(work - 3)):
                break
        Hnm1 = _he_pair(n, x)[1]
        xs.append(x)
        ws.append(fact * sq2pi / (n * n * Hnm1 ** 2))
    # auto-controle : moments 0 et 2 (doivent valoir sqrt(2pi))
    m0 = mp.fsum(ws)
    m2 = mp.fsum(w * xx ** 2 for w, xx in zip(ws, xs))
    assert abs(m0 - sq2pi) < mp.mpf(10) ** (-(dps)), f"GH moment0 KO: {mp.nstr(m0-sq2pi,3)}"
    assert abs(m2 - sq2pi) < mp.mpf(10) ** (-(dps)), f"GH moment2 KO: {mp.nstr(m2-sq2pi,3)}"
    mp.mp.dps = dps
    return xs, ws

# ----- fonction completee d'Epstein
def Lam(s, KMAX):
    s = mp.mpc(s)
    r, nz = rep_numbers(KMAX)
    sqrtdet = mp.sqrt(mp.mpf(23) / 4)
    cinv = mp.mpf(4) / 23
    S1 = mp.mpf(0); S2 = mp.mpf(0)
    for k in nz:
        S1 += r[k] * (mp.pi * k) ** (-s) * mp.gammainc(s, mp.pi * k, mp.inf)
        bk = mp.pi * cinv * k
        S2 += r[k] * bk ** (s - 1) * mp.gammainc(1 - s, bk, mp.inf)
    return S1 + S2 / sqrtdet + 1 / (sqrtdet * (s - 1)) - 1 / s

def xi(s, KMAX):
    s = mp.mpc(s)
    return (mp.mpf(23) / 4) ** (s / 2) * Lam(s, KMAX)

def Ht(z, t, dps, KMAX, nodes):
    """H_t(z) = E[ xi(1/2 + sqrt(2t) V + i z) ], V~N(0,1), via Gauss-Hermite (mp)."""
    mp.mp.dps = dps
    z = mp.mpc(z)
<<<<<<< HEAD
    t = mp.mpf(t)                      # NE PAS passer par float : casserait les petits pas en t
    if t <= 0:
        return xi(mp.mpf('0.5') + 1j * z, KMAX)
    sig = mp.sqrt(2 * t)
=======
    if float(t) <= 0:
        return xi(mp.mpf('0.5') + 1j * z, KMAX)
    sig = mp.sqrt(2 * mp.mpf(t))
>>>>>>> bb725e4 (add dbn-weil)
    xs, ws = nodes
    half = mp.mpf('0.5')
    acc = mp.mpf(0)
    for xx, ww in zip(xs, ws):
        acc += ww * xi(half + sig * xx + 1j * z, KMAX)
    return acc / mp.sqrt(2 * mp.pi)

# ============================================================================
def mode_validate(dps=40, NGH=32):
    KMAX = _kmax(dps); mp.mp.dps = dps
    print(f"# VALIDATE  (dps={dps}, KMAX={KMAX}, GH={NGH})")
    nodes = gh_nodes_mp(NGH, dps)
    print(f"  noeuds GH en precision arbitraire : OK (moments 0,2 verifies)")
    v0 = xi(mp.mpf('0.5') + 1j * Z0, KMAX)
    print(f"  |H_0(z0)| = |xi(1/2+i z0)| = {mp.nstr(abs(v0), 6)}   (~0 attendu)")
    for (zr, t) in [(mp.mpf('5.0'), 0.1), (mp.mpf('16.29'), 0.2)]:
        h = Ht(mp.mpc(zr, 0), t, dps, KMAX, nodes)
        rel = abs(mp.im(h)) / (abs(h) + mp.mpf(10) ** (-dps))
        print(f"  H_{t}({mp.nstr(zr,6)}) reelle ? Im/|.| = {mp.nstr(rel,4)}")

def track_zero(t, z_seed, dps, KMAX, nodes):
    mp.mp.dps = dps
    f = lambda z: Ht(z, t, dps, KMAX, nodes)
    eps = mp.mpf('1e-4')
    seeds = (z_seed, z_seed + eps, z_seed + 1j * eps)
    return mp.findroot(f, seeds, solver='muller',
                       tol=mp.mpf(10) ** (-(dps - 14)), maxsteps=120)

def mode_flow(T0, T1, NT, dps=50, NGH=48):
    KMAX = _kmax(dps); mp.mp.dps = dps
    nodes = gh_nodes_mp(NGH, dps)
    print(f"# FLOW de Bruijn-Newman, zero hors-droite de Z_Q0")
    print(f"#   t in [{T0},{T1}], {NT} pas, dps={dps}, KMAX={KMAX}, GH={NGH}")
    print(f"#   z0 = {mp.nstr(MU0,14)} + i {mp.nstr(DELTA0,12)}")
    print(f"#   {'t':>9}  {'Re z(t)':>20}  {'delta(t)=Im z(t)':>24}  {'|H|':>10}")
    ts = [mp.mpf(T0) + (mp.mpf(T1) - mp.mpf(T0)) * i / (NT - 1) for i in range(NT)]
    z = mp.mpc(Z0); rows = []
    for t in ts:
        try:
            z = track_zero(float(t), z, dps, KMAX, nodes)
            resid = abs(Ht(z, float(t), dps, KMAX, nodes))
            d = mp.im(z)
            print(f"  {float(t):9.5f}  {mp.nstr(mp.re(z),16):>20}  {mp.nstr(d,16):>24}  {mp.nstr(resid,3):>10}")
            rows.append((float(t), float(mp.re(z)), float(d)))
            if float(d) < 2e-3:
                print("  -> delta < 2e-3 : collision proche, suivi direct arrete.")
                break
        except Exception as e:
            print(f"  {float(t):9.5f}  (echec: {e})"); break
    # extrapolation racine delta^2 ~ a (Lambda - t), sur plusieurs fenetres
    if len(rows) >= 3:
        print("\n  Extrapolation  delta(t)^2 ~ a (Lambda_Q0 - t) :")
        for k in (3, 4, 5, 6):
            if len(rows) >= k:
                sub = rows[-k:]
                tt = np.array([r[0] for r in sub]); dd = np.array([r[2] for r in sub]) ** 2
                A = np.vstack([tt, np.ones_like(tt)]).T
                sl, ic = np.linalg.lstsq(A, dd, rcond=None)[0]
                print(f"    {k} derniers points : a={-sl:.5f}, Lambda_Q0 = {(-ic/sl):.6f}")

def mode_collision(t0, z0, dps=50, NGH=48, maxit=12):
    """Zero double {H_t(z)=0, d_z H_t(z)=0} en (t, z*), z* reel.
    Newton 2D fait main : grille 3x3 en (t,z) -> H, d_zH et Jacobien d'un coup."""
    KMAX = _kmax(dps); mp.mp.dps = dps
    nodes = gh_nodes_mp(NGH, dps)
    h = mp.mpf(10) ** (-(dps // 3))     # pas en z (equilibre troncature/arrondi)
    k = mp.mpf(10) ** (-(dps // 3))     # pas en t
    tol = mp.mpf(10) ** (-(dps // 2 - 2))
    print(f"# COLLISION (zero double)  dps={dps}, KMAX={KMAX}, GH={NGH}")
    print(f"#   graine (t, z*) = ({t0}, {z0})")
    def H(tt, zz):
<<<<<<< HEAD
        return mp.re(Ht(mp.mpc(zz, 0), tt, dps, KMAX, nodes))
=======
        return mp.re(Ht(mp.mpc(zz, 0), float(tt), dps, KMAX, nodes))
>>>>>>> bb725e4 (add dbn-weil)
    t = mp.mpf(str(t0)); z = mp.mpf(str(z0))
    for it in range(maxit):
        H00 = H(t, z)
        Hzp = H(t, z + h); Hzm = H(t, z - h)
        Htp = H(t + k, z); Htm = H(t - k, z)
        Hpp = H(t + k, z + h); Hpm = H(t + k, z - h)
        Hmp = H(t - k, z + h); Hmm = H(t - k, z - h)
        F1 = H00
        F2 = (Hzp - Hzm) / (2 * h)                       # d_z H
        J11 = (Htp - Htm) / (2 * k)                      # d_t H
        J12 = F2                                         # d_z H
        J21 = (Hpp - Hpm - Hmp + Hmm) / (4 * h * k)      # d_t d_z H
        J22 = (Hzp - 2 * H00 + Hzm) / h ** 2             # d_zz H
        det = J11 * J22 - J12 * J21
        dt = -(J22 * F1 - J12 * F2) / det
        dz = -(-J21 * F1 + J11 * F2) / det
        t += dt; z += dz
        print(f"  it {it:2d}: t={mp.nstr(t,14)}  z*={mp.nstr(z,14)}  "
              f"|H|={mp.nstr(abs(F1),2)} |d_zH|={mp.nstr(abs(F2),2)}")
        if abs(dt) < tol and abs(dz) < tol:
            break
    print(f"\n  Lambda_Q0 = t*       = {mp.nstr(t, dps - 14)}")
    print(f"  z* (lieu collision)  = {mp.nstr(z, dps - 14)}")

if __name__ == "__main__":
    a = sys.argv[1:]
    if not a or a[0] == "validate":
        dps = int(a[1]) if len(a) > 1 else 40
        NGH = int(a[2]) if len(a) > 2 else 32
        mode_validate(dps, NGH)
    elif a[0] == "flow":
        T0, T1, NT = float(a[1]), float(a[2]), int(a[3])
        dps = int(a[4]) if len(a) > 4 else 50
        NGH = int(a[5]) if len(a) > 5 else 48
        mode_flow(T0, T1, NT, dps, NGH)
    elif a[0] == "collision":
        t0, z0 = float(a[1]), float(a[2])
        dps = int(a[3]) if len(a) > 3 else 50
        NGH = int(a[4]) if len(a) > 4 else 48
        mode_collision(t0, z0, dps, NGH)
    else:
<<<<<<< HEAD
        print(__doc__)
=======
        print(__doc__)
>>>>>>> bb725e4 (add dbn-weil)
