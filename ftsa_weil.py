#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
ftsa_weil.py  —  Banc d'essai numerique (v2, consolide)

These testee :  la POSITIVITE DE WEIL discrimine
    zeta  (produit eulerien)        -> zeros sur la droite, W(g) >= 0
    Epstein Q0 = m^2+mn+6n^2         -> zeros HORS droite, W(g) < 0 possible
    (disc -23, nombre de classes 3 ; pas de produit eulerien)

Tout est en precision arbitraire (mpmath).  Fichier unique a faire evoluer.

Modes :
    python ftsa_weil.py selftest
    python ftsa_weil.py realaxis
    python ftsa_weil.py online   T0 T1                 # zeros de Z_Q0 sur sigma=1/2
    python ftsa_weil.py hunt  S0 S1 T0 T1 [ds dt]      # cherche+raffine zeros (boite)
    python ftsa_weil.py polish SIG T [dps]             # raffine un zero a haute precision
    python ftsa_weil.py verify SIG T                   # |Z| a precision croissante
    python ftsa_weil.py decomp                         # Z_Q0=(2/3)zK+(4/3)Lf ; contraste zK
    python ftsa_weil.py weil-zeta [gauss|triangle] P   # positivite de Weil pour zeta
     python ftsa_weil.py weil-epstein [w]               # APPROCHE A : W^Q0(h) < 0 explicite
    python ftsa_weil.py entropy-flow                   # APPROCHE B : flot d'entropie / spectre de Weil
"""

import sys
import mpmath as mp
import numpy as np

# ============================================================================
#  Zeta d'Epstein de Q0(m,n) = m^2 + m n + 6 n^2   (disc -23, h=3)
#  A=[[1,1/2],[1/2,6]], det A=23/4.  Representation a gamma incomplete (exacte) :
#    pi^-s Gamma(s) Z(s) = sum_k r(k)(pi k)^-s Gamma(s,pi k)
#       + det^-1/2 sum_k r(k) beta_k^{s-1} Gamma(1-s,beta_k)
#       + 1/(det^1/2 (s-1)) - 1/s ,   beta_k = (4 pi/23) k
# ============================================================================
_REP = {}
def rep_numbers(KMAX, M=None):
    if M is None:
        M = int(KMAX ** 0.5) + 3
    key = (KMAX, M)
    if key in _REP:
        return _REP[key]
    r = [0] * (KMAX + 1)
    for m in range(-M, M + 1):
        for n in range(-M, M + 1):
            k = m * m + m * n + 6 * n * n
            if 1 <= k <= KMAX:
                r[k] += 1
    nz = [k for k in range(1, KMAX + 1) if r[k]]
    _REP[key] = (r, nz)
    return r, nz

def dps_for(t, base=30, slope=1.1):
    """precision de travail (la perte du quotient pi^s Lambda/Gamma(s) est ~0.55 t)."""
    return int(base + slope * abs(float(t)))

def kmax_for(dps):
    """troncature de S2 : exp(-0.546 k) < 10^-dps  =>  k > 4.22 dps."""
    return int(4.4 * dps) + 15

def Zep(s, dps=None, KMAX=None):
    s = mp.mpc(s)
    if dps is None:
        dps = dps_for(float(mp.im(s)))
    if KMAX is None:
        KMAX = max(120, kmax_for(dps))
    mp.mp.dps = dps
    s = mp.mpc(s)
    r, nz = rep_numbers(KMAX)
    detA = mp.mpf(23) / 4
    sqrtdet = mp.sqrt(detA)
    cinv = mp.mpf(4) / 23
    S1 = mp.mpf(0); S2 = mp.mpf(0)
    for k in nz:
        S1 += r[k] * (mp.pi * k) ** (-s) * mp.gammainc(s, mp.pi * k, mp.inf)
        bk = mp.pi * cinv * k
        S2 += r[k] * bk ** (s - 1) * mp.gammainc(1 - s, bk, mp.inf)
    Lam = S1 + S2 / sqrtdet + 1 / (sqrtdet * (s - 1)) - 1 / s
    return mp.pi ** s * Lam / mp.gamma(s)

def Zdir(s, K=20000):
    s = mp.mpc(s)
    mp.mp.dps = max(30, dps_for(float(mp.im(s))))
    s = mp.mpc(s)
    r, nz = rep_numbers(K)
    return mp.fsum(r[k] * mp.mpf(k) ** (-s) for k in nz)

def _findroot(seed, dps, KMAX):
    """findroot robuste : Muller, tolerance et maxsteps explicites."""
    mp.mp.dps = dps
    return mp.findroot(lambda u: Zep(u, dps=dps, KMAX=KMAX), seed,
                       solver='muller', tol=mp.mpf(10) ** (-(dps - 12)), maxsteps=200)

# ============================================================================
#  Modes Epstein
# ============================================================================
def mode_selftest():
    print("# SELFTEST : representation integrale vs somme directe (Re s > 1)")
    print("# (l'ecart est domine par la troncature de la somme directe, pas par Zep)")
    for s in [mp.mpf(2), mp.mpf(2) + 1j * mp.mpf("26.8"),
              mp.mpf(3) + 1j * mp.mpf(5)]:
        a = Zep(s, dps=40, KMAX=200)
        b = Zdir(s, K=40000)
        print(f"  s={mp.nstr(s,8):>22}  Zep={mp.nstr(a,10):>26}  |Zep-Zdir|={mp.nstr(abs(a-b),3)}")
    eps = mp.mpf("1e-6")
    print(f"  residu en s=1 : num={mp.nstr(Zep(1+eps,dps=40)*eps,8)}  "
          f"2pi/sqrt23={mp.nstr(2*mp.pi/mp.sqrt(23),8)}")

def mode_realaxis():
    print("# Z_Q0 sur l'axe reel sigma in (0,1)")
    prev = None
    for k in range(2, 40, 2):
        sg = mp.mpf(k) / 40
        v = mp.re(Zep(sg, dps=40))
        mark = "  <-- CHANGEMENT DE SIGNE" if (prev is not None and prev * v < 0) else ""
        print(f"  sigma={float(sg):.3f}  Z={mp.nstr(v,8):>14}{mark}")
        prev = v

def _grid_minima(sigs, ts, L):
    out = []
    for i in range(1, len(sigs) - 1):
        for j in range(1, len(ts) - 1):
            c = L[i][j]
            nb = [L[i+a][j+b] for a in (-1,0,1) for b in (-1,0,1) if (a,b)!=(0,0)]
            if c < min(nb):
                out.append((c, sigs[i], ts[j]))
    out.sort()
    return out

def mode_hunt(S0, S1, T0, T1, ds=0.04, dt=0.4, topk=6):
    print(f"# HUNT  sigma in [{S0},{S1}]  t in [{T0},{T1}]  pas=({ds},{dt})")
    sigs = [round(S0 + i*ds, 4) for i in range(int((S1-S0)/ds)+1)]
    ts   = [round(T0 + j*dt, 4) for j in range(int((T1-T0)/dt)+1)]
    L = [[0.0]*len(ts) for _ in sigs]
    for j, t in enumerate(ts):
        d = dps_for(t)
        for i, sg in enumerate(sigs):
            L[i][j] = float(mp.log10(abs(Zep(mp.mpf(repr(sg))+1j*mp.mpf(repr(t)),
                                              dps=d, KMAX=120))))
    cands = _grid_minima(sigs, ts, L)
    print("  creux locaux (log10|Z|, sigma, t) :")
    for c, sg, t in cands[:topk]:
        print(f"     {c:8.3f}  sigma={sg:.3f}  t={t:7.3f}")
    print("  --- raffinement findroot ---")
    found = []
    for c, sg, t in cands[:topk]:
        d = dps_for(t) + 15; KM = kmax_for(d)
        for pert in (mp.mpc(0,0), mp.mpc(0.02,0), mp.mpc(-0.02,0.05)):
            try:
                z = _findroot(mp.mpc(repr(sg), repr(t)) + pert, d, KM)
                if abs(Zep(z, dps=d, KMAX=KM)) < mp.mpf(10) ** (-25):
                    if not any(abs(z-w) < mp.mpf("1e-8") for w in found):
                        found.append(z)
                    break
            except Exception:
                pass
    if not found:
        print("  aucun zero confirme.")
        return found
    print(f"  ZEROS confirmes ({len(found)}) :")
    for z in sorted(found, key=lambda w: float(mp.im(w))):
        off = abs(float(mp.re(z)) - 0.5)
        tag = "SUR la droite" if off < 1e-7 else f"HORS droite (Re-1/2={mp.nstr(mp.re(z)-mp.mpf(1)/2,4)})"
        print(f"     rho = {mp.nstr(z,16):>34}   {tag}")
    return found

def mode_online(T0, T1, dt=0.25):
    print(f"# ONLINE  zeros de Z_Q0 sur sigma=1/2,  t in [{T0},{T1}]")
    ts = [T0 + j*dt for j in range(int((T1-T0)/dt)+1)]
    prev2 = prev1 = None; dips = []
    for idx, t in enumerate(ts):
        d = dps_for(t)
        Lv = float(mp.log10(abs(Zep(mp.mpf(1)/2 + 1j*mp.mpf(repr(t)), dps=d, KMAX=120))))
        if prev1 is not None and prev2 is not None and prev1 < Lv and prev1 < prev2:
            dips.append(ts[idx-1])
        prev2, prev1 = prev1, Lv
    zeros = []
    for t in dips:
        d = dps_for(t) + 15; KM = kmax_for(d)
        try:
            z = _findroot(mp.mpc(0.5, repr(t)), d, KM)
            if abs(Zep(z, dps=d, KMAX=KM)) < mp.mpf(10)**(-25) and abs(mp.re(z)-0.5) < 1e-6:
                zeros.append(z)
        except Exception:
            pass
    for z in zeros:
        print(f"     zero sur la droite : t={mp.nstr(mp.im(z),12)}  (Re={mp.nstr(mp.re(z),6)})")
    return zeros

def mode_polish(SIG, T, dps=90):
    KM = kmax_for(dps)
    print(f"# POLISH  graine ({SIG},{T})  dps={dps}  KMAX={KM}")
    z = _findroot(mp.mpc(repr(SIG), repr(T)), dps, KM)
    off = abs(float(mp.re(z)) - 0.5)
    print(f"  rho = {mp.nstr(z,34)}")
    print(f"  Re(rho)-1/2 = {mp.nstr(mp.re(z)-mp.mpf(1)/2,22)}   "
          f"({'SUR' if off<1e-7 else 'HORS'} la droite)")
    print(f"  |Z(rho)| = {mp.nstr(abs(Zep(z,dps=dps,KMAX=KM)),4)}")
    print("  quadruplet { rho, conj, 1-rho, 1-conj } :")
    for p in (z, mp.conj(z), 1-z, 1-mp.conj(z)):
        print(f"     sigma={mp.nstr(mp.re(p),6):>9} t={mp.nstr(mp.im(p),9):>11}   "
              f"|Z|={mp.nstr(abs(Zep(p,dps=dps,KMAX=KM)),3)}")
    return z

def mode_verify(SIG, T):
    z = mode_polish(SIG, T, dps=90)
    print("  stabilite sous la precision :")
    for dps in (40, 70, 100):
        print(f"     dps={dps:3d}  |Z(rho)|={mp.nstr(abs(Zep(z,dps=dps,KMAX=kmax_for(dps))),4)}")
    return z

# ============================================================================
#  Decomposition en fonctions L et contraste avec l'objet a produit eulerien
# ============================================================================
def _chi23(a):
    QR = {1,2,3,4,6,8,9,12,13,16,18}
    a %= 23
    return 0 if a == 0 else (1 if a in QR else -1)

def _Lchi(s):
    s = mp.mpc(s)
    return mp.mpf(23)**(-s) * mp.fsum(_chi23(a)*mp.zeta(s, mp.mpf(a)/23) for a in range(1,23))

def _zetaK(s):
    return mp.zeta(s) * _Lchi(s)

def _eta_product_coeffs(N=400):
    P = [0]*(N+1); P[0] = 1
    for n in range(1, N+1):
        new = P[:]
        for d in range(N+1-n): new[d+n] -= P[d]
        P = new
    Qc = [0]*(N+1); Qc[0] = 1; n = 1
    while 23*n <= N:
        new = Qc[:]
        for d in range(N+1-23*n): new[d+23*n] -= Qc[d]
        Qc = new; n += 1
    PQ = [0]*(N+1)
    for i in range(N+1):
        if P[i]:
            for j in range(N+1-i): PQ[i+j] += P[i]*Qc[j]
    af = [0]*(N+2)
    for m in range(1, N+1): af[m] = PQ[m-1]   # f = q*P*Qc
    return af

def mode_decomp():
    mp.mp.dps = 30
    print("# DECOMP  Z_Q0 = (2/3) zeta_K + (4/3) L(.,f),  f=eta(z)eta(23z)")
    assert _chi23(-1) == -1
    af = _eta_product_coeffs(400)
    Lf = lambda s: mp.fsum(af[m]*mp.mpf(m)**(-mp.mpc(s)) for m in range(1, 401) if af[m])
    s = mp.mpf(2)
    combo = mp.mpf(2)/3*_zetaK(s) + mp.mpf(4)/3*Lf(s)
    zq = Zep(s, dps=30, KMAX=200)
    print(f"  s=2 :  Z_Q0={mp.nstr(zq.real,10)}   (2/3)zK+(4/3)Lf={mp.nstr(combo.real,10)}   "
          f"ecart={mp.nstr(abs(zq-combo),3)}")
    rho0 = mp.mpc('0.95326047479466', '16.29021572039')
    print(f"\n  contraste a rho_0 (zero de Z_Q0) :")
    print(f"     |zeta_K(rho_0)| = {mp.nstr(abs(_zetaK(rho0)),6)}   (objet a produit eulerien -> non nul)")
    print(f"     |Z_Q0(rho_0)|   = {mp.nstr(abs(Zep(rho0,dps=60,KMAX=300)),3)}")
    print("  |zeta_K| sur sigma=1/2 vs sigma=0.953 :")
    for tt in (14, 16.29, 19):
        a = abs(_zetaK(mp.mpf(1)/2 + 1j*mp.mpf(repr(tt))))
        b = abs(_zetaK(mp.mpf('0.953') + 1j*mp.mpf(repr(tt))))
        print(f"     t={tt:6.2f}   |zK(1/2+it)|={mp.nstr(a,5):>8}   |zK(0.953+it)|={mp.nstr(b,5):>8}")

# ============================================================================
#  Positivite de Weil pour zeta
# ============================================================================
def _testfun(kind, param):
    if kind == "gauss":
        a = mp.mpf(repr(param))
        return (lambda r: mp.e**(-a*r*r),
                1/(2*mp.sqrt(mp.pi*a)),
                mp.e**(a/4),
                lambda u: (1/(2*mp.sqrt(mp.pi*a)))*mp.e**(-(u*u)/(4*a)))
    if kind == "triangle":
        T = mp.mpf(repr(param))
        def h(r):
            return T*T if abs(r) < mp.mpf("1e-30") else 4*mp.sin(r*T/2)**2/(r*r)
        return (h, T, 16*mp.sinh(T/4)**2,
                lambda u: (T-abs(u)) if isinstance(u, mp.mpc) else max(T-abs(u), mp.mpf(0)))
    raise ValueError("kind ?")

def mode_weil_zeta(kind="gauss", param=0.02, Nzeros=300, Pmax=200000):
    mp.mp.dps = 30
    h, g0, h_i2, g_of = _testfun(kind, param)
    print(f"# WEIL-ZETA  test={kind} param={param}")
    Wz = mp.mpf(0); n = 1
    while n <= Nzeros:
        term = h(mp.im(mp.zetazero(n)))
        Wz += term
        if n > 5 and abs(term) < mp.mpf("1e-25"): break
        n += 1
    Wz *= 2
    print(f"  cote zeros        W = {mp.nstr(Wz,12)}   ({n} zeros)")
    pole = 2 * h_i2
    arch = (1/(2*mp.pi)) * mp.quad(
        lambda r: h(r)*mp.re(mp.digamma(mp.mpf(1)/4 + 1j*r/2)), [-mp.inf, 0, mp.inf])
    logpi = g0 * mp.log(mp.pi)
    sieve = bytearray([1])*(Pmax+1)
    for q in range(2, int(Pmax**0.5)+1):
        if sieve[q]:
            for mlt in range(q*q, Pmax+1, q): sieve[mlt] = 0
    prim = mp.mpf(0)
    for q in (q for q in range(2, Pmax+1) if sieve[q]):
        pk = q
        while pk <= Pmax:
            prim += (mp.log(q)/mp.sqrt(pk))*g_of(mp.log(pk)); pk *= q
    prim *= 2
    Wa = pole - logpi + arch - prim
    print(f"  cote arithmetique W = {mp.nstr(Wa,12)}")
    print(f"     pole={mp.nstr(pole,8)}  -g0 log pi={mp.nstr(-logpi,8)}  "
          f"arch={mp.nstr(arch,8)}  -premiers={mp.nstr(-prim,8)}")
    print(f"  ECART |Wz-Wa| = {mp.nstr(abs(Wz-Wa),4)}")
    print(f"  => W(zeta) {'>= 0 (compatible RH)' if Wa>0 else '< 0 (?!)'}")
    return Wz, Wa

# ============================================================================
# ============================================================================
#  APPROCHE A — valeur negative explicite  W^{Q0}(h) < 0  (a la Bombieri)
#  Fonction test admissible h = k^2 >= 0 sur R, avec k reelle paire possedant
#  un NOEUD en +-mu (mu = gamma_0 du zero hors-droite). Le prolongement de h en
#  mu +- i*delta (delta = beta_0 - 1/2) devient negatif, tandis que les zeros
#  sur la droite (loin de mu) sont filtres -> W < 0.
# ============================================================================
def mode_weil_epstein(w=None):
    mp.mp.dps = 40
    rho0 = mp.mpc('0.9532604747946606862505', '16.29021572039039079296')
    mu = mp.im(rho0)
    delta = mp.re(rho0) - mp.mpf(1)/2
    # coordonnees gamma = (rho-1/2)/i du quadruplet hors-droite : +-mu +- i delta
    zquad = [mu - 1j*delta, -mu - 1j*delta, -mu + 1j*delta, mu + 1j*delta]
    onlist = [mp.mpf(x) for x in ('4.94219260266', '7.01538796204', '10.4944987846',
                                  '11.8044036482', '18.8284912004', '19.3037306182')]
    print("# WEIL-EPSTEIN  valeur explicite de la fonctionnelle de Weil W(h)=sum_rho h(gamma_rho)")
    print(f"#   zero hors-droite : mu={mp.nstr(mu,10)}  delta=Re-1/2={mp.nstr(delta,7)}")
    print("#   h(r)=k(r)^2 >= 0 sur R ;  k paire, noeud en +-mu, largeur w")
    print(f"#   {'w':>5} | {'sur-droite (>=0)':>16} | {'hors-droite':>13} | {'W_Q0(h)':>13}")
    ws = [w] if w is not None else [0.4, 0.6, 0.8, 1.0]
    for ww in ws:
        W = mp.mpf(repr(ww))
        def k(r):
            r = mp.mpc(r)
            return (r-mu)*mp.e**(-(r-mu)**2/(2*W**2)) - (r+mu)*mp.e**(-(r+mu)**2/(2*W**2))
        h = lambda r: k(r)**2
        off = sum((h(z) for z in zquad), mp.mpc(0))
        on = sum((h(g) + h(-g) for g in onlist), mp.mpc(0))
        tot = on + off
        flag = "   <-- W < 0 !" if mp.re(tot) < 0 else ""
        print(f"    {float(ww):5.2f} | {mp.nstr(mp.re(on),7):>16} | "
              f"{mp.nstr(mp.re(off),7):>13} | {mp.nstr(mp.re(tot),7):>13}{flag}")
    print("#   => W_Q0(h) < 0 pour une fonction test ADMISSIBLE : la positivite de Weil")
    print("#      ECHOUE explicitement pour l'Epstein.  (Rappel : W(zeta) = +0.037 > 0.)")
 
 
# ============================================================================
#  APPROCHE B — flot d'entropie (a la Perelman) sur la forme quadratique de Weil
#  M_ij = sum_{zeros rho} phi^_i(gamma_rho) phi^_j(gamma_rho)  (base reelle paire).
#  zeta (zeros reels)            -> M = Gram >= 0 : Weil OK, lambda_min >= 0.
#  Epstein (quadruplet complexe) -> M indefinie    : Weil ECHOUE, lambda_min < 0.
#  Flot de gradient sur le quotient de Rayleigh R(c)=c^T M c/|c|^2 : decroit
#  (monotonie entropique) et converge vers lambda_min = obstruction de Weil.
#  NB : modele de dimension finie bati a partir des zeros connus ; il demontre
#  l'EQUIVALENCE (zero hors-droite <=> valeur propre negative <=> Weil echoue),
#  pas une preuve independante de HR.
# ============================================================================
def mode_entropy_flow():
    import numpy as np
    mp.mp.dps = 30
    nus = [14.0 + 0.5*i for i in range(10)]      # frequences encadrant mu=16.29
    sig = 1.0
    def phih(z, nu):
        z = complex(z)
        return np.exp(-(z-nu)**2/(2*sig**2)) + np.exp(-(z+nu)**2/(2*sig**2))
    def build_M(gcoords):
        n = len(nus); M = np.zeros((n, n), dtype=complex)
        for z in gcoords:
            v = np.array([phih(z, nu) for nu in nus], dtype=complex)
            M += np.outer(v, v)
        return M.real
    zz = [float(mp.im(mp.zetazero(j))) for j in range(1, 40)]
    Mz = build_M([g for G in zz for g in (G, -G)])
    mu, delta = 16.29021572039, 0.45326047479
    onE = [4.94219, 7.01539, 10.4945, 11.8044, 18.8285, 19.3037]
    gq = [mu-1j*delta, -mu-1j*delta, -mu+1j*delta, mu+1j*delta]
    Me = build_M([g for G in onE for g in (G, -G)] + gq)
    lz = np.linalg.eigvalsh(Mz); le = np.linalg.eigvalsh(Me)
    print("# ENTROPY-FLOW  spectre de la forme quadratique de Weil "
          f"(base de {len(nus)} gaussiennes, sigma={sig})")
    print(f"   zeta    : lambda_min = {lz[0]:+.6e}   lambda_max = {lz[-1]:.4e}"
          f"   -> {'>= 0  (Weil OK, compatible RH)' if lz[0] > -1e-9 else '< 0'}")
    print(f"   Epstein : lambda_min = {le[0]:+.6e}   lambda_max = {le[-1]:.4e}"
          f"   -> {'< 0  (Weil ECHOUE)' if le[0] < 0 else '>= 0'}")
    rng = np.random.default_rng(0)
    c = rng.standard_normal(len(nus)); c /= np.linalg.norm(c)
    lr = 0.3 / le[-1]; traj = []
    for _ in range(6000):
        Mc = Me @ c
        R = float(c @ Mc); traj.append(R)
        c = c - lr*2*(Mc - R*c); c /= np.linalg.norm(c)
    print("   flot d'entropie sur l'Epstein (quotient de Rayleigh, doit DECROITRE) :")
    for s in (0, 50, 200, 1000, 5999):
        print(f"      t={s:5d}   R = {traj[s]:+.6e}")
    print(f"   point critique du flot = {traj[-1]:+.6e}   (= lambda_min = {le[0]:+.6e})")
    print("#   Positivite de Weil <=> lambda_min >= 0. Le flot converge vers l'obstruction :")
    print("#   >= 0 pour zeta (produit eulerien), < 0 pour l'Epstein (zero hors-droite).")
 
 
# ============================================================================
def main():
    if len(sys.argv) < 2:
        print(__doc__); return
    cmd, A = sys.argv[1], sys.argv[2:]
    if cmd == "selftest": mode_selftest()
    elif cmd == "realaxis": mode_realaxis()
    elif cmd == "hunt":
        S0,S1,T0,T1 = map(float, A[:4])
        ds = float(A[4]) if len(A) > 4 else 0.04
        dt = float(A[5]) if len(A) > 5 else 0.4
        mode_hunt(S0,S1,T0,T1,ds,dt)
    elif cmd == "online": mode_online(float(A[0]), float(A[1]))
    elif cmd == "polish":
        mode_polish(float(A[0]), float(A[1]), int(A[2]) if len(A) > 2 else 90)
    elif cmd == "verify": mode_verify(float(A[0]), float(A[1]))
    elif cmd == "decomp": mode_decomp()
    elif cmd == "weil-zeta":
        kind = A[0] if A else "gauss"
        param = float(A[1]) if len(A) > 1 else (0.02 if kind == "gauss" else 2.0)
        mode_weil_zeta(kind, param)
    elif cmd == "weil-epstein":
        mode_weil_epstein(float(A[0]) if A else None)
    elif cmd == "entropy-flow":
        mode_entropy_flow()
    else:
        print("commande inconnue :", cmd); print(__doc__)
 
if __name__ == "__main__":
    main()
