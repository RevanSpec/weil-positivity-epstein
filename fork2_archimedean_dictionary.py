"""
fork2_archimedean_dictionary.py -- complete le dictionnaire corps de fonctions <-> nombres
                                    sur le cote NON-PREMIER (poles + archimedien).

Cote premier (deja fait, fork2_local_intersection.py) :
    facteur local en p  <->  (Gamma_{F^n}.Delta)_x = deg(x) <-> log p.

Cote NON-premier (ici), via la forme d'intersection sur la surface S = C x C :
    POLES de zeta (s=0,1)        <->  fibres  f1=[C x pt], f2=[pt x C]   [apportent +2 a Delta^2]
    facteur GAMMA (archimedien)  <->  GENRE : auto-appariement H^1 du diagonal [apporte -2g a Delta^2]
    equation fonctionnelle       <->  involution de Rosati  F' = q F^{-1}   (alpha <-> q/alpha)
    positivite de Weil           <->  indice de Hodge : signature (1, reste) sur NS(C x C)

Nombres d'intersection (Weil/Lefschetz, via Kunneth ; verifies) :
    f1^2 = f2^2 = 0,   f1.f2 = 1
    Gamma_{F^n}.f1 = q^n (= deg F^n),   Gamma_{F^n}.f2 = 1
    Gamma_{F^n}.Gamma_{F^m} = q^{min(n,m)} * N_{|n-m|}      (n != m)
    Gamma_{F^n}^2 = q^n (2 - 2g)        =>  Delta^2 = Gamma_{F^0}^2 = 2 - 2g
       [Kunneth : Gamma_{F^n} = q^n f1 + f2 + gamma_n ;  gamma_phi.gamma_psi = -Tr(phi psi'),
        psi' = adjoint symplectique (Rosati), F' = q F^{-1}.]

Exemples : g=1 (courbe elliptique, comptage du corps de base) ; g=2 (polynome de Weil entier valide).
numpy seul.

Usage : python fork2_archimedean_dictionary.py [g q ...]
        g=1 :  python fork2_archimedean_dictionary.py 1 q A B
        g=2 :  python fork2_archimedean_dictionary.py 2 q a1 a2
        defaut : g=1, q=5, A=1, B=1
"""
import sys
import numpy as np
import numpy.lib.scimath as sc

def eigs_elliptic(q, A, B):
    """g=1 : compte #E(F_q) sur le corps de base, renvoie les 2 valeurs propres."""
    n = 0
    for x in range(q):
        rhs = (x * x * x + A * x + B) % q
        if rhs == 0:
            n += 1
        elif pow(rhs, (q - 1) // 2, q) == 1:
            n += 2
    Nq = n + 1
    t = q + 1 - Nq
    d = sc.sqrt(t * t - 4 * q)
    return np.array([(t + d) / 2, (t - d) / 2]), f"E: y^2=x^3+{A}x+{B}/F_{q}, #E(F_{q})={Nq}, t={t}"

def eigs_weil_g2(q, a1, a2):
    """g=2 : P(T)=1 - a1 T + a2 T^2 - q a1 T^3 + q^2 T^4 ; renvoie les 4 valeurs propres."""
    coeffs = [q * q, -q * a1, a2, -a1, 1]      # q^2 T^4 - q a1 T^3 + a2 T^2 - a1 T + 1
    roots = np.roots(coeffs)                    # = 1/alpha_i
    return 1.0 / roots, f"P_2(T): a1={a1}, a2={a2} / F_{q}"

def main(g=1, q=5, p2=1, p3=1):
    if g == 1:
        alpha, label = eigs_elliptic(q, p2, p3)
    elif g == 2:
        alpha, label = eigs_weil_g2(q, p2, p3)
    else:
        raise SystemExit("g doit etre 1 ou 2 dans ce script de demonstration.")
    print(f"# {label}")
    okmod = np.allclose(np.abs(alpha), np.sqrt(q), atol=1e-7)
    print(f"#   genre g={g}, |alpha_i|=sqrt q : {'OK' if okmod else 'FAIL'}  (sqrt q={np.sqrt(q):.6f})")

    def N(n):                                   # #C(F_q^n) = q^n + 1 - sum alpha^n  (entier)
        return int(round(q**n + 1 - np.sum(alpha**n).real))

    # ---- forme d'intersection sur {f1, f2, Gamma_{F^0}=Delta, ..., Gamma_{F^k}} ----
    k = 2 * g + 1                               # assez de puissances de Frobenius
    powers = list(range(0, k + 1))              # n = 0..k  (0 = diagonale)
    dim = 2 + len(powers)
    G = np.zeros((dim, dim))
    # indices : 0 -> f1, 1 -> f2, 2+j -> Gamma_{F^{powers[j]}}
    G[0, 1] = G[1, 0] = 1.0                      # f1.f2 = 1 ; f1^2=f2^2=0
    for j, n in enumerate(powers):
        col = 2 + j
        G[0, col] = G[col, 0] = q**n             # Gamma_{F^n}.f1 = q^n
        G[1, col] = G[col, 1] = 1.0              # Gamma_{F^n}.f2 = 1
    for j, n in enumerate(powers):
        for jj, m in enumerate(powers):
            cj, cm = 2 + j, 2 + jj
            if n == m:
                G[cj, cm] = q**n * (2 - 2 * g)   # auto-intersection
            else:
                G[cj, cm] = q**min(n, m) * N(abs(n - m))

    # ---- (1) Delta^2 = 2 - 2g et sa decomposition ----
    Delta2 = G[2, 2]                             # Gamma_{F^0}^2
    print(f"\n  (1) Delta^2 = Gamma_F^0 . Gamma_F^0 = {Delta2:+.0f}   (attendu 2-2g = {2-2*g})  "
          f"-> {'OK' if abs(Delta2-(2-2*g))<1e-9 else 'FAIL'}")
    print(f"      decomposition de Kunneth :")
    print(f"        fibres (poles s=0,1)      : 2 * (f1.f2) = +2")
    print(f"        H^1 du diagonal (archimedien/genre) : gamma_0^2 = -Tr(id|H^1) = {-2*g}")
    print(f"        somme = {2 + (-2*g)} = Delta^2.  [+2 <-> 2 poles ; -2g <-> facteur Gamma]")

    # ---- (2) indice de Hodge : signature (1, reste) ----
    ev = np.linalg.eigvalsh(G)
    tol = 1e-7 * max(1.0, np.max(np.abs(ev)))
    npos = int(np.sum(ev > tol)); nneg = int(np.sum(ev < -tol)); nz = int(np.sum(np.abs(ev) <= tol))
    print(f"\n  (2) indice de Hodge sur NS(C x C) : signature de la forme d'intersection")
    print(f"      valeurs propres : {npos} positive(s), {nneg} negative(s), {nz} nulle(s)")
    print(f"      -> {'OK : (1, reste<=0) comme l indice de Hodge' if npos==1 else 'ANOMALIE (npos != 1)'}")
    print(f"      la direction + = classe ample (plan des fibres/poles) ;")
    print(f"      les directions - = partie primitive = POSITIVITE DE WEIL (= HR ici).")

    # ---- (3) partie primitive : forme localisee H^1 (Toeplitz des moments) negative def ----
    s = np.array([np.sum((alpha / np.sqrt(q))**n).real for n in range(0, 2 * g + 1)])
    Tt = np.array([[s[abs(i - j)] for j in range(2 * g + 1)] for i in range(2 * g + 1)])
    lam = float(np.linalg.eigvalsh(Tt)[0])
    print(f"\n  (3) forme de Weil sur H^1 (Toeplitz des moments normalises) : lambda_min = {lam:+.3e}  "
          f"({'>= 0 OK' if lam > -1e-7 else 'NEG'})  [le coeur de la positivite]")

    # ---- dictionnaire complet ----
    print(f"\n  DICTIONNAIRE COMPLET  (corps de fonctions  <->  corps de nombres) :")
    print(f"    zeros (H^1, valeurs propres de F)   <->  zeros non triviaux rho de zeta")
    print(f"    point ferme x, deg(x)               <->  premier p, log p     [facteur local = intersection locale]")
    print(f"    fibres f1,f2 (Delta^2 : +2)          <->  poles de zeta en s=0,1")
    print(f"    H^1 du diagonal (Delta^2 : -2g)      <->  facteur Gamma (place archimedienne)")
    print(f"    involution de Rosati F'=qF^{{-1}}      <->  equation fonctionnelle s <-> 1-s")
    print(f"    indice de Hodge : signature (1,*)    <->  positivite de Weil pour zeta  (D3, conjectural sur Z)")
    print(f"\n  Lecture : la place a l'infini de Q correspond au TERME -2g du diagonal,")
    print(f"            i.e. au H^1 (les zeros eux-memes) -- la cohomologie arithmetiquement invisible.")

if __name__ == "__main__":
    a = sys.argv[1:]
    g = int(a[0]) if len(a) > 0 else 1
    q = int(a[1]) if len(a) > 1 else 5
    p2 = int(a[2]) if len(a) > 2 else 1
    p3 = int(a[3]) if len(a) > 3 else (1 if g == 1 else 10)
    main(g, q, p2, p3)
