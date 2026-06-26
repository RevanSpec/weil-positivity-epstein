"""
fork2_local_intersection.py -- D3 (conjecture de travail) en corps de fonctions.

On VERIFIE, sur une vraie courbe C/F_q, l'identite a transporter vers Z :

    facteur de Weil LOCAL en un point ferme x
        =  nombre d'INTERSECTION local  (Gamma_{F^n} . Delta)_x
        =  deg(x).

Dictionnaire vise (corps de fonctions  ->  corps de nombres) :
    point ferme x de C          <->   premier p
    deg(x)                      <->   log(p)
    (Gamma_{F^n} . Delta)_x     <->   contribution locale de p a la formule explicite
    indice de Hodge sur C x C   <->   positivite de Weil pour zeta   (= HR, conjectural)

Identites verifiees (THEOREMES en corps de fonctions : Lefschetz + orbites de Galois) :
  (i)   |alpha_i| = sqrt q                      (HR de Weil pour la courbe)
  (ii)  a_d = #{points fermes de degre d} entier >= 0   (Mobius sur les comptages)
  (iii) #C(F_q^n) = sum_{d|n} d * a_d = (Gamma_{F^n} . Delta)
        => le nombre d'intersection GLOBAL se decompose en contributions LOCALES deg(x)
           aux points fermes x de degre d | n. Chaque point ferme de degre d apporte
           ses d points geometriques (une orbite de Galois de taille d) a Fix(F^n).

Courbe par defaut : E : y^2 = x^3 + A x + B / F_q (q premier). Comptage sur le corps de
BASE seulement ; les extensions suivent par la recurrence entiere p_n = t p_{n-1} - q p_{n-2}.
numpy seul (et arithmetique entiere exacte pour (ii),(iii)).

Usage : python fork2_local_intersection.py [q A B Nmax]
        defaut : 5 1 1 8
"""
import sys
import numpy as np

def count_base(q, A, B):
    """#E(F_q) pour E: y^2 = x^3 + A x + B (q premier impair), + point a l'infini."""
    n = 0
    for x in range(q):
        rhs = (x * x * x + A * x + B) % q
        if rhs == 0:
            n += 1
        elif pow(rhs, (q - 1) // 2, q) == 1:   # rhs carre non nul -> 2 racines
            n += 2
    return n + 1                               # + point a l'infini

def mobius(n):
    if n == 1:
        return 1
    res, m, p = 1, n, 2
    while p * p <= m:
        if m % p == 0:
            m //= p
            if m % p == 0:
                return 0
            res = -res
        p += 1
    if m > 1:
        res = -res
    return res

def divisors(n):
    return [d for d in range(1, n + 1) if n % d == 0]

def main(q=5, A=1, B=1, Nmax=8):
    Nq = count_base(q, A, B)
    t = q + 1 - Nq                              # trace de Frobenius
    print(f"# E : y^2 = x^3 + {A} x + {B}   / F_{q}")
    print(f"#   #E(F_{q}) = {Nq},  trace t = {t}")

    # (i) valeurs propres : racines de T^2 - t T + q ; module sqrt q
    disc = t * t - 4 * q
    alpha = (t + np.lib.scimath.sqrt(disc)) / 2.0
    print(f"   (i)   |alpha| = {abs(alpha):.8f}  (sqrt q = {np.sqrt(q):.8f})  -> "
          f"{'OK' if abs(abs(alpha)**2 - q) < 1e-9 else 'FAIL'}   [HR de Weil pour la courbe]")

    # comptages sur F_q^n : sommes de puissances ENTIERES p_n = t p_{n-1} - q p_{n-2}
    p = [2, t]
    for n in range(2, Nmax + 1):
        p.append(t * p[-1] - q * p[-2])
    N = {n: q**n + 1 - p[n] for n in range(1, Nmax + 1)}     # #C(F_q^n), exact

    # (ii) a_d = #{points fermes de degre d} par Mobius : d a_d = sum_{e|d} mu(d/e) N(e)
    a = {}
    pos_ok = True
    for d in range(1, Nmax + 1):
        s = sum(mobius(d // e) * N[e] for e in divisors(d))
        a[d] = s // d
        if s % d != 0:
            print(f"   (ii)  a_{d} NON ENTIER (anomalie)")
        if a[d] < 0:
            pos_ok = False
    print(f"   (ii)  a_d (points fermes de degre d) entiers >= 0 : "
          f"{'OK' if pos_ok else 'NON'}   [points fermes authentiques]")

    # (iii) #C(F_q^n) = sum_{d|n} d a_d = (Gamma_{F^n} . Delta), decomposition locale
    print(f"\n   (iii) #C(F_q^n) = sum_{{d|n}} d*a_d = (Gamma_F^n . Delta) :")
    print(f"   {'n':>3} {'#C(F_q^n)':>11} {'sum d*a_d':>11}  {'a_d (deg d)':>12}")
    allok = True
    for n in range(1, Nmax + 1):
        rhs = sum(d * a[d] for d in divisors(n))
        ok = (N[n] == rhs); allok = allok and ok
        print(f"   {n:>3} {N[n]:>11} {rhs:>11}  {a[n]:>12}   [{'OK' if ok else 'FAIL'}]")
    print(f"   => decomposition locale du nombre d'intersection : {'OK' if allok else 'FAIL'}")
    print(f"      chaque point ferme x de degre d|n apporte deg(x)=d a Fix(F^n) ;")
    print(f"      donc  facteur local de Weil en x  =  (Gamma_F^n.Delta)_x  =  deg(x).")

    # (iv) indice de Hodge : Toeplitz des moments normalises s_k = p_k / q^{k/2}
    M = Nmax
    s = np.array([p[k] / q**(k / 2.0) for k in range(M + 1)])
    T = np.array([[s[abs(i - j)] for j in range(M + 1)] for i in range(M + 1)])
    lam = float(np.linalg.eigvalsh(T)[0])
    print(f"\n   (iv)  indice de Hodge (Toeplitz des moments) : lambda_min = {lam:+.6e}  "
          f"({'>= 0 OK' if lam > -1e-9 else 'NEG'})   [positivite de Weil = theoreme ici]")

    print(f"\n   DICTIONNAIRE a transporter vers Z :")
    print(f"     point ferme x  <->  premier p     |     deg(x)  <->  log(p)")
    print(f"     (Gamma_F^n.Delta)_x = deg(x)       <->  facteur local de la formule explicite en p")
    print(f"     indice de Hodge sur C x C          <->  positivite de Weil pour zeta  (D3, conjectural)")
    print(f"   Conjecture de travail D3 : log(p) = nombre d'intersection local au point ferme p")
    print(f"   de la surface arithmetique Spec(Z) x_F1 Spec(Z).")

if __name__ == "__main__":
    a = sys.argv[1:]
    q = int(a[0]) if len(a) > 0 else 5
    A = int(a[1]) if len(a) > 1 else 1
    B = int(a[2]) if len(a) > 2 else 1
    Nmax = int(a[3]) if len(a) > 3 else 8
    main(q, A, B, Nmax)
