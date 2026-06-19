#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
weil_decomp.py — OU l'absence de multiplicativite injecte la negativite de Weil.

Experience propre : Epstein Z_Q0 et zeta_K = zeta*L(chi_-23) sont tous deux de
degre 2 sur le MEME corps Q(sqrt(-23)) -> MEME pole, MEME terme archimedien
(log(23/4)-2log(pi)+2Re psi(1/2+ir) = log23-2log(2pi)+2Re psi(1/2+ir)).
La SEULE difference est les coefficients de -L'/L :
   Epstein : c(k)              (non multiplicatif : c(6),c(12) != 0)
   zeta_K  : Lambda(n)(1+chi(n)) (porte par les puissances de premiers, Euler)
On compare lambda_min, on decompose, et on "remultiplicativise" l'Epstein.

Requiert weil_local.py (et ftsa_weil.py) dans le meme repertoire ; numpy.
"""
import numpy as np
from weil_local import (build_coeffs, vonmangoldt, psi_grids, trap,
                        u_grid, u_ihalf, Gk, LOGPI, LOGN)

K=4000
r,c,err=build_coeffs(K)
Lam=vonmangoldt(K)
rs,pz,pe=psi_grids()

QR={1,2,3,4,6,8,9,12,13,16,18}
def chi(n):
    n%=23
    return 0 if n==0 else (1 if n in QR else -1)
def ndp(k):                       # nombre de facteurs premiers distincts
    n=k;cnt=0;d=2
    while d*d<=n:
        if n%d==0:
            cnt+=1
            while n%d==0:n//=d
        d+=1
    if n>1:cnt+=1
    return cnt

ks=np.arange(2,K+1)
cE=np.array([c[k] for k in ks])                 # Epstein, non multiplicatif
cK=np.array([Lam[k]*(1+chi(k)) for k in ks])    # zeta_K, Euler (prime powers)
ispp=np.array([1.0 if ndp(k)==1 else 0.0 for k in ks])
cE_pp=cE*ispp; cE_co=cE*(1-ispp)                # decoupage Epstein
print("controle : coeffs zeta_K sur composites = %.1e (doit etre 0)"
      %np.max(np.abs(cK*(1-ispp))))

nus=[14.0+0.5*i for i in range(10)]; a=0.5
Phi = LOGN-2*LOGPI+2*pe          # archimedien degre 2, COMMUN a Epstein et zeta_K
lk=np.log(ks.astype(float)); sk=np.sqrt(ks.astype(float))
U=[u_grid(rs,nu,a) for nu in nus]; uih=[u_ihalf(nu,a) for nu in nus]

def comps(coeff):
    n=len(nus); POLE=np.zeros((n,n));ARCH=np.zeros((n,n));COEF=np.zeros((n,n))
    for i in range(n):
        for j in range(i,n):
            POLE[i,j]=POLE[j,i]=2*uih[i]*uih[j]
            ARCH[i,j]=ARCH[j,i]=trap(U[i]*U[j]*Phi,rs)/(2*np.pi)
            COEF[i,j]=COEF[j,i]=2*float(np.sum(coeff*Gk(nus[i],nus[j],a,lk)/sk))
    return POLE,ARCH,COEF

POLE,ARCH,COEF_E=comps(cE)
_,_,COEF_K =comps(cK)
_,_,COEF_pp=comps(cE_pp)
_,_,COEF_co=comps(cE_co)
M_E=POLE+ARCH-COEF_E
M_K=POLE+ARCH-COEF_K
M_mult=POLE+ARCH-COEF_pp

print("\nMEME pole, MEME archimedien (degre 2, meme corps) ; seuls les COEFFICIENTS different :")
print("  Epstein  c(k) non multiplicatif       : lambda_min = %+.6e"%np.linalg.eigvalsh(M_E)[0])
print("  zeta_K   Lambda(n)(1+chi), Euler       : lambda_min = %+.6e"%np.linalg.eigvalsh(M_K)[0])
print("  Epstein  termes COMPOSITES retires     : lambda_min = %+.6e"%np.linalg.eigvalsh(M_mult)[0])

w,V=np.linalg.eigh(M_E); v=V[:,0]
print("\ndecomposition de lambda_min(Epstein) au vecteur propre optimal :")
print("  pole               : %+.6e"%(v@POLE@v))
print("  archimedien        : %+.6e"%(v@ARCH@v))
print("  -coeff prime-power : %+.6e"%(-(v@COEF_pp@v)))
print("  -coeff COMPOSITE   : %+.6e   <== defaut de multiplicativite"%(-(v@COEF_co@v)))
print("  somme              : %+.6e"%(v@M_E@v))
