#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
weil_local.py — Forme quadratique de Weil reconstruite depuis le COTE ARITHMETIQUE
(places finies : premiers / coefficients de -Z'/Z  +  place archimedienne),
SANS jamais utiliser les zeros.

Etapes :
  1. coefficients c(k) de -Z'/Z par inversion de Dirichlet (verifies par identite exacte)
  2. recensement des zeros SUR la droite (changements de signe de xi reelle) -> pour VALIDER
  3. validation de la formule explicite : cote arithmetique == cote zeros  (zeta ET Epstein)
  4. matrice de Weil M_ij batie depuis l'arithmetique ; lambda_min :
        zeta    -> >= 0  (positivite de Weil OK, compatible RH)
        Epstein -> <  0  (positivite de Weil ECHOUE) -- obtenu SANS les zeros
  5. flot d'entropie (descente du quotient de Rayleigh) -> converge vers lambda_min

Requiert : ftsa_weil.py dans le meme repertoire (pour Zep) ; mpmath ; numpy.
Le recensement est mis en cache dans onz.json.
"""
import os, json, mpmath as mp, numpy as np
from math import log, sqrt
from ftsa_weil import Zep

mp.mp.dps = 30
np.seterr(all='ignore')
Ncond = mp.mpf(23)/4
LOGPI = float(mp.log(mp.pi)); LOGN = float(mp.log(Ncond))

# ----------------------------------------------------------- 1. coefficients
def build_coeffs(K=4000):
    Mlat=int(K**0.5)+3
    r=[0]*(K+1)
    for m in range(-Mlat,Mlat+1):
        for n in range(-Mlat,Mlat+1):
            k=m*m+m*n+6*n*n
            if 1<=k<=K: r[k]+=1
    divs=[[] for _ in range(K+1)]
    for d in range(1,K+1):
        for mlt in range(d,K+1,d): divs[mlt].append(d)
    b=[0.0]*(K+1); b[1]=1.0/r[1]
    for k in range(2,K+1):
        s=0.0
        for d in divs[k]:
            if d<k: s+=r[k//d]*b[d]
        b[k]=-s/r[1]
    c=[0.0]*(K+1)
    for k in range(1,K+1):
        s=0.0
        for d in divs[k]:
            if r[d]: s+=r[d]*log(d)*b[k//d]
        c[k]=s
    # check exact : (r * c)(k) == r(k) log k
    err=max(abs(sum(r[d]*c[k//d] for d in divs[k])-r[k]*log(k)) for k in range(2,300))
    return r,c,err

def vonmangoldt(P=100000):
    sieve=np.ones(P+1,bool); sieve[:2]=False
    for q in range(2,int(P**0.5)+1):
        if sieve[q]: sieve[q*q::q]=False
    Lam=np.zeros(P+1)
    for q in np.nonzero(sieve)[0]:
        pk=int(q); lq=log(int(q))
        while pk<=P: Lam[pk]=lq; pk*=int(q)
    return Lam

# ----------------------------------------------------------- 2. census (zeros sur la droite)
def xi_real(t):
    s=mp.mpf(1)/2+1j*mp.mpf(float(t)); dps=int(22+0.7*abs(t))
    return float(mp.re(Ncond**(s/2)*mp.pi**(-s)*mp.gamma(s)*Zep(s,dps=dps,KMAX=120)))

def census_online(Tmax=22.0, dt=0.1, refine=22, cache="onz.json"):
    if os.path.exists(cache):
        return json.load(open(cache))
    ts=np.arange(0.2,Tmax,dt); prev=pv=None; zeros=[]
    for t in ts:
        v=xi_real(t)
        if pv is not None and pv*v<0:
            lo,hi,plo=prev,t,pv
            for _ in range(refine):
                mid=(lo+hi)/2; vm=xi_real(mid)
                if plo*vm<0: hi=mid
                else: lo,plo=mid,vm
            zeros.append(float((lo+hi)/2))
        prev,pv=t,v
    json.dump(zeros,open(cache,"w"))
    return zeros

# le zero hors-droite (quadruplet en coordonnees gamma = (rho-1/2)/i : +-mu +- i delta)
MU,DELTA = 16.29021572039, 0.45326047479
QUAD = [MU-1j*DELTA, -MU-1j*DELTA, -MU+1j*DELTA, MU+1j*DELTA]

# ----------------------------------------------------------- 3. formule explicite, place archi
def psi_grids(R=30.0, hg=0.01):
    rs=np.arange(-R,R+hg,hg)
    pz=np.array([float(mp.re(mp.digamma(mp.mpf(1)/4+1j*mp.mpf(float(x))/2))) for x in rs])
    pe=np.array([float(mp.re(mp.digamma(mp.mpf(1)/2+1j*mp.mpf(float(x)))))   for x in rs])
    return rs,pz,pe
def trap(y,x): y=np.asarray(y); return float(np.sum((y[1:]+y[:-1])/2.0*np.diff(x)))

def validate(a, rs, pz, pe, Lam, c, K, onz):
    a=float(a); h=lambda x:np.exp(-a*x*x); g0=1/(2*sqrt(np.pi*a))
    g=lambda u:g0*np.exp(-(u*u)/(4*a)); hi2=np.exp(a/4)
    # ZETA
    arch=trap(h(rs)*(-LOGPI+pz),rs)/(2*np.pi)
    P=len(Lam)-1; nz=np.nonzero(Lam)[0]
    prim=2*float(np.sum(Lam[nz]/np.sqrt(nz)*g0*np.exp(-(np.log(nz))**2/(4*a))))
    Wz_a=2*hi2+arch-prim
    Wz_0=0.0;n=1
    while True:
        gm=float(mp.im(mp.zetazero(n))); t=np.exp(-a*gm*gm); Wz_0+=2*t
        if n>5 and t<1e-15: break
        n+=1
    # EPSTEIN
    arch_e=trap(h(rs)*(LOGN-2*LOGPI+2*pe),rs)/(2*np.pi)
    ks=np.arange(2,K+1); cc=np.array([c[k] for k in ks])
    prim_e=2*float(np.sum(cc/np.sqrt(ks)*g0*np.exp(-(np.log(ks))**2/(4*a))))
    We_a=2*hi2+arch_e-prim_e
    We_0=sum(2*np.exp(-a*x*x) for x in onz)+float(np.real(sum(np.exp(-a*z*z) for z in QUAD)))
    return (Wz_a,Wz_0,We_a,We_0)

# ----------------------------------------------------------- 4. matrice de Weil
def u_grid(rs,nu,a): return np.exp(-a*(rs-nu)**2)+np.exp(-a*(rs+nu)**2)
def u_ihalf(nu,a):  return 2*np.exp(-a*(nu*nu-0.25))*np.cos(a*nu)        # u(i/2), reel
def u_cplx(z,nu,a): return np.exp(-a*(z-nu)**2)+np.exp(-a*(z+nu)**2)     # u en point complexe
def Gk(ni,nj,a,lk):
    pref=(1/(2*np.pi))*sqrt(np.pi/(2*a))*np.exp(-(lk*lk)/(8*a))
    return pref*2*(np.exp(-(a/2)*(ni-nj)**2)*np.cos((ni+nj)*lk/2)
                  +np.exp(-(a/2)*(ni+nj)**2)*np.cos((ni-nj)*lk/2))

def M_arith(which, nus, a, rs, psi, coeff, klist):
    n=len(nus); M=np.zeros((n,n)); lk=np.log(klist); sk=np.sqrt(klist)
    U=[u_grid(rs,nu,a) for nu in nus]; uih=[u_ihalf(nu,a) for nu in nus]
    Phi=(-LOGPI+psi) if which=='zeta' else (LOGN-2*LOGPI+2*psi)
    for i in range(n):
        for j in range(i,n):
            pole=2*uih[i]*uih[j]
            arch=trap(U[i]*U[j]*Phi,rs)/(2*np.pi)
            coef=2*float(np.sum(coeff*Gk(nus[i],nus[j],a,lk)/sk))
            M[i,j]=M[j,i]=pole+arch-coef
    return M

def M_zero(gcoords, nus, a):   # cote zeros (pour recoupement) : M_ij = sum_rho u_i(g) u_j(g)
    n=len(nus); M=np.zeros((n,n),complex)
    for z in gcoords:
        v=np.array([u_cplx(z,nu,a) for nu in nus],complex); M+=np.outer(v,v)
    return M.real

# ----------------------------------------------------------- main
def main():
    print("="*74)
    K=4000
    r,c,err=build_coeffs(K)
    print("1. coefficients c(k) de -Z'/Z : identite r*c=r(k)logk  ->  erreur max = %.1e"%err)
    print("   (c(6)=%.3f c(10)=%.3f c(12)=%.3f : composites non nuls => pas de produit eulerien)"
          %(c[6],c[10],c[12]))
    Lam=vonmangoldt(100000)

    print("\n2. recensement des zeros SUR la droite (xi reelle, changements de signe)")
    onz=census_online()
    print("   %d zeros : %s"%(len(onz),[round(x,3) for x in onz]))

    print("\n3. VALIDATION formule explicite (cote arithmetique vs cote zeros)")
    rs,pz,pe=psi_grids()
    for a in (0.05,0.10,0.20):
        Wz_a,Wz_0,We_a,We_0=validate(a,rs,pz,pe,Lam,c,K,onz)
        print("   a=%.2f | ZETA arith=%.7f zeros=%.7f d=%.1e | EPST arith=%.7f zeros=%.7f d=%.1e"
              %(a,Wz_a,Wz_0,abs(Wz_a-Wz_0),We_a,We_0,abs(We_a-We_0)))

    print("\n4. SPECTRE de la forme de Weil bati depuis l'ARITHMETIQUE (base près de mu=16.29)")
    nus=[14.0+0.5*i for i in range(10)]; a=0.5
    ksZ=np.nonzero(Lam)[0]; coeffZ=Lam[ksZ]
    ksE=np.arange(2,K+1); coeffE=np.array([c[k] for k in ksE])
    MzA=M_arith('zeta',nus,a,rs,pz,coeffZ,ksZ.astype(float))
    MeA=M_arith('epstein',nus,a,rs,pe,coeffE,ksE.astype(float))
    lzA=np.linalg.eigvalsh(MzA); leA=np.linalg.eigvalsh(MeA)
    print("   [arith]  zeta    : lambda_min = %+.6e  -> %s"%(lzA[0],">=0 (Weil OK)" if lzA[0]>-1e-6 else "<0"))
    print("   [arith]  Epstein : lambda_min = %+.6e  -> %s"%(leA[0],"<0 (Weil ECHOUE)" if leA[0]<0 else ">=0"))
    # recoupement cote zeros (census complet + quadruplet hors-droite)
    gE=[g for G in onz for g in (G,-G)]+QUAD
    MeZ=M_zero(gE,nus,a); leZ=np.linalg.eigvalsh(MeZ)
    print("   [zeros ] Epstein : lambda_min = %+.6e   (recoupement)"%leZ[0])
    print("   ecart entrywise |M_arith - M_zero| (Epstein) = %.2e"%np.max(np.abs(MeA-MeZ)))

    print("\n5. FLOT D'ENTROPIE (descente du quotient de Rayleigh sur M_arith Epstein)")
    rng=np.random.default_rng(0); v=rng.standard_normal(len(nus)); v/=np.linalg.norm(v)
    lr=0.3/leA[-1]; tr=[]
    for _ in range(6000):
        Mv=MeA@v; R=float(v@Mv); tr.append(R); v=v-lr*2*(Mv-R*v); v/=np.linalg.norm(v)
    for s in (0,50,200,1000,5999): print("   t=%5d  R=%+.6e"%(s,tr[s]))
    print("   point critique = %+.6e  (= lambda_min arith = %+.6e)"%(tr[-1],leA[0]))
    print("="*74)

if __name__=="__main__":
    main()
