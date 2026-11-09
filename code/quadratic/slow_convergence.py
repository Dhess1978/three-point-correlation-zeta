#!/usr/bin/env python3
"""Slow convergence (Section 5): one-level density of quadratic Dirichlet L-functions in [0, 1/2] from the ratios conjecture as a function of
l = ln(|D|/pi), with the same formula as compute_families.py of study III (Conrey-Snaith 2007, Theorem 3.1; for D < 0 with Gamma(3/4 + it/2)),
but with the integral of c(t) up to t(x = 1/2) by Gauss-Legendre quadrature (64 nodes), so that it reaches |D| ~ 1e120 (the grid DT = 0.002 of study III does not).
  density_k(l) = 1 + 2 * int_0^{t_half} c_k(t) dt,  c_k(t) = Re(zeta'/zeta(1 + 2it) + A2(t) - exp(-itl) * Gamma ratio_k * zeta(1 - 2it) * A1(t)) / pi,  theta_k(t_half)/pi = 1/2;
  even (k = 0) and odd (k = 1) characters are weighted 1 : 1 (in a contiguous stretch positive and negative fundamental discriminants are equally frequent).
Check: at the l of the families B10 ... G15 against density_0_05_ratios of study III (mean over the D of the family, grid DT) - difference below 2e-3.
usage: python3 -I slow_convergence.py TRENDS.json OUTPUT.json   (TRENDS.json: reference_g15_trends.json with the families B10 ... G15)"""
import sys, json, math, numpy as np, mpmath as mp, scipy.special as sp
TR, OUT = sys.argv[1:3]
Mp = 200000; s = np.ones(Mp + 1, bool); s[:2] = False
for p in range(2, int(Mp**0.5) + 1):
    if s[p]: s[p*p::p] = False
P = np.nonzero(s)[0].astype(float); lp = np.log(P)
mp.mp.dps = 30
XG, WG = np.polynomial.legendre.leggauss(64)
def c_k(t, l, k):
    z = np.exp(2j*t*lp); A1 = np.exp(np.sum(np.log(1 - (z - 1)/(P*P - 1)))); A2 = np.sum(lp/((P + 1)*(P*z - 1)))
    s1 = mp.mpc(1, 2*t); ZL = complex(mp.zeta(s1, derivative=1)/mp.zeta(s1)); ZR = complex(mp.zeta(mp.mpc(1, -2*t)))
    GR = complex(np.exp(sp.loggamma(0.25 + 0.5*k - 0.5j*t) - sp.loggamma(0.25 + 0.5*k + 0.5j*t)))
    return ((ZL + A2) - np.exp(-1j*t*l)*GR*ZR*A1).real/math.pi
def t_pul(l, k):
    th = lambda t: (0.5*t*l + float(np.imag(sp.loggamma(0.25 + 0.5*k + 0.5j*t))))/math.pi - 0.5
    lo, hi = 0.0, 50.0
    for _ in range(100):
        m = 0.5*(lo + hi)
        if th(m) < 0: lo = m
        else: hi = m
    return 0.5*(lo + hi)
def density(l, k):
    T = t_pul(l, k); t = 0.5*T*(XG + 1); w = 0.5*T*WG
    return 1 + 2*float(sum(wi*c_k(ti, l, k) for ti, wi in zip(t, w)))
USP = 1 - float(mp.quad(lambda x: mp.sin(2*mp.pi*x)/(2*mp.pi*x), [0, 0.5]))/0.5   # USp density in [0, 1/2]: 2 * int_0^1/2 (1 - sinc 2x) dx
T = json.load(open(TR)); D = T['density']
check = {}
for f, q in D.items():
    h = 0.5*(density(q['l'], 0) + density(q['l'], 1)); check[f] = dict(l=q['l'], curve=h, study_iii=q['density_0_05_ratios'], difference=h - q['density_0_05_ratios'])
    print(f, round(q['l'], 3), round(h, 5), q['density_0_05_ratios'], flush=True)
LG = sorted(set([round(v, 3) for v in np.concatenate([np.arange(5, 30, 0.5), np.arange(30, 80, 2.5), np.arange(80, 300, 10)])]))
curve = []
for l in LG:
    h0, h1 = density(l, 0), density(l, 1); curve.append([l, 0.5*(h0 + h1), h0, h1])
curve_l = [(r[0], r[1]) for r in curve]
def l_below(target):   # smallest l with the excess below target (linear interpolation on the grid)
    for (l0, a), (l1, b) in zip(curve_l, curve_l[1:]):
        if a - USP >= target > b - USP: return l0 + (a - USP - target)/((a - USP) - (b - USP))*(l1 - l0)
    return None
below = {}
for target in (0.1, 0.05, 0.02, 0.01):
    l = l_below(target); below[str(target)] = dict(l=l, log10D=(l + math.log(math.pi))/math.log(10) if l else None)
o = dict(usp=USP, nodes=64, primes_up_to=Mp, check=check, check_max=max(abs(v['difference']) for v in check.values()), curve=curve, below=below)
json.dump(o, open(OUT, 'w'), indent=1)
print('USp', USP, 'check max', o['check_max']); print(json.dumps(below))
