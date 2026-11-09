#!/usr/bin/env python3
"""Independent check of Section 5 (shares no code with slow_convergence.py, three_bands.py and collect_quadratic.py):
(1) the one-level density of the ratios conjecture in [0, 1/2] by adaptive quadrature (mpmath.quad, theta root by mpmath.findroot) at the l of B10 and G15
    and at the thresholds 0.05 and 0.01 of the curve; USp in [0, 1/2] = 1 - Si(pi)/pi;
(2) the combination of the three bands (weighted mean, chi^2 among the bands, p by mpmath.gammainc) and Welch's t with p by mpmath.betainc;
(3) the fits by numpy.polyfit with weights.
usage: python3 -I quadratic_check.py reference_quadratic.json OUTPUT.json"""
import sys, json, math, numpy as np, mpmath as mp
R = json.load(open(sys.argv[1])); o = {}
mp.mp.dps = 20
Mp = 200000; s = np.ones(Mp + 1, bool); s[:2] = False
for p in range(2, int(Mp**0.5) + 1):
    if s[p]: s[p*p::p] = False
P = np.nonzero(s)[0].astype(float); lp = np.log(P)
def c_k(t, l, k):
    t = float(t); z = np.exp(2j*t*lp); A1 = complex(np.prod(1 - (z - 1)/(P*P - 1))); A2 = complex(np.sum(lp/((P + 1)*(P*z - 1))))
    s1 = mp.mpc(1, 2*t); a = mp.zeta(s1, derivative=1)/mp.zeta(s1) + A2
    g = mp.exp(mp.loggamma(mp.mpc(0.25 + 0.5*k, -0.5*t)) - mp.loggamma(mp.mpc(0.25 + 0.5*k, 0.5*t)))
    return mp.re(a - mp.exp(mp.mpc(0, -t*l))*g*mp.zeta(mp.mpc(1, -2*t))*A1)/mp.pi
def density(l):
    v = []
    for k in (0, 1):
        T = mp.findroot(lambda t: (0.5*t*l + mp.im(mp.loggamma(mp.mpc(0.25 + 0.5*k, t/2))))/mp.pi - 0.5, 1.5/l)
        v.append(1 + 2*float(mp.quad(lambda t: c_k(t, l, k), [0, T/2, T])))
    return 0.5*(v[0] + v[1])
USP = float(1 - mp.si(mp.pi)/mp.pi); o['usp'] = USP; o['d_usp'] = abs(USP - R['usp'])
K = np.array(R['ratios_curve']); o['points'] = {}
for nm in ('B10', 'G15'):
    l = R['table'][nm]['l']; h = density(l); o['points'][nm] = dict(l=l, own=h, curve=float(np.interp(l, K[:, 0], K[:, 1])), study_iii=R['table'][nm]['ratios'])
for c in ('0.05', '0.01'):
    l = R['thresholds'][c]['l']; o['points']['below_' + c] = dict(l=l, own_excess=density(l) - USP, target=float(c))
B = R['three_bands']; Bd = B['bands']; S = B['summary']; bands = ['G15', 'G15b', 'G15c']
def comb(v, e):
    v, e = np.array(v), np.array(e); w = 1/e**2; m = float(np.sum(w*v)/np.sum(w)); ch = float(np.sum(((v - m)/e)**2)); return m, float(1/math.sqrt(np.sum(w))), ch, float(mp.gammainc(1, ch/2, mp.inf, regularized=True))
m, e, ch, p = comb([Bd[f]['density'] for f in bands], [Bd[f]['density_se'] for f in bands]); o['d_combination'] = max(abs(m - S['density']['mean']), abs(e - S['density']['se']), abs(ch - S['density']['chi2']), abs(p - S['density']['p']))
dw = 0
for k, P_ in S['band_statistics'].items():
    a, c = np.array(P_['bands'], float), np.array(P_['controls'], float); va, vc = a.var(ddof=1)/3, c.var(ddof=1)/3; t = (a.mean() - c.mean())/math.sqrt(va + vc); df = (va + vc)**2/(va**2/2 + vc**2/2)
    pw = float(mp.betainc(df/2, 0.5, 0, df/(df + t*t), regularized=True)); dw = max(dw, abs(t - P_['welch_t']), abs(pw - P_['welch_p']))
o['d_welch'] = dw
T = R['table']; F7 = R['fits_7']; keys = F7['families']
l = np.array([T[k]['l'] for k in keys]); ex = np.array([T[k]['density'] - T[k]['usp'] for k in keys]); sd = np.array([T[k]['se'] for k in keys])
b = np.polyfit(l + math.log(math.pi), np.log(ex), 1, w=ex/sd); q = np.polyfit(1/l, ex, 1, w=1/sd); c = float(np.sum(ex/l/sd**2)/np.sum(1/(l*l*sd*sd)))
o['d_fits'] = max(abs(b[0] - F7['power']['exponent']), abs(c - F7['c_over_l']['c']), abs(q[1] - F7['free_limit']['a']), abs(q[0] - F7['free_limit']['c']))
pt = o['points']; o['ok'] = bool(o['d_usp'] < 1e-12 and all(abs(pt[n]['own'] - pt[n]['curve']) < 3e-4 and abs(pt[n]['own'] - pt[n]['study_iii']) < 2e-3 for n in ('B10', 'G15'))
                                and all(abs(pt['below_' + c]['own_excess'] - float(c)) < 2e-4 for c in ('0.05', '0.01')) and o['d_combination'] < 1e-9 and o['d_welch'] < 1e-9 and o['d_fits'] < 1e-9)
json.dump(o, open(sys.argv[2], 'w'), indent=1); print('ok' if o['ok'] else 'DIFFERENCE', json.dumps({k: v for k, v in o.items() if k != 'points'}), json.dumps(pt))
