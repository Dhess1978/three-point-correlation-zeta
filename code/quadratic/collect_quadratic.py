#!/usr/bin/env python3
"""Section 5 of study IV in one file: Table 5 (families of study III, the three bands at |D| = 1e12 and their combination), the fits of the excess over USp
(c/l, a + c/l, a power of |D|, slope of l x excess against l), the curve of the ratios conjecture and its thresholds, and the three bands against the controls.
usage: python3 -I collect_quadratic.py TRENDS.json FAMILIES_G15bc.json THREE_BANDS.json SLOW_CONVERGENCE.json OUTPUT.json
  TRENDS.json: reference_g15_trends.json (study III scripts, families B10 ... G15); FAMILIES_G15bc.json: reference_families_G15b_G15c.json"""
import sys, json, math, numpy as np
from scipy import stats
TR, FB, TB, SC, OUT = sys.argv[1:6]
T = json.load(open(TR)); F = json.load(open(FB))['families']; B = json.load(open(TB)); C = json.load(open(SC))
D = T['density']; FAM = ['B10', 'C100', 'D100', 'E100', 'G15']
tab = {k: dict(l=D[k]['l'], n=D[k]['n'], density=D[k]['density_0_05'], se=D[k]['density_0_05_se'], ratios=D[k]['density_0_05_ratios'], usp=D[k]['density_0_05_usp']) for k in FAM}
for f in ('G15b', 'G15c'):
    A = F[f]['all']; tab[f] = dict(l=A['mean_ln_D'] - math.log(math.pi), n=A['n'], density=A['density_0_05'], se=A['density_0_05_se'], ratios=A['density_0_05_ratios'], usp=A['density_0_05_usp'])
H = B['summary']['density']; lG = float(np.mean([B['bands'][f]['l'] for f in ('G15', 'G15b', 'G15c')]))
combined = dict(l=lG, n=int(sum(B['bands'][f]['n'] for f in ('G15', 'G15b', 'G15c'))), density=H['mean'], se=H['se'], ratios=H['ratios'], usp=H['usp'])
for v in list(tab.values()) + [combined]: v.update(l_excess=v['l']*(v['density'] - v['usp']), l_excess_se=v['l']*v['se'], l_excess_ratios=v['l']*(v['ratios'] - v['usp']))
def fits(keys):
    l = np.array([tab[k]['l'] for k in keys]); e = np.array([tab[k]['density'] - tab[k]['usp'] for k in keys]); s = np.array([tab[k]['se'] for k in keys]); x = l + math.log(math.pi); n = len(keys)
    w = (e/s)**2; A = np.vstack([np.ones_like(x), x]).T; Cm = np.linalg.inv(A.T @ (A*w[:, None])); b = Cm @ A.T @ (w*np.log(e))
    wc = 1/s**2; c = float(np.sum(wc*e/l)/np.sum(wc/l**2)); c_se = float(1/math.sqrt(np.sum(wc/l**2)))
    Q = np.vstack([np.ones_like(l), 1/l]).T; Cq = np.linalg.inv(Q.T @ (Q*wc[:, None])); q = Cq @ Q.T @ (wc*e)
    chi_c = float(np.sum(((e - c/l)/s)**2)); chi_m = float(np.sum(((e - np.exp(b[0] + b[1]*x))/s)**2)); chi_q = float(np.sum(((e - q[0] - q[1]/l)/s)**2)); l120 = 120*math.log(10) - math.log(math.pi)
    return dict(families=keys, c_over_l=dict(c=c, se=c_se, chi2=chi_c, dof=n - 1, p=float(stats.chi2.sf(chi_c, n - 1))),
                power=dict(exponent=float(b[1]), se=float(math.sqrt(Cm[1, 1])), log_amplitude=float(b[0]), chi2=chi_m, dof=n - 2, p=float(stats.chi2.sf(chi_m, n - 2))),
                free_limit=dict(a=float(q[0]), a_se=float(math.sqrt(Cq[0, 0])), c=float(q[1]), c_se=float(math.sqrt(Cq[1, 1])), chi2=chi_q, dof=n - 2, p=float(stats.chi2.sf(chi_q, n - 2))),
                at_1e120=dict(l=l120, c_over_l=c/l120, power=float(math.exp(b[0] + b[1]*120*math.log(10)))))
o = dict(table=tab, combined_1e12=combined, fits_5=fits(FAM), fits_7=fits(FAM[:4] + ['G15', 'G15b', 'G15c']), usp=C['usp'], ratios_curve=C['curve'], thresholds=C['below'], curve_check=C['check'],
         three_bands=B)
json.dump(o, open(OUT, 'w'), indent=1)
rows = ['family;l;log10_D;n;density_0_half;se;ratios;usp;l_excess;l_excess_se;l_excess_ratios']
for k, v in list(tab.items()) + [('G15-G15c combined', combined)]:
    rows.append(';'.join([k] + ['%.10g' % x for x in (v['l'], (v['l'] + math.log(math.pi))/math.log(10))] + [str(v['n'])] + ['%.10g' % v[q] for q in ('density', 'se', 'ratios', 'usp', 'l_excess', 'l_excess_se', 'l_excess_ratios')]))
open(OUT.replace('.json', '') + '_table5.csv', 'w').write('\n'.join(rows) + '\n')
f7 = o['fits_7']; print('c/l: c = %.3f +- %.3f, chi2 %.1f/%d (p %.2f); power %.3f +- %.3f, chi2 %.1f/%d (p %.4f); a = %.4f +- %.4f; thresholds 0.05 at 10^%.2f, 0.01 at 10^%.2f' % (
    f7['c_over_l']['c'], f7['c_over_l']['se'], f7['c_over_l']['chi2'], f7['c_over_l']['dof'], f7['c_over_l']['p'], f7['power']['exponent'], f7['power']['se'], f7['power']['chi2'], f7['power']['dof'], f7['power']['p'],
    f7['free_limit']['a'], f7['free_limit']['a_se'], C['below']['0.05']['log10D'], C['below']['0.01']['log10D']))
