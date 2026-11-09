#!/usr/bin/env python3
"""
Low-lying zeros of quadratic Dirichlet L-functions, follow-up of study III (study IV): the family G15 (all fundamental discriminants
10^12 < |D| <= 10^12 + 15 000, about 9 100 L-functions; "much smaller sets of zeros much further up", S. J. Miller) set into the trends over the
conductor of Sections 4.3-4.8, as a prediction test: the fits of Section 4.8 (with E100, without G15) predict G15, then the fits are repeated with G15.

usage:  python3 -I g15_trends.py DATA_DIR G15_DIR OUT_DIR      (a few seconds; writes reference_g15_trends.json)
        DATA_DIR = data of study III (reference_*.json and reference_*_E100.json); G15_DIR = reference_families_G15.json, reference_two_level_G15.json
        and, if present, reference_recomputation_G15.json and reference_explicit_weights_G15.json (from the scripts with --families G15)
        DATA_DIR holds reference_families.json, reference_two_level.json, reference_recomputation.json, reference_explicit_weights.json and their
        counterparts for E100 (..._E100.json from compute_families.py, two_level_first_zero.py, layered_recomputation.py and explicit_weights.py
        with --families E100)

(1) the shift of the mean first zero of even characters against USp(2N_eff) (N_eff on [0, 1/2], Section 4.5) as c/l, c/l^2 and c l^-a with free a,
    for the five families of Section 4.5 (B10, C100, R10, D100, E10) and with E100 in place of E10, and the deviation of E100 from the fits without it;
(2) the same for the deviation of the shape, dF(1) = F(1) - USp(2N_eff) at x_1 = 1 (Section 4.7, two-level part);
(3) the exact ratio of one stretch over all starting positions to random samples, the reduction by D mod 840 and the ratio of the exact reduction
    by the stratifications to the line through two points (Section 4.7), and the slope of the explicit formula with and without the dual term.
Weighted least squares with the statistical errors; l = ln(|D|/pi) of the family.
"""
import sys, os, json, math, numpy as np
from scipy.stats import chi2 as CHI2
DATA, GD, OUT_DIR = sys.argv[1], sys.argv[2], sys.argv[3]; os.makedirs(OUT_DIR, exist_ok=True)
ld = lambda n: json.load(open(os.path.join(DATA, n)))
lg = lambda n: json.load(open(os.path.join(GD, n)))
has = lambda n: os.path.exists(os.path.join(GD, n))
FAM = ld('reference_families.json')['families']; FAM.update(ld('reference_families_E100.json')['families']); FAM.update(lg('reference_families_G15.json')['families'])
TL = ld('reference_two_level.json')['families']; TL.update(ld('reference_two_level_E100.json')['families']); TL.update(lg('reference_two_level_G15.json')['families'])
RC = ld('reference_recomputation.json'); RE = ld('reference_recomputation_E100.json'); RG = lg('reference_recomputation_G15.json') if has('reference_recomputation_G15.json') else None
EW = ld('reference_explicit_weights.json')['families']; EW.update(ld('reference_explicit_weights_E100.json')['families'])
if has('reference_explicit_weights_G15.json'): EW.update(lg('reference_explicit_weights_G15.json')['families'])
def l_of(f): return FAM[f]['plus']['mean_ln_D'] - math.log(math.pi)
def fit_power(l, y, s, a):
    x = l**-a; c = float(np.sum(x*y/s**2)/np.sum(x*x/s**2)); return c, float(1/math.sqrt(np.sum(x*x/s**2))), float(np.sum(((y - c*x)/s)**2))
def fits(rows):
    l = np.array([r[1] for r in rows]); y = np.array([r[2] for r in rows]); s = np.array([r[3] for r in rows]); out = {}
    for a in (1, 2):
        c, ce, ch = fit_power(l, y, s, a); out['a%d' % a] = dict(c=c, se=ce, chi2=ch, dof=len(rows) - 1, p=float(CHI2.sf(ch, len(rows) - 1)))
    grid = np.round(np.arange(0.2, 6.0001, 0.005), 3); chs = np.array([fit_power(l, y, s, a)[2] for a in grid]); i = int(np.argmin(chs)); ok = grid[chs <= chs[i] + 1]
    c, ce, ch = fit_power(l, y, s, grid[i])
    out['free'] = dict(a=float(grid[i]), a_lo=float(ok.min()), a_hi=float(ok.max()), c=c, chi2=ch, dof=len(rows) - 2, p=float(CHI2.sf(ch, max(len(rows) - 2, 1))))
    return out
def pred(F, l): return {k: v['c']*l**-(1 if k == 'a1' else 2 if k == 'a2' else v['a']) for k, v in F.items()}
OUT = dict(study='III', follow_up='G15 (study IV)', families_base=['B10', 'C100', 'R10', 'D100', 'E100'], families_new=['B10', 'C100', 'R10', 'D100', 'E100', 'G15'])
NEW = 'G15'
# ---- (1) mean first zero, even and odd, N_eff on [0, 1/2]
sh = {}
for f in ('B10', 'C10', 'C100', 'R10', 'D10', 'D100', 'E10', 'E100', NEW):
    sh[f] = dict(l=l_of(f), **{sg: dict(d=FAM[f][sg]['x1_mean'] - FAM[f][sg]['x1_usp2n_eff'], se=FAM[f][sg]['x1_se'], ks=FAM[f][sg]['ks_usp2n_eff'], ks_critical=FAM[f][sg]['ks_critical'], n=FAM[f][sg]['n'])
                                for sg in ('plus', 'minus')})
OUT['shift'] = sh
R0 = [(f, sh[f]['l'], sh[f]['plus']['d'], sh[f]['plus']['se']) for f in OUT['families_base']]; R1 = R0 + [(NEW, sh[NEW]['l'], sh[NEW]['plus']['d'], sh[NEW]['plus']['se'])]
F0 = fits(R0); F1 = fits(R1); lN = sh[NEW]['l']; P0 = pred(F0, lN)
OUT['shift_fit'] = dict(without_G15=F0, with_G15=F1, G15=dict(value=sh[NEW]['plus']['d'], se=sh[NEW]['plus']['se'], predicted_without_G15=P0, z={k: (sh[NEW]['plus']['d'] - v)/sh[NEW]['plus']['se'] for k, v in P0.items()}),
                        odd_G15=dict(value=sh[NEW]['minus']['d'], se=sh[NEW]['minus']['se']))
# ---- (2) shape: dF(1) of even characters (two-level part)
dF = {f: dict(l=l_of(f), plus=[TL[f]['plus']['rows']['1.0']['dF'], TL[f]['plus']['rows']['1.0']['F'][1]], minus=[TL[f]['minus']['rows']['1.0']['dF'], TL[f]['minus']['rows']['1.0']['F'][1]]) for f in TL}
OUT['dF1'] = dF
S0 = [(f, dF[f]['l'], dF[f]['plus'][0], dF[f]['plus'][1]) for f in ('B10', 'C100', 'R10', 'D100', 'E100')]; S1 = S0 + [(NEW, dF[NEW]['l'], dF[NEW]['plus'][0], dF[NEW]['plus'][1])]
G0 = fits(S0); G1 = fits(S1); Q0 = pred(G0, lN)
OUT['dF1_fit'] = dict(without_G15=G0, with_G15=G1, G15=dict(value=dF[NEW]['plus'][0], se=dF[NEW]['plus'][1], predicted_without_G15=Q0, z={k: (dF[NEW]['plus'][0] - v)/dF[NEW]['plus'][1] for k, v in Q0.items()}))
OUT['calibration_G15'] = dict(plus=TL[NEW]['plus'].get('calibration'), minus=TL[NEW]['minus'].get('calibration'), even_minus_odd=TL[NEW].get('even_minus_odd_by_calibration'))
# ---- (3) density: [0, 1/2] above the USp limit, chi^2 against the ratios conjecture and the null models, for all families
OUT['density'] = {f: {k: FAM[f]['all'][k] for k in ('n', 'chi2_005_ratios', 'chi2_005_without_AD', 'chi2_005_usp2n', 'chi2_005_usp', 'density_0_05', 'density_0_05_se', 'density_0_05_ratios', 'density_0_05_usp', 'x1_mean', 'x1_se', 'x1_usp2n_eff', 'N_eff')}
                  for f in ('B10', 'C100', 'R10', 'D100', 'E100', NEW) if f in FAM}
for f in OUT['density']: OUT['density'][f]['l'] = l_of(f); d = OUT['density'][f]; d['excess_over_usp_z'] = (d['density_0_05'] - d['density_0_05_usp'])/d['density_0_05_se']; d['ratios_z'] = (d['density_0_05'] - d['density_0_05_ratios'])/d['density_0_05_se']
# ---- (4) consecutive stretch and explicit formula (only if the G15 counterparts exist; the stretch needs at least 6 074 L)
cons = {}
for f, J in (('B10', RC), ('C100', RC), ('D100', RC), ('E100', RE)) + (((NEW, RG),) if RG else ()):
    try:
        F = J['families'][f]; D = F['designs']; r = D['random']['estimate']['low']; m8 = D['mod840']['estimate']['low']
        cons[f] = dict(l=F['l'], n=F['n'], ratio_stretch_random=J['z_44'][f]['stretch_all_starts_ratio'], reduction_mod840=(r[0] - m8[0])/r[0], z_mod840=J['z_44'][f]['mod840'],
                       slope_diagonal=EW[f]['regression']['slope_chi'] if f in EW else None, slope_dual=EW[f]['regression']['slope_chi_with_dual'] if f in EW else None,
                       stretch_explained_by_primes=EW[f]['stretch'].get('reduction_sum') if f in EW else None, stretch_reduction_exact=EW[f]['stretch'].get('reduction_exact') if f in EW else None)
    except Exception as e: cons[f] = dict(error=str(e)[:200])
OUT['conductor'] = cons
json.dump(OUT, open(os.path.join(OUT_DIR, 'reference_g15_trends.json'), 'w'), indent=1)
print('shift of even characters: G15 %+.4f +- %.4f; predicted from the fits with E100 (1/l, 1/l^2, free): %s; z %s' % (sh[NEW]['plus']['d'], sh[NEW]['plus']['se'],
      ', '.join('%.4f' % v for v in P0.values()), ', '.join('%.1f' % v for v in OUT['shift_fit']['G15']['z'].values())))
print('odd characters: G15 %+.4f +- %.4f' % (sh[NEW]['minus']['d'], sh[NEW]['minus']['se']))
print('dF(1) of even characters: G15 %+.4f +- %.4f; predicted %s' % (dF[NEW]['plus'][0], dF[NEW]['plus'][1], ', '.join('%.4f' % v for v in Q0.values())))
for f, d in OUT['density'].items():
    print('%-5s l %.2f n %6d chi2 ratios %.3f usp %.1f | [0,1/2] %.4f +- %.4f (ratios %.4f, USp %.4f; %.1f errors above USp, %.1f from ratios) | N_eff %.2f' % (
          f, d['l'], d['n'], d['chi2_005_ratios'], d['chi2_005_usp'], d['density_0_05'], d['density_0_05_se'], d['density_0_05_ratios'], d['density_0_05_usp'], d['excess_over_usp_z'], d['ratios_z'], d['N_eff']))
for f, c in cons.items(): print(f, c)
