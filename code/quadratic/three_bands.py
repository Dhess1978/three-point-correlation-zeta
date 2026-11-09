#!/usr/bin/env python3
"""Three bands of the same width at |D| = 1e12 (G15, G15b, G15c; about 9,100 L-functions each) against three control stretches of the same width at 1e10 (E15a, E15b, E15c).
(1) low zeros: density in [0, 1/2], chi^2 against the ratios conjecture, first zero and shift of the even characters by bands and combined (weighted mean);
    combined shift against the predictions of the fits 1/l, 1/l^2 and a free power without G15 (study III, Section 4.8; g15_trends.py);
(2) explicit formula: slope with the dual term by bands;
(3) band statistics (stretch/random over all starts, mod 840 z, share of the reduction explained by the primes, window test chi^2 with the dual term and jackknife):
    mean and spread of the three bands at 1e12 against the three controls at 1e10 (Welch t; with n = 3 and 3 only indicative).
usage: python3 -I three_bands.py FOLDER_G15 FOLDER_G15bc FOLDER_E15 OUTPUT.json   (outputs of the scripts of study III and of band_windows_jackknife.py)"""
import sys, os, json, math, numpy as np
from scipy import stats
G0, GB, E0, OUT = sys.argv[1:5]
ld = lambda d, n: json.load(open(os.path.join(d, n)))
FAM = ld(G0, 'reference_families_G15.json')['families']; FAM.update(ld(GB, 'reference_families_G15b_G15c.json')['families'])
EW = ld(G0, 'reference_explicit_weights_G15.json')['families']; EW.update(ld(GB, 'reference_explicit_weights_G15b_G15c.json')['families'])
RC = {'G15': ld(G0, 'reference_recomputation_G15.json'), 'G15b': ld(GB, 'reference_recomputation_G15b_G15c.json')}; RC['G15c'] = RC['G15b']
JK = {'G15': ld(G0, 'windows_jackknife_G15.json')}
for f in ('G15b', 'G15c'):
    p = os.path.join(GB, 'windows_jackknife_%s.json' % f); JK[f] = json.load(open(p)) if os.path.exists(p) else None
EWc = ld(E0, 'reference_explicit_weights_E15a_E15b_E15c.json')['families']; RCc = ld(E0, 'reference_recomputation_E15a_E15b_E15c.json')
TR = ld(G0, 'reference_g15_trends.json'); F0 = TR['shift_fit']['without_G15']
PAS, KON = ['G15', 'G15b', 'G15c'], ['E15a', 'E15b', 'E15c']
def band_statistics(f, E, R):
    S = E[f]['stretch']; z = R['z_44'][f]
    return dict(stretch_random=z['stretch_all_starts_ratio'], z840=z['mod840'], share_primes=S['reduction_sum']/S['reduction_exact'], windows_dual=E[f]['windows_summary']['chi2_with_dual'])
o = dict(bands={}, controls={}, summary={})
for f in PAS:
    A, P, M = FAM[f]['all'], FAM[f]['plus'], FAM[f]['minus']; l = P['mean_ln_D'] - math.log(math.pi)
    b = dict(l=l, n=A['n'], n_even=P['n'], n_odd=M['n'], density=A['density_0_05'], density_se=A['density_0_05_se'], density_ratios=A['density_0_05_ratios'], density_usp=A['density_0_05_usp'],
             z_usp=(A['density_0_05'] - A['density_0_05_usp'])/A['density_0_05_se'], z_ratios=(A['density_0_05'] - A['density_0_05_ratios'])/A['density_0_05_se'],
             chi2_ratios=A['chi2_005_ratios'], chi2_ratios_even=P['chi2_005_ratios'], chi2_ratios_odd=M['chi2_005_ratios'], chi2_usp=A['chi2_005_usp'],
             x1=A['x1_mean'], x1_se=A['x1_se'], x1_usp2n_eff=A['x1_usp2n_eff'], N_eff=A['N_eff'],
             shift_even=[P['x1_mean'] - P['x1_usp2n_eff'], P['x1_se']], shift_odd=[M['x1_mean'] - M['x1_usp2n_eff'], M['x1_se']],
             dual=EW[f]['regression']['slope_chi_with_dual'], diagonal=EW[f]['regression']['slope_chi'])
    b.update(band_statistics(f, EW, RC[f]))
    if JK.get(f): b.update(windows_jackknife=JK[f]['chi2_jackknife'], windows_error_ratio=JK[f]['ratio_se_median'], windows_max_z=JK[f]['max_z_jackknife'])
    o['bands'][f] = b
for f in KON: o['controls'][f] = dict(n=EWc[f]['n'], l=EWc[f]['l'], **band_statistics(f, EWc, RCc))
B = o['bands']; w = lambda v, s: (float(np.sum(np.array(v)/np.array(s)**2)/np.sum(1/np.array(s)**2)), float(1/math.sqrt(np.sum(1/np.array(s)**2))))
def combine(v, s):
    m, e = w(v, s); ch = float(np.sum(((np.array(v) - m)/np.array(s))**2)); return dict(mean=m, se=e, chi2=ch, dof=len(v) - 1, p=float(stats.chi2.sf(ch, len(v) - 1)))
S = o['summary']
S['density'] = combine([B[f]['density'] for f in PAS], [B[f]['density_se'] for f in PAS])
H = S['density']; S['density']['usp'] = float(np.mean([B[f]['density_usp'] for f in PAS])); S['density']['ratios'] = float(np.mean([B[f]['density_ratios'] for f in PAS]))
S['density']['z_usp'] = (H['mean'] - H['usp'])/H['se']; S['density']['z_ratios'] = (H['mean'] - H['ratios'])/H['se']
S['x1'] = combine([B[f]['x1'] for f in PAS], [B[f]['x1_se'] for f in PAS]); S['x1']['usp2n_eff'] = float(np.mean([B[f]['x1_usp2n_eff'] for f in PAS])); S['x1']['z_usp2n_eff'] = (S['x1']['mean'] - S['x1']['usp2n_eff'])/S['x1']['se']
S['shift_even'] = combine([B[f]['shift_even'][0] for f in PAS], [B[f]['shift_even'][1] for f in PAS]); S['shift_odd'] = combine([B[f]['shift_odd'][0] for f in PAS], [B[f]['shift_odd'][1] for f in PAS])
lm = float(np.mean([B[f]['l'] for f in PAS])); pred = {k: v['c']*lm**-(1 if k == 'a1' else 2 if k == 'a2' else v['a']) for k, v in F0.items()}
S['shift_even']['l'] = lm; S['shift_even']['predictions_without_G15'] = pred; S['shift_even']['z'] = {k: (S['shift_even']['mean'] - v)/S['shift_even']['se'] for k, v in pred.items()}
S['shift_even']['z_zero'] = S['shift_even']['mean']/S['shift_even']['se']
S['dual'] = combine([B[f]['dual'][0] for f in PAS], [B[f]['dual'][1] for f in PAS]); S['dual']['z_from_1'] = (S['dual']['mean'] - 1)/S['dual']['se']
S['chi2_ratios'] = float(np.mean([B[f]['chi2_ratios'] for f in PAS]))
S['band_statistics'] = {}
for k in ('stretch_random', 'z840', 'share_primes', 'windows_dual'):
    a = np.array([B[f][k] for f in PAS]); c = np.array([o['controls'][f][k] for f in KON]); t = stats.ttest_ind(a, c, equal_var=False)
    S['band_statistics'][k] = dict(bands=a.tolist(), controls=c.tolist(), mean_1e12=float(a.mean()), sd_1e12=float(a.std(ddof=1)), mean_1e10=float(c.mean()), sd_1e10=float(c.std(ddof=1)),
                           welch_t=float(t.statistic), welch_p=float(t.pvalue), overlap=bool(a.max() >= c.min() and c.max() >= a.min()))
json.dump(o, open(OUT, 'w'), indent=1)
print(json.dumps(S, indent=1)[:4000])
