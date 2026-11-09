#!/usr/bin/env python3
"""Figures 7 and 8 of study IV (Section 5; light background, English and Czech, PNG 300 dpi and PDF) from reference_quadratic.json.
usage: python3 -I figures_quadratic.py reference_quadratic.json OUTPUT_FOLDER"""
import sys, json, math, os, numpy as np
import matplotlib; matplotlib.use('Agg'); import matplotlib.pyplot as plt
from matplotlib.ticker import FuncFormatter, ScalarFormatter
R = json.load(open(sys.argv[1])); OUT = sys.argv[2]; os.makedirs(OUT, exist_ok=True)
plt.rcParams.update({'font.family': 'DejaVu Sans', 'font.size': 9, 'axes.titlesize': 10, 'axes.labelsize': 9, 'legend.fontsize': 8, 'axes.edgecolor': '#555555',
                     'axes.linewidth': 0.8, 'xtick.color': '#333333', 'ytick.color': '#333333', 'axes.grid': True, 'grid.color': '#e3e3e3', 'grid.linewidth': 0.6,
                     'figure.facecolor': 'white', 'axes.facecolor': 'white', 'savefig.facecolor': 'white', 'legend.frameon': False})
INK, MUT, RAT, USP, BND, CTL, POW, CL = '#222222', '#8a8a8a', '#1f5fa8', '#d07a1a', '#c0392b', '#7b5ea7', '#b2589a', '#2e8b57'
BC = {'G15': '#c0392b', 'G15b': '#e67e22', 'G15c': '#17a589'}
CZF = FuncFormatter(lambda v, p: ('%g' % round(v, 6)).replace('.', ',').replace('-', '−'))
def comma(fig):
    for ax in fig.axes:
        if ax.get_yscale() == 'linear' and isinstance(ax.yaxis.get_major_formatter(), ScalarFormatter): ax.yaxis.set_major_formatter(CZF)
def save(fig, name):
    if name.endswith('_cz'): comma(fig)
    fig.savefig(os.path.join(OUT, name + '.png'), dpi=300, bbox_inches='tight'); fig.savefig(os.path.join(OUT, name + '.pdf'), bbox_inches='tight', metadata={'CreationDate': None}); plt.close(fig)
T = R['table']; Cb = R['combined_1e12']; F7 = R['fits_7']; U = R['usp']; K = np.array(R['ratios_curve']); TH = R['thresholds']; TB = R['three_bands']
lg = lambda l: (l + math.log(math.pi))/math.log(10)
FAM = ['B10', 'C100', 'D100', 'E100']; BANDS = ['G15', 'G15b', 'G15c']; CON = ['E15a', 'E15b', 'E15c']
for en in (False, True):
    sfx = '_en' if en else '_cz'; w = lambda c, e: e if en else c; f = lambda x, d: ('%.*f' % (d, x)) if en else ('%.*f' % (d, x)).replace('.', ',').replace('-', '−')
    # 7) slow convergence
    fig, axs = plt.subplots(1, 2, figsize=(10.0, 3.8))
    ax = axs[0]; Kd = lg(K[:, 0]); m = (Kd >= 3) & (Kd <= 13.2)
    ax.plot(Kd[m], (K[m, 0]*(K[m, 1] - U)), color=RAT, lw=1.4, label=w('domněnka poměrů', 'ratios conjecture'))
    ax.axhline(F7['c_over_l']['c'], color=CL, lw=1, ls='--', label=w('fit c/l: c = %s' % f(F7['c_over_l']['c'], 2), 'fit c/l: c = %s' % f(F7['c_over_l']['c'], 2)))
    for k in FAM: ax.errorbar(lg(T[k]['l']), T[k]['l_excess'], yerr=T[k]['l_excess_se'], fmt='o', ms=4, color=INK, lw=0.9)
    for j, k in enumerate(BANDS): ax.errorbar(lg(T[k]['l']) + 0.12*(j - 1), T[k]['l_excess'], yerr=T[k]['l_excess_se'], fmt='o', ms=4, color=BC[k], lw=0.9, label=k)
    ax.errorbar(lg(Cb['l']) + 0.45, Cb['l_excess'], yerr=Cb['l_excess_se'], fmt='D', ms=5, color=INK, mfc='white', lw=1, label=w('tři pásma spolu', 'three bands combined'))
    for k in FAM: ax.annotate(k, (lg(T[k]['l']), T[k]['l_excess'] + T[k]['l_excess_se'] + 0.03), ha='center', fontsize=7, color=MUT)
    ax.set_xlabel('$\\log_{10}|D|$'); ax.set_ylabel(w('(hustota − USp) × l', '(density − USp) × l')); ax.set_ylim(1.9, 3.3); ax.legend(loc='lower left', fontsize=7, ncol=2)
    ax.set_title(w('Přebytek nad USp krát l = ln(|D|/π): skoro konstantní', 'Excess over USp times $l = \\ln(|D|/\\pi)$: nearly constant'), fontsize=9)
    ax = axs[1]; d = np.linspace(2.5, 120, 600); l = d*math.log(10) - math.log(math.pi); ok = l > 3
    ax.plot(Kd, K[:, 1] - U, color=RAT, lw=1.6, label=w('domněnka poměrů (bez volného parametru)', 'ratios conjecture (no free parameter)'))
    ax.plot(d[ok], F7['c_over_l']['c']/l[ok], color=CL, lw=1.1, ls='--', label='c/l')
    ax.plot(d, np.exp(F7['power']['log_amplitude'] + F7['power']['exponent']*d*math.log(10)), color=POW, lw=1.1, ls=':', label=w('mocnina |D|, exponent %s' % f(F7['power']['exponent'], 3), 'power of |D|, exponent %s' % f(F7['power']['exponent'], 3)))
    for k in FAM + BANDS: ax.plot(lg(T[k]['l']), T[k]['density'] - T[k]['usp'], 'o', ms=3.5, color=BC.get(k, INK))
    for c in ('0.05', '0.01'):
        ax.axhline(float(c), color=MUT, lw=0.7, ls='-.'); x0 = TH[c]['log10D']; ax.plot(x0, float(c), 'o', ms=6, mfc='none', color=RAT)
        e = math.floor(x0); ax.annotate('$%s\\cdot10^{%d}$' % (f(10**(x0 - e), 1), e), (x0, float(c)), xytext=(4, 5), textcoords='offset points', fontsize=7.5, color=RAT)
    ax.set_xlim(0, 120); ax.set_ylim(0, 0.32); ax.set_xlabel('$\\log_{10}|D|$'); ax.set_ylabel(w('přebytek hustoty v [0, ½] nad USp', 'excess of the density in $[0, 1/2]$ over USp')); ax.legend(loc='upper right', fontsize=7)
    ax.set_title(w('Za daty se popisy rozcházejí o řády', 'Beyond the data the descriptions part by orders of magnitude'), fontsize=9)
    fig.tight_layout(); save(fig, '07_slow_convergence' + sfx)
    # 8) three bands at 1e12 against the controls at 1e10
    fig, axs = plt.subplots(2, 3, figsize=(10.4, 5.6)); S = TB['summary']; Bd = TB['bands']; Ct = TB['controls']
    ax = axs[0, 0]; x = np.arange(4)
    for j, k in enumerate(BANDS): ax.errorbar(j, Bd[k]['density'], yerr=Bd[k]['density_se'], fmt='o', ms=5, color=BC[k], lw=1)
    ax.errorbar(3, S['density']['mean'], yerr=S['density']['se'], fmt='D', ms=5, color=INK, mfc='white', lw=1)
    ax.axhline(S['density']['ratios'], color=RAT, lw=1.2, label=w('poměry', 'ratios')); ax.axhline(S['density']['usp'], color=USP, lw=1.2, ls='--', label='USp')
    ax.set_xticks(x); ax.set_xticklabels(BANDS + [w('spolu', 'comb.')]); ax.set_ylabel(w('hustota v [0, ½]', 'density in $[0, 1/2]$')); ax.legend(loc='lower left', fontsize=7)
    ax.set_title(w('%s chyby nad USp (spolu)' % f(S['density']['z_usp'], 1), '%s errors above USp (combined)' % f(S['density']['z_usp'], 1)), fontsize=9)
    ax = axs[0, 1]
    for j, k in enumerate(BANDS): ax.errorbar(j, Bd[k]['x1'], yerr=Bd[k]['x1_se'], fmt='o', ms=5, color=BC[k], lw=1)
    ax.errorbar(3, S['x1']['mean'], yerr=S['x1']['se'], fmt='D', ms=5, color=INK, mfc='white', lw=1); ax.axhline(S['x1']['usp2n_eff'], color=USP, lw=1.2, label='USp$(2N_{\\rm eff})$')
    ax.set_xticks(x); ax.set_xticklabels(BANDS + [w('spolu', 'comb.')]); ax.set_ylabel(w('průměrná první nula', 'mean first zero')); ax.legend(loc='lower left', fontsize=7)
    ax.set_title(w('První nula: %s chyby od USp(2N_ef)' % f(S['x1']['z_usp2n_eff'], 1), 'First zero: %s errors from USp$(2N_{\\rm eff})$' % f(S['x1']['z_usp2n_eff'], 1)), fontsize=9)
    stat = [('stretch_random', w('úsek / náhodné', 'stretch / random')), ('z840', w('mod 840: z proti náhodným', 'mod 840: z against random')), ('share_primes', w('podíl vysvětlený prvočísly', 'share explained by the primes')), ('windows_dual', w('test oken χ² (30 koef.)', 'window test $\\chi^2$ (30 coef.)'))]
    for (key, tit), ax in zip(stat, [axs[0, 2], axs[1, 0], axs[1, 1], axs[1, 2]]):
        P = S['band_statistics'][key]; v = P['bands'] + P['controls']; cols = [BC[k] for k in BANDS] + [CTL]*3
        ax.bar(np.arange(6), v, color=cols, width=0.65); ax.axhline(0, color=MUT, lw=0.6)
        ax.plot([-0.4, 2.4], [P['mean_1e12']]*2, color=INK, lw=0.9, ls=':'); ax.plot([2.6, 5.4], [P['mean_1e10']]*2, color=CTL, lw=0.9, ls=':')
        ax.set_xticks(np.arange(6)); ax.set_xticklabels(BANDS + CON, rotation=40, fontsize=7)
        ax.set_title(tit + ' · t = %s, p = %s' % (f(P['welch_t'], 2), f(P['welch_p'], 2)), fontsize=8.5)
    fig.suptitle(w('Tři pásma u |D| = 10¹² (G15, G15b, G15c) proti třem kontrolním úsekům u 10¹⁰ (E15a, E15b, E15c)', 'Three bands at $|D| = 10^{12}$ (G15, G15b, G15c) against three control stretches at $10^{10}$ (E15a, E15b, E15c)'), fontsize=10)
    fig.tight_layout(); save(fig, '08_three_bands' + sfx)
print('ok', sorted(os.listdir(OUT))[-8:])
