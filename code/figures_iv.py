#!/usr/bin/env python3
"""Figures of study IV (light background, English and Czech), PNG 300 dpi and PDF, from reference_study_iv.json.
usage: python3 -I figures_iv.py DATA_DIR OUTPUT_FOLDER  (DATA_DIR with reference_study_iv.json and arrays_study_iv.npz)"""
import sys, json, math, os, numpy as np
import matplotlib; matplotlib.use('Agg'); import matplotlib.pyplot as plt
from matplotlib.ticker import FuncFormatter, ScalarFormatter
from scipy import stats
R = json.load(open(os.path.join(sys.argv[1], 'reference_study_iv.json'))); AR = np.load(os.path.join(sys.argv[1], 'arrays_study_iv.npz')); OUT = sys.argv[2]; os.makedirs(OUT, exist_ok=True)
plt.rcParams.update({'font.family': 'DejaVu Sans', 'font.size': 9, 'axes.titlesize': 10, 'axes.labelsize': 9, 'legend.fontsize': 8, 'axes.edgecolor': '#555555',
                     'axes.linewidth': 0.8, 'xtick.color': '#333333', 'ytick.color': '#333333', 'axes.grid': True, 'grid.color': '#e3e3e3', 'grid.linewidth': 0.6,
                     'figure.facecolor': 'white', 'axes.facecolor': 'white', 'savefig.facecolor': 'white', 'legend.frameon': False})
INK, MUT, RAT, GUE, NOI, CS, ODL, PAI, ZET, CUE = '#222222', '#8a8a8a', '#1f5fa8', '#d07a1a', '#7fb3e6', '#c0392b', '#7b5ea7', '#2e8b57', '#b2589a', '#6b6b6b'
def cz(s): return s.replace('.', ',')
CZF = FuncFormatter(lambda v, p: ('%g' % round(v, 6)).replace('.', ',').replace('-', '−'))
def comma(fig):
    for ax in fig.axes:
        if ax.get_yscale() == 'linear' and isinstance(ax.yaxis.get_major_formatter(), ScalarFormatter): ax.yaxis.set_major_formatter(CZF)
        if ax.get_xscale() == 'linear' and isinstance(ax.xaxis.get_major_formatter(), ScalarFormatter): ax.xaxis.set_major_formatter(CZF)
def save(fig, name):
    if name.endswith('_cz'): comma(fig)
    fig.savefig(os.path.join(OUT, name + '.png'), dpi=300, bbox_inches='tight'); fig.savefig(os.path.join(OUT, name + '.pdf'), bbox_inches='tight', metadata={'CreationDate': None}); plt.close(fig)
B = R['blocks']; K = list(B); X = np.arange(len(K)); EXP = R['expectation_chi2']; VT = R['var_t2']
S = R['blocks_summary']; HB = B['H3e10']; NERR = round((HB['g3_55'] - HB['r3_55'][0]) / HB['r3_55'][1]); PCT = round(100 * (1 - HB['r3_55'][0] / HB['g3_55']))
LMID = np.array([B[k]['L_med'] for k in K])   # median height of the zeros of the block (in the lowest block most zeros lie near its top)
for en in (False, True):
    sfx = '_en' if en else '_cz'; w = lambda c, e: e if en else c; f = lambda x, d: ('%.*f' % (d, x)) if en else cz('%.*f' % (d, x))
    # 1) Conrey–Snaith 2008 repeated: size of the differences against pure noise
    C = R['cs2008']; rep, our, noi = C['reported'], C['models']['heights_averaged'], C['noise']
    fig, ax = plt.subplots(figsize=(6.4, 3.3)); lab = [w('směrodatná odchylka', 'standard deviation'), w('střední |rozdíl|', 'mean |difference|'), w('maximum', 'maximum')]
    vals = [[rep['sd'], rep['mean_abs'], rep['max']], [our['sd'], our['mean_abs'], our['max']], [noi['sd'][0], noi['mean_abs'][0], noi['max'][0]]]
    col = [CS, INK, NOI]; nm = [w('Conrey a Snaith 2008 (uvedeno)', 'Conrey and Snaith 2008 (reported)'), w('naše data − vzorec (průměr přes výšky)', 'our data − formula (averaged over heights)'), w('čistý šum (model + chyba koše, 2 000 simulací)', 'pure noise (model + error of the bin, 2,000 simulations)')]
    for j in range(3): ax.bar(np.arange(3) + (j - 1)*0.26, vals[j], width=0.24, color=col[j], label=nm[j])
    ax.errorbar([2 + 0.26], [noi['max'][0]], yerr=[[noi['max'][0] - noi['max'][1]], [noi['max'][2] - noi['max'][0]]], fmt='none', ecolor=INK, lw=1, capsize=3)
    ax.axhline(rep['expected_formula_error'], color=RAT, lw=1, ls='--', label=w('čekaná chyba vzorce $T^{-1/2}$', 'expected error of the formula $T^{-1/2}$'))
    ax.set_xticks(np.arange(3)); ax.set_xticklabels(lab); ax.set_ylabel(w('rozdíl data − vzorec', 'difference data − formula')); ax.legend(loc='upper left', fontsize=7.5)
    ax.set_title(w('Prvních 100 000 nul: rozdíly mají velikost šumu (χ² na koš %s, čeká se %s)' % (f(our['chi2_ab02'], 3), f(EXP, 3)), 'First 100,000 zeros: the differences have the size of the noise ($\\chi^2$ per bin %s, expected %s)' % (f(our['chi2_ab02'], 3), f(EXP, 3))), fontsize=9)
    save(fig, '01_conrey_snaith_2008' + sfx)
    # 2) chi^2 per bin in 14 blocks: formula (old and new) and null models
    fig, axs = plt.subplots(1, 2, figsize=(9.6, 3.6))
    ax = axs[0]; q = 3*math.sqrt(VT/784); ax.axhspan(EXP - q, EXP + q, color=RAT, alpha=0.08, lw=0); ax.axhline(EXP, color=RAT, lw=1, ls='--')
    ax.plot(X - 0.12, [B[k]['chi3_old'] for k in K], 'o', ms=3.5, color=MUT, label=w('chyba = rozptyl kusů', 'error = scatter of the pieces')); ax.plot(X + 0.12, [B[k]['chi3'] for k in K], 'o', ms=4, color=INK, label=w('chyba proti modelu ve výšce kusu, a, b ≥ 0,2', 'error against the model at the height of the piece, a, b ≥ 0.2'))
    ax.set_xticks(X); ax.set_xticklabels(K, rotation=60, fontsize=7); ax.set_ylim(0.75, 1.35); ax.set_ylabel(w('χ² na koš R₃ proti vzorci', '$\\chi^2$ per bin of $R_3$ against the formula')); ax.legend(loc='upper left', fontsize=7)
    ax.set_title(w('Vzorec Conreyho a Snaithové: očekávání %s ± 3σ' % f(EXP, 3), 'Formula of Conrey and Snaith: expectation %s ± 3σ' % f(EXP, 3)), fontsize=9)
    ax = axs[1]
    for key, c, l in (('chi3_cue_neff', CUE, w('CUE(N_ef), 1. řád', 'CUE($N_{\\rm eff}$), first order')), ('chi3_gue', GUE, 'GUE'), ('chi3_zeta_only', ZET, w('vzorec bez aritmetiky', 'formula without arithmetic')), ('chi3_cue_fit', '#999999', w('CUE s volným N', 'CUE with free $N$')),
                      ('chi3_pairs_only', PAI, w('jen páry', 'pairs only')), ('chi3', INK, w('vzorec', 'formula'))):
        ax.plot(X, [B[k][key] for k in K], 'o-', ms=3, lw=0.9, color=c, label=l)
    ax.set_yscale('log'); ax.axhline(EXP, color=RAT, lw=0.8, ls='--'); ax.set_xticks(X); ax.set_xticklabels(K, rotation=60, fontsize=7); ax.set_ylabel(w('χ² na koš R₃', '$\\chi^2$ per bin of $R_3$')); ax.legend(loc='upper left', fontsize=6.8, ncol=2)
    ax.set_title(w('Nulové modely a části vzorce (stejná chyba a maska)', 'Null models and parts of the formula (same errors and mask)'), fontsize=9)
    fig.tight_layout(); save(fig, '02_chi2_blocks' + sfx)
    # 3) R3 in the bin (0.5-0.6)^2 over the height
    fig, ax = plt.subplots(figsize=(6.2, 3.4)); L = LMID
    ax.errorbar(L, [B[k]['r3_55'][0] for k in K], yerr=[B[k]['r3_55'][1] for k in K], fmt='o', ms=3.5, color=INK, lw=0.8, label=w('nuly (Platt, LMFDB)', 'zeros (Platt, LMFDB)'))
    ax.plot(L, [B[k]['m3_55'] for k in K], 's', ms=4, mfc='none', color=RAT, label=w('vzorec (domněnka poměrů)', 'formula (ratios conjecture)')); ax.axhline(B[K[0]]['g3_55'], color=GUE, lw=1.2, label='GUE')
    ax.set_xlabel(w('výška L = ln(t/2π) (medián bloku)', 'height $L = \\ln(t/2\\pi)$ (median of the block)')); ax.set_ylabel('$R_3$ ' + w('v koši', 'in the bin') + ' $(0.5\\text{–}0.6)^2$')
    ax.legend(loc='lower right', fontsize=7.5)
    ax.set_title(w('Tři blízké nuly: u t = 3·10¹⁰ stále asi %d chyb (%d %%) pod GUE' % (NERR, PCT), 'Three close zeros: at $t = 3\\cdot10^{10}$ still about %d errors (%d %%) below GUE' % (NERR, PCT)), fontsize=9)
    fig.tight_layout(); save(fig, '03_r3_close_triples' + sfx)
    # 4) amplitude lambda over the height, with Odlyzko's sets
    fig, ax = plt.subplots(figsize=(7.0, 3.5)); O = R['odlyzko']; Oc = R['odlyzko_combined']
    ax.axhline(1, color=RAT, lw=1, ls='--', label=w('domněnka poměrů (λ = 1)', 'ratios conjecture ($\\lambda = 1$)')); ax.axhline(0, color=GUE, lw=1.2, label='GUE ($\\lambda = 0$)')
    ax.errorbar(LMID, [B[k]['lam3']['lam'] for k in K], yerr=[B[k]['lam3']['err'] for k in K], fmt='o', ms=3.5, color=INK, lw=0.8, label=w('Platt, 14 bloků', 'Platt, 14 blocks'))
    LO = [math.log(O[k]['t']/(2*math.pi)) for k in O]; ax.errorbar(LO, [O[k]['lam3'][0] for k in O], yerr=[O[k]['lam3'][1] for k in O], fmt='D', ms=4, color=ODL, lw=0.8, label=w('Odlyzko, po 10⁴ nulách', 'Odlyzko, $10^4$ zeros each'))
    ax.errorbar([max(LO) + 1.2], [Oc['lam3'][0]], yerr=[Oc['lam3'][1]], fmt='*', ms=8, color=ODL, lw=0.8, label=w('Odlyzko spolu', 'Odlyzko together'))
    ax.set_xlabel(w('výška L = ln(t/2π)', 'height $L = \\ln(t/2\\pi)$')); ax.set_ylabel(w('amplituda λ₃', 'amplitude $\\lambda_3$')); ax.set_ylim(-2.3, 2.6); ax.legend(loc='upper center', fontsize=7, ncol=3)
    ins = ax.inset_axes([0.27, 0.08, 0.42, 0.33]); ins.axhline(1, color=RAT, lw=0.8, ls='--')
    ins.errorbar(LMID, [B[k]['lam3']['lam'] for k in K], yerr=[B[k]['lam3']['err'] for k in K], fmt='o', ms=2.2, color=INK, lw=0.6); ins.set_ylim(0.975, 1.025); ins.tick_params(labelsize=6.5); (ins.xaxis.set_major_formatter(CZF), ins.yaxis.set_major_formatter(CZF)) if not en else None; ins.set_title(w('zvětšeno: %s–%s' % (f(S['lam3'][0], 3), f(S['lam3'][1], 3)), 'zoom: %s–%s' % (f(S['lam3'][0], 3), f(S['lam3'][1], 3))), fontsize=7)
    ax.set_title(w('Odchylka od GUE jedním číslem: λ₃ = 1 na %s–%s %% do 3·10¹⁰' % (f(100 * S['lam3_err'][0], 1), f(100 * S['lam3_err'][1], 1)), 'The deviation from GUE as one number: $\\lambda_3 = 1$ within %s–%s %% up to $3\\cdot10^{10}$' % (f(100 * S['lam3_err'][0], 1), f(100 * S['lam3_err'][1], 1))), fontsize=9)
    fig.tight_layout(); save(fig, '04_amplitude' + sfx)
    # 5) height inside a block: pieces against the prediction
    fig, axs = plt.subplots(1, 3, figsize=(9.6, 3.0), sharey=False)
    for ax, k in zip(axs, ('L1e5', 'Z6a', 'S1')):
        b = B[k]; xk = (np.arange(40) - 19.5)/40; l3 = AR[k + '_l3k']; p3 = AR[k + '_p3k']; la = b['lam3']['lam']
        ax.plot(xk, l3, 'o', ms=3, color=INK, label=w('kusy (40)', 'pieces (40)')); ax.plot(xk, la*p3, '-', color=RAT, lw=1.4, label=w('předpověď ve výšce kusu', 'prediction at the height of the piece')); ax.axhline(np.mean(l3), color=MUT, lw=0.8, ls=':', label=w('průměr kusů', 'mean of the pieces'))
        ax.set_title('%s: %s %s ± %s, %s %s' % (k, w('sklon', 'slope'), f(b['lam3']['slope_obs'], 3), f(b['lam3']['slope_err'], 3), w('předpověď', 'predicted'), f(b['lam3']['slope_pred'], 3)), fontsize=8.5); ax.set_xlabel(w('poloha kusu v bloku', 'position of the piece in the block'))
    axs[0].set_ylabel(w('λ₃ kusu', '$\\lambda_3$ of the piece')); axs[0].legend(loc='upper right', fontsize=6.8); fig.tight_layout(); save(fig, '05_height_inside_block' + sfx)
    # 6) bin deviations: map for the highest block and the distribution of all z against t39
    fig, axs = plt.subplots(1, 2, figsize=(8.6, 3.6)); z = np.array(R['z_maps']['H3e10'], float).reshape(30, 30); zz = np.where(np.arange(30)[:, None] >= 2, 1, np.nan) * np.where(np.arange(30)[None, :] >= 2, 1, np.nan) * z
    im = axs[0].imshow(zz.T, origin='lower', extent=[0, 3, 0, 3], cmap='RdBu_r', vmin=-4, vmax=4); fig.colorbar(im, ax=axs[0], fraction=0.046, pad=0.04, label='(data − ' + w('vzorec', 'formula') + ')/' + w('chyba', 'error'))
    axs[0].set_xlabel('$a$'); axs[0].set_ylabel('$b$'); axs[0].set_title(w('H3e10 (t = 3,06·10¹⁰): odchylky košů', 'H3e10 ($t = 3.06\\cdot10^{10}$): deviations of bins'), fontsize=9); axs[0].grid(False)
    zl = z[2:, 2:].ravel(); zl = zl[np.isfinite(zl)]; z1 = np.array(R['z_maps']['L1e5'], float).reshape(30, 30)[2:, 2:].ravel(); z1 = z1[np.isfinite(z1)]; zall = np.concatenate([zl, z1])
    hb = np.linspace(-5, 5, 41); axs[1].hist(zall, bins=hb, density=True, color=NOI, label=w('L1e5 a H3e10 (1 568 košů)', 'L1e5 and H3e10 (1,568 bins)')); xx = np.linspace(-5, 5, 400); axs[1].plot(xx, stats.t.pdf(xx, 39), color=INK, lw=1.2, label=w('t-rozdělení, 39 st. v.', '$t$ distribution, 39 d.o.f.'))
    axs[1].set_yscale('log'); axs[1].set_ylim(1e-4, 1); axs[1].set_xlabel('$z$'); axs[1].legend(loc='upper right', fontsize=7); axs[1].set_title(w('Rozdělení odchylek košů', 'Distribution of the deviations of bins'), fontsize=9)
    fig.tight_layout(); save(fig, '06_bin_deviations' + sfx)
print('ok', sorted(os.listdir(OUT))[:20])
