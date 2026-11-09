"""CSV exports of study IV from reference_study_iv.json and arrays_study_iv.npz.

Usage:  python3 -I export_csv.py DATA_DIR
writes DATA_DIR/blocks.csv, bins_r3.csv, bins_r2.csv, pieces.csv, odlyzko.csv,
cs2008_bins.csv and grid_models.csv.
"""
import os
import sys
import csv
import json
import numpy as np

D = sys.argv[1]
R = json.load(open(os.path.join(D, 'reference_study_iv.json')))
A = dict(np.load(os.path.join(D, 'arrays_study_iv.npz')))
g = lambda x: '%.10g' % x
gx, gw = np.polynomial.legendre.leggauss(2)
u = (np.arange(30) * 0.1)[:, None] + 0.05 * (gx[None, :] + 1)
GUE3 = ((1 - np.sinc(u[:, None, :, None])**2 - np.sinc(u[None, :, None, :])**2 - np.sinc(u[:, None, :, None] + u[None, :, None, :])**2
         + 2 * np.sinc(u[:, None, :, None]) * np.sinc(u[None, :, None, :]) * np.sinc(u[:, None, :, None] + u[None, :, None, :]))
        * (gw / 2)[:, None] * (gw / 2)[None, :]).sum((2, 3)).ravel()
x4, w4 = np.polynomial.legendre.leggauss(4)
GUE2 = ((1 - np.sinc((np.arange(60) * 0.05)[:, None] + 0.025 * (x4[None, :] + 1))**2) * w4 / 2).sum(1)


def write(name, head, rows):
    with open(os.path.join(D, name), 'w', newline='') as fh:
        w = csv.writer(fh)
        w.writerow(head)
        w.writerows(rows)
    print(name, len(rows))


B = R['blocks']
cols = ['n', 't_first', 't_last', 't_med', 'Neff', 'chi3', 'chi3_old', 'chi3_gue', 'chi3_zeta_only', 'chi3_pairs_only', 'chi3_cue_neff', 'chi3_cue_fit', 'N_fit_r3', 'N_fit_r2',
        'chi2', 'chi2_old', 'chi2_gue', 'error_ratio_rms', 'z3_max', 'z3_over3', 'r2_ratio_055']
write('blocks.csv', ['block'] + cols + ['L_min', 'L_max', 'r3_55', 'r3_55_err', 'm3_55', 'gue3_55', 'lambda3', 'lambda3_err', 'lambda2', 'lambda2_err',
                                         'slope_pieces_obs', 'slope_pieces_err', 'slope_pieces_pred', 'lambda3_err_scatter_of_pieces'],
      [[k] + [g(b[c]) if isinstance(b[c], float) else b[c] for c in cols] + [g(b['L'][0]), g(b['L'][1]), g(b['r3_55'][0]), g(b['r3_55'][1]), g(b['m3_55']), g(b['g3_55']),
                                                                             g(b['lam3']['lam']), g(b['lam3']['err']), g(b['lam2']['lam']), g(b['lam2']['err']),
                                                                             g(b['lam3']['slope_obs']), g(b['lam3']['slope_err']), g(b['lam3']['slope_pred']), g(b['lam3']['err_pieces_scatter'])] for k, b in B.items()])
rows = []
for k in B:
    for i in range(900):
        a, b = divmod(i, 30)
        rows.append([k, '%.1f' % (a * 0.1), '%.1f' % (b * 0.1), int(a >= 2 and b >= 2), g(A[k + '_y3'][i]), g(A[k + '_e3'][i]), g(A[k + '_e3_pieces'][i]), g(A[k + '_m3'][i]),
                     g(A[k + '_m3_q5'][i]), g(A[k + '_m3_zeta_q5'][i]), g(A[k + '_m3_pairs_q5'][i]), g(GUE3[i])])
write('bins_r3.csv', ['block', 'a_from', 'b_from', 'in_chi2_mask', 'R3', 'error', 'error_scatter_of_pieces', 'formula_height_of_pieces', 'formula_5_quantiles', 'formula_without_arithmetic',
                      'pairs_only', 'GUE'], rows)
rows = [[k, '%.2f' % (c * 0.05), g(A[k + '_y2'][c]), g(A[k + '_e2'][c]), g(A[k + '_e2_pieces'][c]), g(A[k + '_m2'][c]), g(A[k + '_m2_q10'][c]), g(GUE2[c])] for k in B for c in range(60)]
write('bins_r2.csv', ['block', 'd_from', 'R2', 'error', 'error_scatter_of_pieces', 'formula_height_of_pieces', 'formula_10_quantiles', 'GUE'], rows)
rows = [[k, j + 1, g(A[k + '_l3k'][j]), g(A[k + '_p3k'][j]), g(A[k + '_l2k'][j]), g(A[k + '_p2k'][j])] for k in B for j in range(40)]
write('pieces.csv', ['block', 'piece', 'lambda3', 'prediction3', 'lambda2', 'prediction2'], rows)
O = R['odlyzko']
write('odlyzko.csv', ['set', 'first_zero_number', 'zeros', 'median_height', 'L', 'lambda3', 'lambda3_err', 'lambda2', 'lambda2_err'],
      [[k, v['first_zero_number'], v['zeros'], g(v['t']), g(v['L']), g(v['lam3'][0]), g(v['lam3'][1]), g(v['lam2'][0]), g(v['lam2'][1])] for k, v in O.items()])
rows = []
for i in range(900):
    a, b = divmod(i, 30)
    rows.append(['%.1f' % (a * 0.1), '%.1f' % (b * 0.1), int(a >= 1 and b >= 1), g(A['cs2008_y'][i]), g(A['cs2008_e'][i]), g(A['cs2008_e_raw'][i]), g(A['cs2008_m_avg'][i]), g(A['cs2008_m_T'][i])])
write('cs2008_bins.csv', ['a_from', 'b_from', 'in_comparison', 'R3', 'error', 'error_scatter_of_pieces', 'formula_averaged_over_heights', 'formula_at_T'], rows)
rows = [[g(A['grid_L'][j]), 'R3', '%.1f' % (divmod(i, 30)[0] * 0.1), '%.1f' % (divmod(i, 30)[1] * 0.1), g(A['grid_r3'][j, i])] for j in range(len(A['grid_L'])) for i in range(900)]
rows += [[g(A['grid_L'][j]), 'R2', '%.2f' % (c * 0.05), '', g(A['grid_r2'][j, c])] for j in range(len(A['grid_L'])) for c in range(60)]
write('grid_models.csv', ['L', 'function', 'a_or_d_from', 'b_from', 'formula'], rows)
