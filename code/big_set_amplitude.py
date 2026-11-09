"""Amplitude and chi^2 for a large set of zeros at one height (output of stream_histograms.py).

The height range of such a set is negligible (1e9 zeros near t = 1e22 span
ln(t/2 pi) by about 1e-14), so the formula is evaluated at the middle height;
the error of a bin is the scatter of the K pieces about their mean / sqrt(K).
Same estimators as compute_study_iv.py (step odlyzko): weights (m - g)/m,
R3 bins with a, b >= 0.1 in the amplitude, chi^2 of R3 on a, b >= 0.2.

Usage:  python3 -I big_set_amplitude.py HIST.npz TABLES.npz OUT.json
"""
import os
import sys
import json
import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from correlations import gue_r3, gue_r2, chi2_per_bin, MASK_01, MASK_02
from ratios_formula import Tables, Primes, r3_bins, r2_bins

H = np.load(sys.argv[1])
tab, pr = Tables(sys.argv[2]), Primes(30_000)
t = float(H['t_median_estimate'])
m3, m2 = r3_bins(tab, pr, t), r2_bins(tab, t)
g3, g2 = gue_r3(2), gue_r2(4)
r3k, r2k = H['r3k'], H['r2k']
K = r3k.shape[0]
M1 = MASK_01.reshape(30, 30)
d3 = m3 - g3
w3 = np.where(M1 & (m3 > 0), d3 / np.maximum(m3, 1e-9), 0)
d2 = m2 - g2
w2 = np.where(m2 > 0, d2 / np.maximum(m2, 1e-9), 0)
l3 = ((r3k - g3) * w3).sum((1, 2)) / (w3 * d3).sum()
l2 = ((r2k - g2) * w2).sum(1) / (w2 * d2).sum()
y3, e3 = r3k.mean(0).ravel(), r3k.std(0, ddof=1).ravel() / np.sqrt(K)
y2, e2 = r2k.mean(0), r2k.std(0, ddof=1) / np.sqrt(K)
me = lambda l: [float(l.mean()), float(l.std(ddof=1) / np.sqrt(K))]
out = dict(zeros=int(H['zeros']), left_ends=int(H['n']), t=t, L=float(np.log(t / (2 * np.pi))), overflow=int(H['overflow']),
           lam3=me(l3), lam2=me(l2), chi3_formula=chi2_per_bin(y3, m3.ravel(), e3, MASK_02), chi3_gue=chi2_per_bin(y3, g3.ravel(), e3, MASK_02),
           chi2_formula=chi2_per_bin(y2, m2, e2, np.ones(60, bool)), chi2_gue=chi2_per_bin(y2, g2, e2, np.ones(60, bool)),
           r3_55=[float(y3[155]), float(e3[155])], formula_55=float(m3.ravel()[155]), gue_55=float(g3.ravel()[155]), expectation_chi2=(K - 1) / (K - 3))
json.dump(out, open(sys.argv[3], 'w'), indent=1)
print(json.dumps(out, indent=1))
