#!/usr/bin/env python3
"""Window test of the explicit formula (study III, Section 4.7) with the error from a grouped jackknife. In a narrow family (one band of width
15,000 in |D|) neighbouring D lie in one band, and the error of a regression coefficient assumes independent L-functions. Here: 50 contiguous groups
by |D|, jackknife leaving out one group at a time.
usage: python3 -I band_windows_jackknife.py STUDY_III_CODE ZEROS_ROOT REFERENCE_EXPLICIT_WEIGHTS.json FAMILY [OUTPUT.json]
  STUDY_III_CODE: the folder with read_zeros.py of study III; ZEROS_ROOT: the store of zeros (ks_nuly.py format, with manifest.json)"""
import sys, os, json, math, numpy as np
sys.path.insert(0, os.path.abspath(sys.argv[1])); import read_zeros as RZ
ROOT, REF, FAM = sys.argv[2:5]; OUT = sys.argv[5] if len(sys.argv) > 5 else None
XW = [(0.1, 0.4), (0.4, 0.8), (0.8, 1.5), (1.5, 2.5), (2.5, 4.0)]; PW = [2, 3, 5, 7, 11, 13]; G = 50
def kron(D, p):
    if p == 2: m = np.mod(D, 8); return np.where(D % 2 == 0, 0, np.where((m == 1) | (m == 7), 1, -1)).astype(float)
    r = np.arange(1, p, dtype=np.int64); tb = -np.ones(p); tb[(r * r) % p] = 1.0; tb[0] = 0.0; return tb[np.mod(D, p)]
Ds, CW = [], []
for D, g in RZ.study_zeros(FAM, ROOT):
    x = np.array([RZ.theta(D, t) / math.pi for t in g]); Ds.append(D); CW.append([np.sum((x > a) & (x <= b)) for a, b in XW])
Ds = np.array(Ds, dtype=np.int64); CW = np.array(CW, float); n = len(Ds); L = np.log(np.abs(Ds) / math.pi); Sg = (Ds < 0).astype(float)
o = np.argsort(np.abs(Ds), kind='stable'); grp = np.empty(n, int); grp[o] = (np.arange(n) * G) // n
R = json.load(open(REF))['families'][FAM]['windows']
def fit(idx, p):
    x = kron(Ds[idx], p); cols = [np.ones(len(idx)), x, x * x] + ([L[idx]] if np.std(L) > 1e-3 else []) + [Sg[idx]]
    F = np.stack(cols, 1); cf, *_ = np.linalg.lstsq(F, CW[idx], rcond=None); rs = CW[idx] - F @ cf
    se = np.sqrt(np.diag(np.linalg.inv(F.T @ F))[1] * (rs * rs).sum(0) / (len(idx) - F.shape[1])); return cf[1], se
allidx = np.arange(n); out = dict(family=FAM, n=n, groups=G, p={}); zn, zj = [], []
for p in PW:
    c0, se0 = fit(allidx, p); jk = np.array([fit(np.nonzero(grp != k)[0], p)[0] for k in range(G)]); sej = np.sqrt((G - 1) / G * ((jk - jk.mean(0))**2).sum(0))
    pr = np.array(R[str(p)]['predicted_chi_with_dual']); zn += list((c0 - pr) / se0); zj += list((c0 - pr) / sej)
    out['p'][str(p)] = dict(fit=c0.tolist(), se_regression=se0.tolist(), se_jackknife=sej.tolist(), ratio=(sej / se0).tolist(), prediction=pr.tolist())
zn, zj = np.array(zn), np.array(zj)
out.update(chi2_regression=float((zn**2).sum()), chi2_jackknife=float((zj**2).sum()), coefficients=len(zn), ratio_se_median=float(np.median([r for v in out['p'].values() for r in v['ratio']])), max_z_jackknife=float(np.abs(zj).max()))
print(FAM, 'n', n, 'chi2 (30 coef.) regression %.1f, jackknife %.1f; median se_jk/se_reg %.2f; max |z| jk %.1f' % (out['chi2_regression'], out['chi2_jackknife'], out['ratio_se_median'], out['max_z_jackknife']))
if OUT: json.dump(out, open(OUT, 'w'), indent=1)
