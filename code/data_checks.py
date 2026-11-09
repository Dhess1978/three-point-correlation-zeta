"""Checks of the zeros of zeta used in study IV (Section 2).

1. md5 of every Platt file against md5.txt of the LMFDB;
2. consecutive files continue exactly (N1 of a file = N0 of the next) and the
   spacing across the boundary is an ordinary spacing;
3. the zero counting function at the first and last zero of every file:
   N(t) = theta(t)/pi + 1 + S(t) with N = n - 1/2 at the n-th zero, |S| < 1;
4. against Odlyzko's table of the first 2,001,052 zeros (zeros6, 9 decimals);
5. against the first 1e5 zeros of the LMFDB in text (30 digits);
6. against mpmath.zetazero at n = 1, 1e5 and 2,001,052.

Usage:
  python3 -I data_checks.py PLATT_DIR OUT.json [ZEROS6.txt] [LMFDB_1e5.txt]
"""
import os
import sys
import json
import struct
import hashlib
import numpy as np
import mpmath as mp

PD, OUT = sys.argv[1], sys.argv[2]
Z6 = sys.argv[3] if len(sys.argv) > 3 else None
LM = sys.argv[4] if len(sys.argv) > 4 else None
mp.mp.dps = 30


def read(path):
    """Zeros of a file as (t0 + offset) with the numbers N0 + 1 ... N1 and
    the exact offsets kept apart for differences."""
    raw = open(path, 'rb').read()
    nb, = struct.unpack_from('<Q', raw, 0)
    pos = 8
    T0, OFF, N = [], [], []
    for _ in range(nb):
        t0, t1, n0, n1 = struct.unpack_from('<ddQQ', raw, pos)
        pos += 32
        n = n1 - n0
        r = np.frombuffer(raw, np.uint8, 13 * n, pos).reshape(n, 13)
        pos += 13 * n
        inc = r[:, :8].copy().view('<u8')[:, 0] * 2.0**-101 + r[:, 8:12].copy().view('<u4')[:, 0] * 2.0**-37 + r[:, 12] * 2.0**-5
        T0.append(np.full(n, t0))
        OFF.append(np.cumsum(inc))
        N.append(np.arange(n0 + 1, n1 + 1))
    return np.concatenate(T0), np.concatenate(OFF), np.concatenate(N)


md5 = {}
for line in open(os.path.join(PD, 'md5.txt')):
    h, f = line.split()
    md5[f.lstrip('*')] = h
files = sorted([f for f in os.listdir(PD) if f.startswith('zeros_') and f.endswith('.dat')], key=lambda f: int(f[6:-4]))
o = dict(files=len(files), md5_ok=0, boundaries=[], S_max=0.0, zeros=0)
prev = None
first = []          # zeros 1 ... 2,001,052 for the comparisons
for f in files:
    raw = open(os.path.join(PD, f), 'rb').read()
    o['md5_ok'] += int(hashlib.md5(raw).hexdigest() == md5[f])
    del raw
    t0, off, N = read(os.path.join(PD, f))
    o['zeros'] += len(N)
    for i in (0, len(N) - 1):
        t = mp.mpf(t0[i]) + mp.mpf(off[i])
        S = float(N[i] - mp.mpf(1) / 2 - (mp.siegeltheta(t) / mp.pi + 1))
        o['S_max'] = max(o['S_max'], abs(S))
    if prev is not None and prev[1] == N[0] - 1:
        gap = (t0[0] - prev[2]) + (off[0] - prev[3])
        tm = t0[0] + off[0]
        o['boundaries'].append(dict(file=f, unfolded_spacing=float(gap * np.log(tm / (2 * np.pi)) / (2 * np.pi))))
    if N[0] <= 2_001_052:
        k = N <= 2_001_052
        first.append((t0[k], off[k], N[k]))
    prev = (f, int(N[-1]), float(t0[-1]), float(off[-1]))
    print(f, len(N), 'S_max so far %.3f' % o['S_max'], flush=True)
t = np.concatenate([a + b for a, b, c in first])
n = np.concatenate([c for a, b, c in first])
assert n[0] == 1 and np.all(np.diff(n) == 1)
o['contiguous_from_1'] = int(n[-1])
o['boundary_spacing_range'] = [min(b['unfolded_spacing'] for b in o['boundaries']), max(b['unfolded_spacing'] for b in o['boundaries'])]
if Z6:
    z6 = np.loadtxt(Z6)
    o['odlyzko_zeros6'] = dict(zeros=len(z6), max_abs_diff=float(np.abs(z6 - t[:len(z6)]).max()))
    print('zeros6', o['odlyzko_zeros6'], flush=True)
if LM:
    lm = [mp.mpf(l.split()[1]) for l in open(LM) if l.strip()]
    o['lmfdb_1e5'] = dict(zeros=len(lm), max_abs_diff=float(max(abs(lm[i] - mp.mpf(t[i])) for i in range(len(lm)))))
    print('lmfdb', o['lmfdb_1e5'], flush=True)
o['mpmath_zetazero'] = {}
for k in (1, 100_000, 2_001_052):
    if k <= len(t):
        o['mpmath_zetazero'][str(k)] = float(abs(mp.zetazero(k).imag - mp.mpf(t[k - 1])))
print(json.dumps({k: v for k, v in o.items() if k != 'boundaries'}))
json.dump(o, open(OUT, 'w'), indent=1)
