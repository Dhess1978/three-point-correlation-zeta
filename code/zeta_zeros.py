"""Zeros of the Riemann zeta function used in study IV: reading, unfolding, blocks.

Platt's zeros (LMFDB, files zeros_<N>.dat).  A file starts with the number of
blocks (uint64).  Every block has a header (t0, t1: float64; N0, N1: uint64)
and N1 - N0 records of 13 bytes (uint64, uint32, uint8): the increment from
the previous zero (from t0 for the first zero of the block) in units of
2^-101, split into a low, middle and high part.  The zeros of a block are
t0 + cumulative sum of the increments; their numbers are N0 + 1 ... N1.  Only
the declared blocks are read (zeros_2546000.dat has padding after them).

Odlyzko's tables zeros3/zeros4/zeros5: a text header that gives the base
('Values of gamma - <base>') and the first zero number ('numbered 10^<k>'),
then one zero per line as an offset from the base.

Unfolding: a spacing dt between consecutive zeros becomes
s = dt * ln(t_m / 2 pi) / (2 pi), with t_m the midpoint of the spacing.

Command line:
    python3 -I zeta_zeros.py PLATT_DIR WORK_DIR
writes WORK_DIR/zeros/zeros_<N>_s.npy (unfolded spacings), zeros_<N>_t.npy
(height of the left zero of every spacing) and info.json, after checking the
md5 sum of every file against PLATT_DIR/md5.txt.
"""
import os
import re
import sys
import json
import struct
import hashlib
import numpy as np

REC = np.dtype([('a', '<u8'), ('b', '<u4'), ('c', 'u1')])
TWO_PI = 2 * np.pi


def read_platt(path):
    """Return (t0, offset, blocks): t0 of the block of every zero, the offset of
    the zero from t0 and the list of (N0, N1) of the blocks; t = t0 + offset.
    Keeping t0 and the offset apart keeps differences of heights exact."""
    T0, OFF, NS = [], [], []
    with open(path, 'rb') as f:
        nb = struct.unpack('<Q', f.read(8))[0]
        pos = 8
        for _ in range(nb):
            f.seek(pos)
            t0, t1, n0, n1 = struct.unpack('<ddQQ', f.read(32))
            n = n1 - n0
            r = np.frombuffer(f.read(13 * n), dtype=REC)
            inc = (r['a'].astype(np.float64) * 2.0**-101 + r['b'].astype(np.float64) * 2.0**-37
                   + r['c'].astype(np.float64) * 2.0**-5)
            T0.append(np.full(n, t0))
            OFF.append(np.cumsum(inc))
            NS.append((n0, n1))
            pos += 32 + 13 * n
    return np.concatenate(T0), np.concatenate(OFF), NS


def md5_table(platt_dir):
    out = {}
    for line in open(os.path.join(platt_dir, 'md5.txt')):
        h, f = line.split()
        out[f.lstrip('*')] = h
    return out


def platt_files(platt_dir):
    fs = [f for f in os.listdir(platt_dir) if f.startswith('zeros_') and f.endswith('.dat')]
    return sorted(fs, key=lambda f: int(f[6:-4]))


def prepare(platt_dir, out_dir, verify=True):
    """Unfold every Platt file; write *_s.npy, *_t.npy and info.json."""
    os.makedirs(out_dir, exist_ok=True)
    md5 = md5_table(platt_dir) if verify else {}
    info = {}
    for f in platt_files(platt_dir):
        path = os.path.join(platt_dir, f)
        h = hashlib.md5(open(path, 'rb').read()).hexdigest()
        if verify and h != md5[f]:
            raise ValueError('md5 mismatch: ' + f)
        t0, off, ns = read_platt(path)
        dt = np.diff(t0) + np.diff(off)                 # exact differences of heights
        tm = t0[:-1] + off[:-1] + 0.5 * dt              # midpoint of the spacing
        s = dt * np.log(tm / TWO_PI) / TWO_PI
        assert np.all(s > 0), f
        stem = os.path.join(out_dir, f[:-4])
        np.save(stem + '_s.npy', s.astype(np.float64))
        np.save(stem + '_t.npy', (t0[:-1] + off[:-1]).astype(np.float64))
        info[f] = dict(N0=int(ns[0][0]), N1=int(ns[-1][1]), zeros=int(ns[-1][1] - ns[0][0]),
                       t_first=float(t0[0] + off[0]), t_last=float(t0[-1] + off[-1]), md5=h, mean_s=float(s.mean()))
        print(f, info[f]['zeros'], round(info[f]['mean_s'], 6), flush=True)
    json.dump(info, open(os.path.join(out_dir, 'info.json'), 'w'), indent=1)
    return info


# Blocks of study IV.  The contiguous range (zeros 1 ... 68,276,991, files below
# t = 3e7) is cut by zero number; spacing k joins the zeros number k+1 and k+2.
CONTIGUOUS = (('L1e5', 1, 10**5), ('Z6a', 10**5, 10**6), ('Z6b', 10**6, 2 * 10**6), ('P5e6', 2 * 10**6, 5 * 10**6))
S_START, S_LEN, S_COUNT = 5 * 10**6 - 1, 10_500_000, 6
HIGH = (('H1e8', ('zeros_99146000', 'zeros_101246000', 'zeros_103346000')),
        ('H1e9', ('zeros_997946000', 'zeros_1000046000', 'zeros_1002146000')),
        ('H1e10', ('zeros_9998546000', 'zeros_10000646000', 'zeros_10002746000')),
        ('H3e10', ('zeros_30603746000', 'zeros_30605846000', 'zeros_30607946000')))
BLOCK_NAMES = tuple([b[0] for b in CONTIGUOUS] + ['S%d' % (j + 1) for j in range(S_COUNT)] + [b[0] for b in HIGH])


def load_files(prep_dir, files, info=None):
    """Spacings and left heights of consecutive files; the spacing across the
    boundary of two files (the files continue exactly) is inserted."""
    if info is None:
        info = json.load(open(os.path.join(prep_dir, 'info.json')))
    S, T = [], []
    for i, f in enumerate(files):
        if i:
            ta, tb = info[files[i - 1]]['t_last'], info[f]['t_first']
            assert info[files[i - 1]]['N1'] == info[f]['N0'], (files[i - 1], f)
            S.append(np.array([(tb - ta) * np.log(0.5 * (ta + tb) / TWO_PI) / TWO_PI]))
            T.append(np.array([ta]))
        S.append(np.load(os.path.join(prep_dir, f[:-4] + '_s.npy')))
        T.append(np.load(os.path.join(prep_dir, f[:-4] + '_t.npy')))
    return np.concatenate(S), np.concatenate(T)


def blocks(prep_dir, names=None):
    """List of (name, s, t) for the 14 blocks (or the named ones)."""
    info = json.load(open(os.path.join(prep_dir, 'info.json')))
    order = sorted(info, key=lambda f: int(f[6:-4]))
    want = set(BLOCK_NAMES if names is None else names)
    out = []
    if want & set(BLOCK_NAMES[:4 + S_COUNT]):
        cont = [f for f in order if int(f[6:-4]) < 30_000_000]
        s_all, t_all = load_files(prep_dir, cont, info)
        assert len(s_all) == info[cont[-1]]['N1'] - info[cont[0]]['N0'] - 1
        for nm, lo, hi in CONTIGUOUS:
            if nm in want:
                out.append((nm, s_all[lo - 1:hi - 1], t_all[lo - 1:hi - 1]))
        k = S_START
        for j in range(S_COUNT):
            k1 = min(len(s_all), k + S_LEN)
            if 'S%d' % (j + 1) in want:
                out.append(('S%d' % (j + 1), s_all[k:k1], t_all[k:k1]))
            k = k1
        del s_all, t_all
    for nm, pre in HIGH:
        if nm in want:
            s, t = load_files(prep_dir, [p + '.dat' for p in pre], info)
            out.append((nm, s, t))
    return out


def first_zeros(prep_dir, count=100_000):
    """Spacings and heights of the first zeros (block of Conrey and Snaith)."""
    info = json.load(open(os.path.join(prep_dir, 'info.json')))
    order = sorted(info, key=lambda f: int(f[6:-4]))
    files = []
    for f in order:
        files.append(f)
        if info[f]['N1'] > count + 100:
            break
    s, t = load_files(prep_dir, files, info)
    return s[0:count - 1], t[0:count - 1]


def read_odlyzko(path):
    """Return (base, k, offsets): zeros base + offsets numbered 10^k + 1, ..."""
    txt = open(path).read()
    base = int(re.search(r'gamma - (\d+)', txt).group(1))
    k = int(re.search(r'numbered 10\^(\d+)', txt).group(1))
    g = np.array([float(l) for l in txt.split('\n') if re.fullmatch(r'\s*\d+\.\d+\s*', l)])
    return base, k, g


def unfold_odlyzko(base, g):
    """Unfolded spacings and the median height of an Odlyzko set."""
    tm = base + 0.5 * (g[:-1] + g[1:])
    s = np.diff(g) * np.log(tm / TWO_PI) / TWO_PI
    return s, float(base + np.median(g))


if __name__ == '__main__':
    prepare(sys.argv[1], os.path.join(sys.argv[2], 'zeros'))
