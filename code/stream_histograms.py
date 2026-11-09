"""Histograms of R2 and R3 in K pieces for very large sets of zeros, read in chunks.

For sets that do not fit in memory (e.g. 1e9 zeros near the zero number 1e23).
The result is identical to correlations.histograms_in_pieces on the whole
array: the pieces are defined by the total number of left ends, which a first
pass counts, and the last 2 * NEIGH spacings of every chunk are carried over
to the next one.

Input: one or more text files (plain or .gz) with one zero per line, in
Odlyzko's form 'offset' with the base given in a header line
'... gamma - <base> ...' or by --base.  Heights near 1e22 do not fit in a
float64 with 1e-9 precision, so every line is split at the decimal point:
the integer and the fractional part are kept apart and differences are exact
to the accuracy of the table.

Usage:
  python3 -I stream_histograms.py OUT.npz FILE [FILE ...] [--base B] [--K 40] [--chunk 20000000]
Writes r2k[K, 60], r3k[K, 30, 30], count[K], n, base, t_first, t_last, t_median_estimate.
Memory: about 0.2 GB per 1e7 zeros of chunk.
"""
import gzip
import re
import sys
import numpy as np

sys.path.insert(0, __file__.rsplit('/', 1)[0] if '/' in __file__ else '.')
from correlations import NEIGH, AMAX, N2, D2, N3, D3

TWO_PI = 2 * np.pi


def opener(path):
    return gzip.open(path, 'rt') if path.endswith('.gz') else open(path)


def header_base(path):
    with opener(path) as fh:
        head = fh.read(2000)
    m = re.search(r'gamma - (\d+)', head)
    return int(m.group(1)) if m else None


def zero_lines(paths, chunk):
    """Yield (integer parts int64, fractional parts float64) in chunks."""
    buf = []
    for p in paths:
        with opener(p) as fh:
            for line in fh:
                tok = line.strip()
                if not tok or not re.fullmatch(r'\d+\.\d+', tok):
                    continue
                buf.append(tok)
                if len(buf) == chunk:
                    yield split(buf)
                    buf = []
    if buf:
        yield split(buf)


def split(tokens):
    a = np.char.partition(np.array(tokens), '.')
    ip = a[:, 0].astype(np.int64)
    fp = np.char.add('0.', a[:, 2]).astype(np.float64)
    return ip, fp


def count(paths, chunk):
    return sum(len(ip) for ip, fp in zero_lines(paths, chunk))


def accumulate(paths, base, K=40, chunk=20_000_000):
    nz = count(paths, chunk)
    n = (nz - 1) - 2 * NEIGH                       # left ends, as in histograms_in_pieces
    h2 = np.zeros(K * N2)
    h3 = np.zeros(K * N3 * N3)
    over = 0
    carry_s = np.zeros(0)
    start = 0                                        # global index of the first left end in the buffer
    last = None                                      # (integer part, fraction) of the last zero read
    t_first = t_last = None
    for ip, fp in zero_lines(paths, chunk):
        if last is None:
            t_first = float(base) + float(ip[0]) + float(fp[0])
            dt = (np.diff(ip) + np.diff(fp)).astype(np.float64)
            tm = float(base) + ip[:-1] + fp[:-1] + 0.5 * dt
        else:
            ipx = np.concatenate([[last[0]], ip])
            fpx = np.concatenate([[last[1]], fp])
            dt = (np.diff(ipx) + np.diff(fpx)).astype(np.float64)
            tm = float(base) + ipx[:-1] + fpx[:-1] + 0.5 * dt
        last = (ip[-1], fp[-1])
        t_last = float(base) + float(ip[-1]) + float(fp[-1])
        s = np.concatenate([carry_s, dt * np.log(tm / TWO_PI) / TWO_PI])
        m = len(s) - 2 * NEIGH                      # left ends with all neighbours in the buffer
        if m > 0:
            X = np.concatenate([[0.0], np.cumsum(s)])
            gi = start + np.arange(m)
            piece = np.where(gi < n, (gi * K) // max(n, 1), -1)
            for p in range(1, NEIGH + 1):
                a = X[p:p + m] - X[:m]
                mk = (a < AMAX) & (piece >= 0)
                i = np.nonzero(mk)[0]
                aa = a[mk]
                h2 += np.bincount(piece[i] * N2 + (aa / D2).astype(int), minlength=K * N2)
                ia = (aa / D3).astype(int)
                for q in range(1, NEIGH + 1):
                    b = X[i + p + q] - X[i + p]
                    mb = b < AMAX
                    h3 += np.bincount(piece[i[mb]] * N3 * N3 + ia[mb] * N3 + (b[mb] / D3).astype(int), minlength=K * N3 * N3)
                    if q == NEIGH:
                        over += int(mb.sum())
                if p == NEIGH:
                    over += int(mk.sum())
            start += m
            carry_s = s[m:]
        else:
            carry_s = s
    assert start >= n, (start, n)
    cnt = np.bincount((np.arange(n) * K) // n, minlength=K).astype(float)
    r2k = h2.reshape(K, N2) / (cnt[:, None] * D2)
    r3k = h3.reshape(K, N3, N3) / (cnt[:, None, None] * D3 * D3)
    return dict(r2k=r2k, r3k=r3k, count=cnt, n=n, zeros=nz, overflow=over, base=base, t_first=t_first, t_last=t_last,
                t_median_estimate=0.5 * (t_first + t_last))


if __name__ == '__main__':
    args = sys.argv[1:]
    opts = {}
    for k in ('--base', '--K', '--chunk'):
        if k in args:
            i = args.index(k)
            opts[k] = args[i + 1]
            del args[i:i + 2]
    out, files = args[0], args[1:]
    base = int(opts['--base']) if '--base' in opts else header_base(files[0])
    assert base is not None, 'no base: give --base'
    res = accumulate(files, base, int(opts.get('--K', 40)), int(opts.get('--chunk', 20_000_000)))
    np.savez_compressed(out, **res)
    print('zeros', res['zeros'], 'left ends', res['n'], 'overflow', res['overflow'], 't', res['t_first'], res['t_last'])
