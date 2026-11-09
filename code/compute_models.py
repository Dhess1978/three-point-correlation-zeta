"""Models of study IV from the ratios conjecture (the expensive part, about 30 minutes on 2 cores).

1. tables on the line 1 + iy (ratios_formula.build_tables), if missing;
2. models of every block at 5 height quantiles of the block (as in study IV,
   Section 4.2 'old' chi^2 and the weights of the amplitude): R3 full formula,
   R3 without arithmetic ('zeta'), R3 'pairs only'; R2 at 10 height quantiles;
3. the height grid: R3 (full) and R2 at the 56 heights L = 0.75 ... 23.

Usage:  python3 -I compute_models.py WORK_DIR [PROCESSES]
reads WORK_DIR/zeros (zeta_zeros.py), writes WORK_DIR/tables.npz,
WORK_DIR/models_blocks.npz and WORK_DIR/models_grid.npz.
"""
import os
import sys
import time
import numpy as np
import multiprocessing as mpc

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from zeta_zeros import blocks
from correlations import gue_r3
from ratios_formula import build_tables, Tables, Primes, r3_bins, r3_pairs_bins, r2_bins, height_quantiles, HEIGHT_GRID

W = sys.argv[1].rstrip('/') + '/'
NPROC = int(sys.argv[2]) if len(sys.argv) > 2 else 2
T0 = time.time()
log = lambda *a: print('[%5.0f s]' % (time.time() - T0), *a, flush=True)
if not os.path.exists(W + 'tables.npz'):
    log('tables', build_tables(W + 'tables.npz'))
TAB = Tables(W + 'tables.npz')
PR = Primes(30_000)


def job(arg):
    kind, t = arg
    if kind == 'full':
        return r3_bins(TAB, PR, t, 'full')
    if kind == 'zeta':
        return r3_bins(TAB, PR, t, 'zeta')
    if kind == 'pairs':
        return r3_pairs_bins(TAB, t)
    if kind == 'r2':
        return r2_bins(TAB, t)
    raise ValueError(kind)


if __name__ == '__main__':
    if not os.path.exists(W + 'models_blocks.npz'):
        q5, q10 = {}, {}
        for nm, s, t in blocks(W + 'zeros'):
            q5[nm] = height_quantiles(t, 5)
            q10[nm] = height_quantiles(t, 10)
        log('height quantiles of', len(q5), 'blocks')
        jobs = [(nm, kind, tt) for nm in q5 for kind in ('full', 'zeta', 'pairs') for tt in q5[nm]]
        jobs += [(nm, 'r2', tt) for nm in q10 for tt in q10[nm]]
        with mpc.get_context('fork').Pool(NPROC) as pool:
            res = pool.map(job, [(k, tt) for nm, k, tt in jobs], chunksize=1)
        out = dict(gue=gue_r3(2))
        for nm in q5:
            for kind in ('full', 'zeta', 'pairs', 'r2'):
                out[nm + '_' + kind] = np.mean([r for (n_, k_, tt), r in zip(jobs, res) if n_ == nm and k_ == kind], axis=0)
            out[nm + '_q5'] = q5[nm]
            out[nm + '_q10'] = q10[nm]
        np.savez_compressed(W + 'models_blocks.npz', **out)
        log('block models')
    if not os.path.exists(W + 'models_grid.npz'):
        ts = [2 * np.pi * np.exp(L) for L in HEIGHT_GRID]
        with mpc.get_context('fork').Pool(NPROC) as pool:
            g3 = pool.map(job, [('full', t) for t in ts], chunksize=1)
            g2 = pool.map(job, [('r2', t) for t in ts], chunksize=1)
        np.savez_compressed(W + 'models_grid.npz', L=HEIGHT_GRID, G3=np.array([g.ravel() for g in g3]), G2=np.array(g2))
        log('height grid')
    log('done')
