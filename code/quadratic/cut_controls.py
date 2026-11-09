#!/usr/bin/env python3
"""Control stretches for Section 5: three parts of the family E100 of study III (|D| near 1e10) with the same width as the bands at 1e12 (15,000).
E15a: 1e10 + 20,000 <= |D| < 1e10 + 35,000; E15b: 1e10 + 95,000 ... + 110,000; E15c: 1e10 + 185,000 ... + 200,000.
The blocks of E100 (200 discriminants each) are copied unchanged with their SHA-256 checksums; each control gets its own store (manifest.json + folder),
which the scripts of study III then process like any family (--families E15a, E15b, E15c).
usage: python3 -I cut_controls.py E100_STORE OUTPUT_FOLDER      (E100_STORE: the folder with manifest.json and E100/ of study III)"""
import sys, os, json, copy, shutil, hashlib
SRC, OUT = sys.argv[1:3]
CONTROLS = {'E15a': (10000020000, 10000035000), 'E15b': (10000095000, 10000110000), 'E15c': (10000185000, 10000200000)}
M = json.load(open(os.path.join(SRC, 'manifest.json'))); E = M['rodiny']['E100']
for name, (lo, hi) in CONTROLS.items():
    blocks = {k: v for k, v in E['bloky'].items() if lo <= int(k) < hi}
    assert len(blocks) == (hi - lo)//E['spec']['blok'], (name, len(blocks))
    m = copy.deepcopy({k: v for k, v in M.items() if k != 'rodiny'})
    m['rodiny'] = {name: dict(spec=dict(typ='rozsah', lo=lo, hi=hi, blok=E['spec']['blok'], popis='control stretch of E100, width 15 000'), bloky=blocks)}
    d = os.path.join(OUT, name + '_set'); os.makedirs(os.path.join(d, name), exist_ok=True)
    for k, b in blocks.items():
        s = os.path.join(SRC, 'E100', b['soubor']); t = os.path.join(d, name, b['soubor'])
        assert hashlib.sha256(open(s, 'rb').read()).hexdigest() == b['sha256'], s
        shutil.copyfile(s, t)
    json.dump(m, open(os.path.join(d, 'manifest.json'), 'w'), ensure_ascii=False, indent=1)
    print(name, len(blocks), 'blocks', sum(b['nL'] for b in blocks.values()), 'L-functions')
