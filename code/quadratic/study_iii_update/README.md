# Scripts of study III used in Section 5

Section 5 of study IV uses the code of study III (doi:10.5281/zenodo.23991165) without change, except for the four files in this folder. Copy them over the code of study III and run that code as described in its README.

| file | change against study III |
|---|---|
| `ks_nuly.py` | the generator of the zeros (PARI/GP) with the families G15 ($10^{12} < \|D\| \le 10^{12} + 15\,000$), G15b ($10^{12} + 10^5 < \|D\| \le 10^{12} + 1.15\cdot10^5$) and G15c ($10^{12} + 2\cdot10^5 < \|D\| \le 10^{12} + 2.15\cdot10^5$); `--zasobnik` (largest PARI stack per process in GB) and a retry with a doubled stack after an overflow; `--rodiny g15` computes G15b and G15c |
| `compute_families.py` | G15, G15b and G15c are treated as consecutive families (like B10 … E100), so the band statistics are computed for them |
| `explicit_weights.py` | in a narrow family the logarithm of the conductor is practically constant; its column is left out of the regression when its standard deviation is below 10⁻³ (otherwise it is collinear with the constant) |
| `g15_trends.py` | new: sets G15 into the trends over the conductor of Sections 4.3–4.8 of study III and writes `reference_g15_trends.json` |

Commands (about 12 h per band on 9 cores with up to 4 GB of PARI stack per process):

```
python3 ks_nuly.py --rodiny G15 --jadra 9 --zasobnik 4 --vystup STORE
python3 ks_nuly.py --rodiny g15 --jadra 9 --zasobnik 4 --vystup STORE
python3 -I ks_kontrola.py STORE                                   # SHA-256 of every block against manifest.json
python3 -I compute_families.py ROOT OUT --families G15b,G15c
python3 -I two_level_first_zero.py ROOT OUT/reference_families_G15b_G15c.json OUT --families G15b,G15c
python3 -I explicit_weights.py ROOT OUT --families G15b,G15c
python3 -I layered_recomputation.py ROOT DATA_III OUT --families G15b,G15c
python3 -I g15_trends.py DATA_III OUT_G15 OUT_G15                  # after the same four steps for --families G15
```

`ROOT` is the folder that holds the store, `DATA_III` the data folder of study III. The comments of the generator are in Czech, as in study III.
