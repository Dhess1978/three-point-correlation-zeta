# Three-point correlation of the Riemann zeros at finite height: the ratios conjecture bin by bin on 147 million zeros up to t = 3·10¹⁰ (study IV)

**David Hess (MAYAN ALFA V1008 LAB v.34.72) · 2026** · [Česká verze](README.cs.md)

Zenodo: DOI assigned at release.

Fourth study of the series, after *Berry saturation and the statistics of the Riemann zeros* (Zenodo, [doi:10.5281/zenodo.23815328](https://doi.org/10.5281/zenodo.23815328)), *Beyond the diagonal approximation* (study II, [doi:10.5281/zenodo.23928017](https://doi.org/10.5281/zenodo.23928017)) and *Low-lying zeros of 795,381 quadratic Dirichlet L-functions* (study III, [doi:10.5281/zenodo.23991165](https://doi.org/10.5281/zenodo.23991165)). It measures the three-point correlation R₃ of the zeros of the Riemann zeta function at finite height and compares it with the formula of Conrey and Snaith (2008) from the ratios conjecture, which contains all lower-order terms down to the constant term. The zeros are those computed by D. Platt and distributed by the LMFDB (147.2 million zeros in 14 blocks from the first zero to t = 3.06·10¹⁰) and three sets of 10⁴ zeros of A. Odlyzko near t = 10¹¹, 10²⁰ and 10²¹. The error of every bin comes from 40 consecutive pieces of each block.

| result | value |
|---|---|
| χ² per bin of R₃ against the formula, 14 blocks (784 bins each, a, b ≥ 0.2; expected (K−1)/(K−3) = 1.054) | 0.905 – 1.145, mean 1.038, scatter 0.066; no block beyond 3σ |
| the same for R₂ (60 bins) | 0.662 – 1.290, mean 1.016 |
| R₃ against GUE / formula without arithmetic factors / pairs only / CUE(N_eff) to first order / CUE with free N | 6.4 – 117 / 5.0 – 90 / 2.7 – 37 / 27 – 446 / 4.7 – 84 |
| R₃ in the bin (0.5–0.6)² at t = 3.06·10¹⁰: data / formula / GUE | 0.2533 ± 0.0009 / 0.2519 / 0.2806 |
| amplitude λ₃ of the deviation from GUE, R_data − R_GUE = λ (R_ratios − R_GUE), 14 blocks | 0.993 – 1.010, errors 0.2 – 0.8 %; χ² = 13.9 for 14 against λ = 1; slope (5.4 ± 3.4)·10⁻⁴ per unit of ln(t/2π) |
| fall of the amplitude across the lowest block (L1e5): observed / predicted by the formula at the height of each piece | −0.391 ± 0.053 / −0.401; χ² = 7.1 for the 14 slopes |
| the first 10⁵ zeros (comparison of Conrey and Snaith 2008): standard deviation of data − formula / pure noise / their reported value | 0.0271 / 0.0265 ± 0.0007 / about 0.03; χ² per bin 1.050 |
| Odlyzko's sets of 10⁴ zeros near 2.7·10¹¹, 1.4·10²⁰, 1.4·10²¹, together | λ₃ = 0.62 ± 0.38 (consistent with the formula and with GUE) |
| null models of the pipeline: CUE(10), 10⁷ eigenvalues, χ² against exact CUE₁₀ / GUE | 1.04 / 7.27 |
| 200 replicas of CUE(100) with 10⁴ values: λ₃/σ mean and standard deviation | −0.02, 0.98 |
| independent recomputation (no shared code): histograms, errors and amplitudes / every χ² / the formula from mpmath without tables at ln(t/2π) ≥ 6 | to rounding (10⁻¹⁴) / within 7·10⁻⁶ / within 4·10⁻⁶ (`data/independent_check.json`) |
| Section 5 (addendum to study III): quadratic L-functions at \|D\| = 10¹², three bands G15, G15b, G15c (27,355 L-functions); χ² per bin against the ratios conjecture / against USp | 0.996, 0.869, 1.211 / 4.48, 3.96, 4.52 |
| combined density in [0, ½]: data / ratios conjecture / USp (Katz–Sarnak limit) | 0.5040 ± 0.0061 / 0.5054 / 0.4105 (15.4 errors above USp) |
| excess over USp times l = ln(\|D\|/π) from B10 (5·10⁴) to the three bands (10¹²) | 2.52, 2.58, 2.57, 2.59, 2.48 ± 0.16; fit to 7 points: c/l p = 0.67, power of \|D\| p = 0.007, a + c/l: a = 0.002 ± 0.005; ratios conjecture: excess below 0.05 only near \|D\| ≈ 3.6·10²¹, below 0.01 near 2.5·10¹⁰¹ |
| band statistics, stretch/random (rule fixed before G15c) | bands 0.40, 0.81, 0.68 against the controls at 10¹⁰ 0.80, 0.69, 0.91: undecided |

The formula of Conrey and Snaith and the GUE limit are known theory; the contribution is a bin-by-bin test with an error for every bin on 1.5·10⁸ zeros at 14 heights, the differences reported by Conrey and Snaith shown to have the size of the counting noise of the triples, the amplitude of the deviation from GUE in each of the 14 blocks with its prediction inside blocks, the separation of the arithmetic three-point part from a product of pair corrections and from an effective matrix size, and null models of the whole pipeline. Agreement is numerical support for a conjecture, not a proof; nothing here bears on the Riemann hypothesis.

## Contents

| folder | what |
|---|---|
| `paper/` | the article `three_point_correlation_en` (PDF, DOCX, Markdown, LaTeX `main.tex`); the Czech version in `paper/cs/` |
| `figures/` | figures in English (`_en`) and Czech (`_cz`), PNG 300 dpi and PDF |
| `data/` | `reference_study_iv.json` (every number of Sections 4.1–4.5), `arrays_study_iv.npz` (histograms, errors and models of every block, amplitudes of the pieces, the height grid of models), CSV exports (`blocks.csv`, `bins_r3.csv`, `bins_r2.csv`, `pieces.csv`, `odlyzko.csv`, `cs2008_bins.csv`, `grid_models.csv`), `data_checks.json` (Section 2), `independent_check.json`; Section 5: `reference_quadratic.json` (Table 5, the fits, the curve of the ratios conjecture, the three bands and their controls), `reference_quadratic_table5.csv`, `three_bands.json`, `slow_convergence.json`, `quadratic_check.json`, and `quadratic/` (outputs of the scripts of study III for G15, G15b, G15c and the control stretches E15a–E15c, the window jackknife, and the manifests with the SHA-256 checksum of every block of zeros) |
| `code/` | `zeta_zeros.py` (reading Platt's files, unfolding, the 14 blocks, Odlyzko's tables), `correlations.py` (histograms of R₂ and R₃ in 40 pieces, GUE and CUE references, random unitary matrices), `ratios_formula.py` (tables on the line 1 + iy, R₂ and R₃ from the ratios conjecture, height grid), `compute_models.py`, `compute_study_iv.py`, `data_checks.py`, `independent_check.py`, `export_csv.py`, `figures_iv.py`, `stream_histograms.py` and `big_set_amplitude.py` (histograms read in chunks and the amplitude for a very large set of zeros at one height, e.g. 10⁹ zeros near the zero number 10²³, should such a set become available); Section 5 in `code/quadratic/`: `cut_controls.py` (the control stretches E15a, E15b, E15c of E100), `band_windows_jackknife.py`, `three_bands.py`, `slow_convergence.py` (the ratios conjecture near the central point over the conductor), `collect_quadratic.py`, `quadratic_check.py` (independent check), `figures_quadratic.py`, and in `study_iii_update/` the four scripts of study III with the families G15, G15b, G15c |

## Data

The zeros are not included. Platt's zeros: <https://beta.lmfdb.org/data/riemann-zeta-zeros/> (files `zeros_<N>.dat` and `md5.txt`); the study uses the 30 files `zeros_14.dat` … `zeros_27746000.dat` (contiguous, zeros 1 … 68,276,991) and the triples near t ≈ 10⁸, 10⁹, 10¹⁰ and 3.06·10¹⁰ listed in `code/zeta_zeros.py` (about 1.8 GB). Odlyzko's tables `zeros3`, `zeros4`, `zeros5` (and `zeros6` for a check): <https://www-users.cse.umn.edu/~odlyzko/zeta_tables/>. The first 10⁵ zeros of the LMFDB in text are used only by `data_checks.py`.

## Reproducing

```
pip install -r requirements.txt
python3 -I code/zeta_zeros.py PLATT_DIR WORK                     # md5 of every file, unfolding; about a minute
python3 -I code/compute_models.py WORK 2                          # tables on 1 + iy and the models; about 17 minutes on 2 cores
python3 -I code/compute_study_iv.py WORK ODLYZKO_DIR data         # Sections 4.1-4.5; about 25 minutes, 4 GB RAM
python3 -I code/data_checks.py PLATT_DIR data/data_checks.json ZEROS6 LMFDB_1E5     # Section 2; about 5 minutes
python3 -I code/independent_check.py PLATT_DIR ODLYZKO_DIR data data/independent_check.json   # about 6 minutes
python3 -I code/export_csv.py data
python3 -I code/figures_iv.py data figures
```

Section 5 needs the code and the data of study III (doi:10.5281/zenodo.23991165) with the four files of `code/quadratic/study_iii_update/` copied over it, PARI/GP for the zeros, and the store of the family E100 of study III for the controls (commands for the zeros and the scripts of study III in `code/quadratic/study_iii_update/README.md`; about 12 h per band on 9 cores). Then:

```
python3 -I code/quadratic/cut_controls.py E100_STORE CONTROLS                                   # E15a, E15b, E15c; then the scripts of study III on them
python3 -I code/quadratic/band_windows_jackknife.py STUDY_III_CODE ROOT_G15bc OUT_G15bc/reference_explicit_weights_G15b_G15c.json G15b OUT_G15bc/windows_jackknife_G15b.json   # the same for G15c, and for G15 into OUT_G15
python3 -I code/quadratic/three_bands.py OUT_G15 OUT_G15bc OUT_E15 data/three_bands.json
python3 -I code/quadratic/slow_convergence.py OUT_G15/reference_g15_trends.json data/slow_convergence.json          # Gauss-Legendre, 64 nodes; a few minutes
python3 -I code/quadratic/collect_quadratic.py OUT_G15/reference_g15_trends.json OUT_G15bc/reference_families_G15b_G15c.json data/three_bands.json data/slow_convergence.json data/reference_quadratic.json
python3 -I code/quadratic/quadratic_check.py data/reference_quadratic.json data/quadratic_check.json
python3 -I code/quadratic/figures_quadratic.py data/reference_quadratic.json figures
```

Without the zeros, Section 5 is reproduced from the intermediate outputs in `data/quadratic/` (a few minutes; `three_bands.json`, `reference_quadratic.json` and the table are reproduced byte for byte):

```
python3 -I code/quadratic/three_bands.py data/quadratic/G15 data/quadratic/G15b_G15c data/quadratic/E15 data/three_bands.json
python3 -I code/quadratic/slow_convergence.py data/quadratic/G15/reference_g15_trends.json data/slow_convergence.json
python3 -I code/quadratic/collect_quadratic.py data/quadratic/G15/reference_g15_trends.json data/quadratic/G15b_G15c/reference_families_G15b_G15c.json data/three_bands.json data/slow_convergence.json data/reference_quadratic.json
```

`ODLYZKO_DIR` holds the three tables (any file names; a table is recognised by its header). The random seeds are fixed, so the published numbers are reproduced exactly.

The computations were developed and verified in the protected MAYAN ALFA V1008 laboratory (internal, version 34.72; mayanalfa.com), which is not part of this repository.

## Licence

Text, figures and data: CC BY 4.0 (`LICENSE-CC-BY-4.0.md`). Code: MIT (`LICENSE`).

## Citation

See `CITATION.cff`. Previous studies: Hess, D. (2026), Berry saturation and the statistics of the Riemann zeros from t ≈ 10⁴ to 10²¹: a reproducible numerical study, Zenodo, doi:10.5281/zenodo.23815328; Hess, D. (2026), Beyond the diagonal approximation: pair correlation, form factor and resonances of the Riemann zeros from t ≈ 10⁴ to 10²¹ — a reproducible numerical study II, Zenodo, doi:10.5281/zenodo.23928017; Hess, D. (2026), Low-lying zeros of 795,381 quadratic Dirichlet L-functions up to conductor 10¹⁰: the ratios conjecture, consecutive families and the first zero — a reproducible numerical study III, Zenodo, doi:10.5281/zenodo.23991165.
