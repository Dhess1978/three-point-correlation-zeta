---
title: "Three-point correlation of the Riemann zeros at finite height: the ratios conjecture bin by bin on 147 million zeros up to $t = 3\\cdot10^{10}$, the amplitude of the deviation from GUE, and smaller sets further up — a reproducible numerical study IV"
author: "David Hess (MAYAN ALFA V1008 LAB v.34.72)"
date: "9.11.2026"
abstract: |
  The fourth study of this series measures the three-point correlation $R_3$ of the zeros of the Riemann zeta function at finite height and compares it with the formula that Conrey and Snaith (2008) derived from the ratios conjecture, which contains all lower-order terms down to the constant term. We use 147.2 million zeros computed by D. Platt (LMFDB), arranged in 14 blocks from the first zero to $t = 3.06\cdot10^{10}$, and estimate the error of every bin from 40 consecutive pieces of each block. (i) Conrey and Snaith compared their formula with the first $10^5$ zeros and reported differences with a standard deviation of about 0.03 and a maximum of 0.164, larger than the error they expected for the formula, $T^{-1/2} \approx 0.003$. Repeated with the error of every bin, the differences have the size of the counting noise of the triples (standard deviation 0.027 against 0.0265 ± 0.0007 for pure noise), and the formula agrees bin by bin, $\chi^2$ per bin 1.050 against the expected 1.054. (ii) In all 14 blocks the formula agrees bin by bin ($\chi^2$ per bin 0.905–1.145, mean 1.038 for 784 bins each), while GUE gives 6.4–117, the formula without its arithmetic factors 5.0–90, a model built from pairs alone 2.7–37, and unitary matrices of an effective size fixed without a free parameter 27–446. At $t = 3.06\cdot10^{10}$, $R_3$ in the bin $(0.5\text{–}0.6)^2$ is still about 31 standard errors below GUE and within 1.6 errors of the formula. (iii) Writing the deviation from GUE as one amplitude, $R_{\mathrm{data}} - R_{\mathrm{GUE}} = \lambda\,(R_{\mathrm{ratios}} - R_{\mathrm{GUE}})$, we find $\lambda_3 = 0.993$–$1.010$ in all blocks with errors of 0.2–0.8 %, no trend with height ($5.4\cdot10^{-4} \pm 3.4\cdot10^{-4}$ per unit of $\ln(t/2\pi)$), and the same for the pair correlation. Inside the lowest blocks the height changes enough that the amplitude of consecutive pieces falls, as the formula averaged over the heights of each piece predicts (slope $-0.391 \pm 0.053$ against $-0.401$ in the lowest block); ignoring this inflates the error of $\lambda$ up to 2.8 times, but changes the error of single bins by at most 0.2 %. (iv) Following a suggestion of S. J. Miller, we test smaller sets much further up: three sets of $10^4$ zeros of A. Odlyzko near $t = 2.7\cdot10^{11}$, $1.4\cdot10^{20}$ and $1.4\cdot10^{21}$ give $\lambda_3 = 0.62 \pm 0.38$ together, consistent with both the formula and GUE; about $1.9\cdot10^5$ zeros near $10^{21}$ would be needed for an error of 0.2. As an addendum to study III, the quadratic Dirichlet L-functions with $10^{12} < |D| \le 10^{12} + 15{,}000$ and two further bands of the same width ($27{,}355$ L-functions, computed for this study) follow the ratios conjecture ($\chi^2$ per bin 0.87–1.21) and lie 15.4 standard errors above the Katz–Sarnak limit near the central point (combined); over eight decades of the conductor the excess falls only as $1/\ln(|D|/\pi)$, as the ratios conjecture predicts, a slowness of the kind Watkins found for cubic twists. All numbers are checked by null models of the pipeline, a sensitivity analysis and an independent recomputation that shares no code with the main one; the agreement is numerical support for a conjecture, not a proof of it.
---

# 1. Introduction

Studies I–III of this series [1–3] measured the pair correlation of the Riemann zeros and its saturation (Berry [4]), the off-diagonal terms of Bogomolny and Keating [5, 6], and the one-level density of quadratic Dirichlet L-functions against the ratios conjecture and the Katz–Sarnak limit. The natural next statistic for the zeros of a single L-function is the three-point correlation $R_3$.

In the limit of large height, $R_3$ of the unfolded zeros is the three-point function of the Gaussian unitary ensemble (GUE), the determinant of the sine kernel; for test functions of restricted support this was proved by Hejhal [7] and by Rudnick and Sarnak [8]. At finite height there are lower-order terms that carry the primes. Bogomolny and Keating [9, 6] derived them semiclassically from the Hardy–Littlewood conjecture; Conrey and Snaith [10] derived them from the ratios conjecture of Conrey, Farmer and Zirnbauer [11] (see also [12]) down to the constant term, showed agreement with Bogomolny and Keating, and compared the result with the first $10^5$ zeros. Bogomolny and Keating later gave a determinantal form that can be evaluated at finite height [13]. The formula of Conrey and Snaith has therefore already been compared with zeros; this study cannot claim a first test.

What we add is scale and a statistical framework. We use 147 million zeros computed by Platt [14] and distributed by the LMFDB [15], in 14 blocks from the first zero up to $t = 3.06\cdot10^{10}$, attach an error to every bin, and compare with $\chi^2$ and with null models that remove parts of the formula. We repeat the comparison of Conrey and Snaith on the first $10^5$ zeros with the error of every bin (Section 4.1), test the formula in all blocks (4.2), summarise the deviation from GUE by one amplitude and follow it over the height (4.3), and test small sets much further up, as S. J. Miller suggested to us (4.4); he also warned that convergence in such problems can be very slow and pointed us to the study of Watkins [21], which Section 5 takes up. Section 4.5 collects the null models and the sensitivity analysis. Section 5 is an addendum to study III in the same spirit: quadratic L-functions at $|D| = 10^{12}$. Section 6 lists the statements and their status, Section 7 separates what is new from what is confirmed, Section 8 discusses the limitations.

# 2. Data

**Zeros of zeta.** D. Platt computed the zeros of $\zeta(\tfrac12 + it)$ up to $t \approx 3.06\cdot10^{10}$ (103.8 billion zeros) with rigorous interval arithmetic [14]; the LMFDB distributes them in binary files with an accuracy of $\pm 2^{-102}$ [15]. We downloaded 30 files (147,474,957 zeros; md5 checksums of all files verified): the zeros number 1 to 68,276,991 without gaps ($t \le 2.98\cdot10^7$) and three consecutive files each near $t \approx 10^8$, $10^9$, $10^{10}$ and $3.06\cdot10^{10}$; the 14 blocks of Table 1 are cut from them. We check the files against each other (consecutive files continue exactly), against the zero counting function ($N(t) = \theta(t)/\pi + 1 + S(t)$ with $|S| < 1$ at every file boundary), against Odlyzko's table of the first 2,001,052 zeros (differences below $3\cdot10^{-9}$) and against the first $10^5$ zeros of the LMFDB (below $1.2\cdot10^{-11}$).

| block | zeros | median $t$ | $N_{\mathrm{eff}}$ |
|------|-------:|------:|------:|
| L1e5 | 99,983 | $4.04\cdot10^4$ | 2.02 |
| Z6a | 899,984 | $3.48\cdot10^5$ | 2.51 |
| Z6b | 999,984 | $8.70\cdot10^5$ | 2.72 |
| P5e6 | 2,999,984 | $1.89\cdot10^6$ | 2.90 |
| S1–S6 | 6 × 10,499,984 | $5.1\cdot10^6$ – $2.8\cdot10^7$ | 3.13–3.52 |
| H1e8 | 16,649,768 | $1.02\cdot10^8$ | 3.82 |
| H1e9 | 18,937,009 | $1.00\cdot10^9$ | 4.35 |
| H1e10 | 21,244,828 | $1.00\cdot10^{10}$ | 4.88 |
| H3e10 | 22,366,293 | $3.06\cdot10^{10}$ | 5.13 |

*Table 1. The 14 blocks of Platt's zeros; S1–S6 are consecutive blocks of 10.5 million spacings each. The column "zeros" counts the left ends of triples: the last 17 zeros of a block, which lack 16 following zeros, are not counted (L1e5 thus uses the first $10^5$ zeros).*

The block L1e5 consists of the same first $10^5$ zeros as used by Conrey and Snaith. $N_{\mathrm{eff}} = \ln(t/2\pi)/\sqrt{12\Lambda}$, $\Lambda = 1.57314\ldots$, is the effective size of a unitary matrix that matches the first-order correction to the spacing distribution [16]; it is used only for a null model.

**Smaller sets further up.** A. Odlyzko published sets of $10^4$ consecutive zeros near the zeros number $10^{12}$, $10^{21}$ and $10^{22}$ (heights $2.7\cdot10^{11}$, $1.4\cdot10^{20}$, $1.4\cdot10^{21}$) [17]; we read them as data and use them in Section 4.4.

# 3. Methods

**Unfolding and histograms.** A spacing $\Delta t$ between consecutive zeros is unfolded as $s = \Delta t\,\ln(t_m/2\pi)/2\pi$ with $t_m$ the midpoint. For every zero (the left end of a triple) we take a second zero at the unfolded distance $a < 3$ and a third zero at the distance $b < 3$ beyond the second, each among the next 8 zeros (in no block is the ninth following zero closer than 3), so that $a$ and $b$ are the two gaps of the triple. $R_3(a, b)$ is the density of triples in bins of $0.1 \times 0.1$ on $[0, 3)^2$, $R_2$ the density of pairs in bins of 0.05 on $[0, 3)$; both are normalised so that independent points give 1.

**Error of a bin.** Each block is divided into $K = 40$ consecutive pieces. The error of a bin is the standard deviation of the 40 piece values about the model of each piece, divided by $\sqrt{K}$. Since this error is itself estimated from 40 values, a bin deviation divided by it follows a $t$ distribution with 39 degrees of freedom, and the expected $\chi^2$ per bin is $(K-1)/(K-3) = 1.054$, not 1; the variance of $z^2$ is 2.41 instead of 2.

**The formula.** $R_3$ and $R_2$ at height $t$ are evaluated from Conrey–Snaith [10], eqs. (2.1)–(2.7), which we re-derived from the ratios recipe of [11]; the arithmetic factors are computed with the primes up to $10^5$ ($A$, $B$) and $3\cdot10^4$ ($Q$, $S$), and the values of $\zeta$ and its logarithmic derivatives on the line $1 + iy$ are tabulated for $y \le 52$ (the lowest height of the grid needs $y$ up to 50). Bin averages use $2 \times 2$ Gauss points for $R_3$ and 4 points for $R_2$. The formula depends on the height; a block spans a range of heights, and inside the lowest blocks this range is large ($\ln(t/2\pi)$ from 0.8 to 9.4 in L1e5). We therefore compute the model on a grid of 56 heights ($\ln(t/2\pi) = 0.75$ to 23) and average it, for every piece, over the heights of the zeros of that piece by linear interpolation; the model of a block is the mean of its pieces. The same is done for the formula of Conrey and Snaith, which is itself integrated over $0 \le t \le T$.

**Mask.** In the bins with $a$ or $b$ between 0.1 and 0.2, $R_3$ is of order $10^{-4}$; in the lowest block most pieces contain no triple there, and the error from the pieces does not hold (Poisson regime). $\chi^2$ is therefore given for $a, b \ge 0.2$ (784 bins); the bins with $a$ or $b$ below 0.1 are empty.

**Amplitude.** To follow the deviation from GUE with one number, we write $R_{\mathrm{data}} - R_{\mathrm{GUE}} = \lambda\,(R_{\mathrm{ratios}} - R_{\mathrm{GUE}})$ and estimate $\lambda$ as a weighted linear statistic over the bins, with the weight $(m - g)/m$ ($m$ the formula, $g$ GUE), separately for each piece. $\lambda = 0$ is GUE, $\lambda = 1$ the formula. Because the formula changes with height inside a block, the piece values are compared with the prediction $p_k$ of the formula for the heights of piece $k$: $\lambda = \sum_k l_k p_k / \sum_k p_k^2$, with the error from the residuals $l_k - \lambda p_k$.

**Null models and checks.** (a) The same code on eigenvalues of random unitary matrices CUE(10) ($10^7$ values), compared with the exact CUE$_{10}$ correlations; (b) the spacings of one block shuffled, which keeps $R_2$ of neighbours but destroys genuine three-point structure; (c) the formula without arithmetic factors ($A = 1$, $B = Q = S = 0$) and a model built from pairs alone ($R_{\mathrm{GUE}}$ times the product of the ratios $R_2^{\mathrm{ratios}}/R_2^{\mathrm{GUE}}$ for the three pairs); (d) unitary matrices of size $N_{\mathrm{eff}}$ to first order; (e) for the amplitude: 200 replicas of CUE(100) with $10^4$ values each and 2,236 sub-blocks of $10^4$ zeros of the highest block; (f) sensitivity to shifts of the zeros, to exact unfolding, to the quadrature and the number of primes, and to $K$. (g) An independent recomputation in Python reads the binary files and Odlyzko's tables itself, unfolds, builds the histograms and evaluates the formula with mpmath without tables; it shares no code with the main computation.

# 4. Results

## 4.1 The comparison of Conrey and Snaith repeated

Conrey and Snaith [10, Section 2.1] compared their formula with the first $10^5$ zeros in unscaled variables. They report for the difference data minus formula a mean of $-0.00127$, a standard deviation of about 0.03, a mean absolute value of 0.0257 and a maximum of 0.164, and expected an error of the formula of order $T^{-1/2} \approx 0.003$; they attributed the larger difference to the small height.

We repeat the comparison on the same zeros, unfolded, in bins of $0.1 \times 0.1$ (841 bins with $a, b \ge 0.1$), with the formula averaged over the heights of the zeros. The difference has a standard deviation of 0.0271, a mean absolute value of 0.0207 and a maximum of 0.090 (Figure 1). Pure noise — the formula plus a normal deviate with the error of each bin, 2,000 simulations — gives $0.0265 \pm 0.0007$, $0.0200$ and a maximum of 0.088–0.129 in 90 % of the simulations. The error of a bin estimated from the pieces (rms 0.0266) agrees with the Poisson error of the number of triples (0.0271), an independent route.

With the error of every bin the formula agrees: $\chi^2$ per bin 1.050 for the 784 bins with $a, b \ge 0.2$, against the expected 1.054; GUE gives 6.37. The differences that Conrey and Snaith saw therefore have the size of the counting noise of the triples in $10^5$ zeros — about ten times the error of the formula. Their maximum 0.164 lies above the maximum of our noise; their bins are not known to us, and finer bins give larger noise and a larger maximum, so we can say only that the size matches.

![Figure 1. The first $10^5$ zeros: standard deviation, mean absolute value and maximum of the difference data − formula as reported by Conrey and Snaith [10] (red), in our repetition with the formula averaged over the heights of the zeros (black) and for pure noise, the formula plus a normal deviate with the error of each bin in 2,000 simulations (light blue; the bar marks the range of the maximum in 90 % of the simulations). Dashed: the expected error of the formula, $T^{-1/2}$.](../figures/01_conrey_snaith_2008_en.png){width=80%}


Two methodological points follow. Evaluating the formula at one height $T$ instead of averaging over the heights of the zeros changes the model by a standard deviation of 0.0068 (maximum 0.028) and raises $\chi^2$ to 1.174; this is not what Conrey and Snaith did — their formula is an integral over heights — and our first guess that it explains their differences is rejected. Including the nearly empty bins with $a$ or $b$ between 0.1 and 0.2 raises $\chi^2$ to 1.31 with the error from the pieces and leaves it at 1.05 with the scatter of the pieces about their mean: the result there depends on how the error is estimated, hence the mask.

## 4.2 The formula in 14 blocks up to $3.06\cdot10^{10}$

| block | $\chi^2_3$ ratios | $\chi^2_3$ GUE | $\chi^2_2$ ratios | $R_3(0.55, 0.55)$ data | formula | GUE |
|------|------:|------:|------:|------:|------:|------:|
| L1e5 | 1.050 | 6.4 | 0.775 | $0.1580 \pm 0.0119$ | 0.1819 | 0.2806 |
| Z6a | 1.093 | 22.1 | 1.209 | $0.2049 \pm 0.0043$ | 0.2027 | 0.2806 |
| Z6b | 1.078 | 19.2 | 1.015 | $0.2091 \pm 0.0043$ | 0.2103 | 0.2806 |
| P5e6 | 0.905 | 44.6 | 1.016 | $0.2177 \pm 0.0023$ | 0.2153 | 0.2806 |
| S1 | 1.089 | 116.9 | 1.238 | $0.2202 \pm 0.0016$ | 0.2211 | 0.2806 |
| S2 | 0.989 | 101.3 | 0.873 | $0.2253 \pm 0.0012$ | 0.2248 | 0.2806 |
| S3 | 0.977 | 89.1 | 1.079 | $0.2262 \pm 0.0014$ | 0.2268 | 0.2806 |
| S4 | 1.078 | 85.7 | 1.235 | $0.2295 \pm 0.0015$ | 0.2281 | 0.2806 |
| S5 | 0.988 | 84.4 | 0.876 | $0.2295 \pm 0.0012$ | 0.2291 | 0.2806 |
| S6 | 1.052 | 78.2 | 0.900 | $0.2286 \pm 0.0011$ | 0.2298 | 0.2806 |
| H1e8 | 0.994 | 94.8 | 0.662 | $0.2347 \pm 0.0012$ | 0.2354 | 0.2806 |
| H1e9 | 0.990 | 69.8 | 1.119 | $0.2433 \pm 0.0013$ | 0.2432 | 0.2806 |
| H1e10 | 1.108 | 55.9 | 0.943 | $0.2499 \pm 0.0012$ | 0.2494 | 0.2806 |
| H3e10 | 1.145 | 51.0 | 1.290 | $0.2533 \pm 0.0009$ | 0.2519 | 0.2806 |

*Table 2. $\chi^2$ per bin for $R_3$ (784 bins, $a, b \ge 0.2$; expected 1.054, three standard deviations $\pm 0.17$ for independent bins) and $R_2$ (60 bins; $\pm 0.60$), and $R_3$ in the bin $(0.5\text{–}0.6)^2$.*

![Figure 2. Left: $\chi^2$ per bin of $R_3$ against the formula in the 14 blocks, with the error of a bin from the scatter of the pieces (grey) and against the model at the height of each piece with the mask $a, b \ge 0.2$ (black); band: the expectation 1.054 with three standard deviations for independent bins. Right: with the same errors and mask, GUE, the formula without its arithmetic factors, the model from pairs alone, CUE($N_{\mathrm{eff}}$) to first order, CUE with a free size and the formula (logarithmic scale).](../figures/02_chi2_blocks_en.png){width=95%}


The formula agrees in every block (Table 2, Figure 2 left); the mean of the 14 values of $\chi^2_3$ is 1.038 with a scatter of 0.066, no block is further than three standard deviations from the expectation, and $R_2$ gives a mean of 1.016 with the scatter expected for 60 bins. The largest single deviation in a block, $|z| = 3.0$–4.9, matches the maximum of 784 values of the $t$ distribution with 39 degrees of freedom (median 3.60, 99 % quantile 4.99 of the exact distribution of the maximum), and the share of bins with $|z| > 3$ is 0.37 % against 0.47 % expected; the large deviations are scattered without a pattern (Figure 6).

GUE fails everywhere ($\chi^2_3$ 6.4–117, $\chi^2_2$ 9.9–207). The deviation from GUE is large and shrinks slowly: in the bin $(0.5\text{–}0.6)^2$, where both gaps of a triple are about half a mean spacing, $R_3$ rises from 0.158 at $t \approx 4\cdot10^4$ to 0.2533 at $3.06\cdot10^{10}$, still about 31 standard errors (10 %) below GUE (0.2806) and 1.6 errors from the formula (Figure 3). For comparison, $R_2$ at the same distance (bin 0.55–0.60) is 3–5 % below GUE in the lowest blocks and 1.2 % at $3.06\cdot10^{10}$. Close triples are much more sensitive to the lower-order terms than close pairs, because $R_3$ near the origin grows as the product of the squares of the distances.

![Figure 3. $R_3$ in the bin $(0.5\text{–}0.6)^2$, both gaps of the triple about half a mean spacing, against the median height of each block: zeros (black, with errors), the formula (blue squares) and GUE (orange). At $t = 3.06\cdot10^{10}$ the data are still about 31 errors (10 %) below GUE.](../figures/03_r3_close_triples_en.png){width=85%}


The parts of the formula matter (Figure 2, right). Without its arithmetic factors the formula gives $\chi^2_3$ 5.0–90; a model built from the pair correlation alone (GUE times the ratio of $R_2$ for each of the three pairs) gives 2.7–37 (same errors and mask as above). The arithmetic in the three-point function is therefore not captured by a product of the pair corrections. Unitary matrices of the effective size $N_{\mathrm{eff}}$ to first order [16] overshoot and are worse than GUE ($\chi^2_3$ 27–446); fitting the matrix size freely, $R_3$ requires a larger size (4.4 → 8.7 over the blocks) than $R_2$ (3.4 → 6.7) and the fit remains poor ($\chi^2_3$ 4.7–84). A finite unitary matrix is not a substitute for the arithmetic lower-order terms of $R_3$ at these heights.

## 4.3 One amplitude for the deviation from GUE

| block | $\lambda_3$ | $\lambda_2$ |
|------|------:|------:|
| L1e5 | $0.998 \pm 0.008$ | $0.998 \pm 0.017$ |
| Z6a | $0.995 \pm 0.005$ | $0.997 \pm 0.008$ |
| Z6b | $1.004 \pm 0.004$ | $1.003 \pm 0.010$ |
| P5e6 | $0.998 \pm 0.004$ | $0.996 \pm 0.008$ |
| S1 | $0.998 \pm 0.003$ | $1.001 \pm 0.004$ |
| S2 | $0.998 \pm 0.002$ | $0.998 \pm 0.004$ |
| S3 | $1.000 \pm 0.003$ | $1.000 \pm 0.005$ |
| S4 | $1.003 \pm 0.003$ | $1.004 \pm 0.005$ |
| S5 | $1.002 \pm 0.003$ | $1.002 \pm 0.006$ |
| S6 | $0.993 \pm 0.004$ | $0.990 \pm 0.007$ |
| H1e8 | $0.999 \pm 0.003$ | $1.002 \pm 0.005$ |
| H1e9 | $0.998 \pm 0.004$ | $0.999 \pm 0.008$ |
| H1e10 | $0.999 \pm 0.005$ | $1.005 \pm 0.010$ |
| H3e10 | $1.010 \pm 0.004$ | $1.018 \pm 0.009$ |

*Table 3. Amplitude of the deviation from GUE in the shape of the formula; 0 = GUE, 1 = ratios conjecture.*

In all 14 blocks $\lambda_3$ lies between 0.993 and 1.010 (Table 3, Figure 4) with errors of 0.2–0.8 %; GUE is excluded by 120 to 436 errors. The largest deviation from 1 is 2.3 errors (H3e10, $1.0096 \pm 0.0041$); $\chi^2$ of all blocks against $\lambda = 1$ is 13.9 for 14. A straight line in $L = \ln(t/2\pi)$ gives a slope of $5.4\cdot10^{-4} \pm 3.4\cdot10^{-4}$, i.e. no trend within the errors. The pair correlation gives the same picture, $\lambda_2 = 0.990$–1.018, largest deviation 2.0 errors.

![Figure 4. Amplitude $\lambda_3$ of the deviation from GUE in the shape of the formula (1: ratios conjecture, 0: GUE) for the 14 Platt blocks (black; enlarged in the inset) and for Odlyzko's three sets of $10^4$ zeros (diamonds) with their combination (star, drawn to the right of the highest set).](../figures/04_amplitude_en.png){width=90%}


**Height inside a block.** The amplitude of consecutive pieces falls inside the lowest blocks — in L1e5 by $0.39 \pm 0.05$ across the block (Figure 5). This is not a failure of the formula but its prediction: the model of a block is an average over a wide range of heights, and the deviation from GUE is larger at the lower end. Averaging the formula over the heights of the zeros of each piece predicts the fall piece by piece: slopes $-0.391 \pm 0.053$ observed against $-0.401$ predicted (L1e5), $-0.227 \pm 0.020$ against $-0.229$ (Z6a), $-0.111 \pm 0.009$ against $-0.109$ (S1), $\chi^2 = 7.1$ for the 14 slopes. Computed as the scatter of the pieces about their mean, the error of $\lambda$ would be inflated by up to 2.8 times (L1e5). For $\chi^2$ of single bins the same correction is negligible: the ratio of the rms errors with and without it is 1.000–1.002 in all blocks. The amplitude adds about 800 bins with the same sign, so the slow change of the model with height adds coherently while the noise adds in quadrature; in a single bin the change is small against the noise.

![Figure 5. Amplitude $\lambda_3$ of the 40 consecutive pieces of the blocks L1e5, Z6a and S1 against the position of the piece in its block, with the prediction of the formula averaged over the heights of the zeros of each piece (blue) and the mean of the pieces (dotted).](../figures/05_height_inside_block_en.png){width=95%}


## 4.4 Smaller sets further up

S. J. Miller suggested to us to test "much smaller sets of zeros much further up", because convergence in such problems can be very slow and lower-order effects can hide at small conductor. For the zeros of zeta the public data far above $3\cdot10^{10}$ are Odlyzko's sets of $10^4$ zeros. With the amplitude estimated as above (Table 4):

| height | zeros | $\lambda_3$ | $\lambda_2$ |
|------|------:|------:|------:|
| $2.7\cdot10^{11}$ | 9,983 | $0.92 \pm 0.46$ | $0.95 \pm 0.78$ |
| $1.4\cdot10^{20}$ | 9,983 | $-0.82 \pm 1.07$ | $-1.75 \pm 1.96$ |
| $1.4\cdot10^{21}$ | 9,983 | $0.50 \pm 0.86$ | $0.77 \pm 1.90$ |
| together | | $0.62 \pm 0.38$ | $0.60 \pm 0.68$ |

*Table 4. Amplitude of the deviation from GUE in Odlyzko's sets of $10^4$ zeros; the formula is evaluated at the median height of each set.*

The three sets are consistent with the formula and with GUE (Figure 4); $10^4$ zeros do not decide. The null models confirm that the errors are not underestimated: 200 replicas of CUE(100) with $10^4$ values each give $\lambda_3/\sigma$ with mean $-0.02$ and standard deviation 0.98 (an independent repetition with another seed $-0.00$ and 1.01), and the highest Platt block divided into 2,236 sub-blocks of $10^4$ zeros gives a mean of 1.010 (the whole block: 1.010) with a scatter of 0.30 against the estimated error of 0.35, so for $10^4$ zeros the estimated error is about 15 % conservative. Near $10^{21}$ an error of 0.2 would need about 186,000 zeros, an error of 0.1 about 745,000. Odlyzko's computation of more than $10^9$ zeros near the zero number $10^{23}$, used by Forrester and Mays [18] and Bornemann, Forrester and Mays [19], would give an error of about 0.003; it is not public.

## 4.5 Null models of the pipeline and sensitivity

The same code on $10^7$ eigenvalues of CUE(10) gives $\chi^2_3 = 1.04$ against the exact CUE$_{10}$ and 7.27 against GUE ($R_2$: 0.82 and 4.25; another seed 0.98 and 7.32): the pipeline reproduces a known three-point structure and distinguishes it from GUE. Shuffling the spacings of block S1 destroys the genuine three-point part: $D = R_3 - R_2 \cdot R_2$ of neighbours is zero within the errors ($\chi^2 = 1.06$ against 1,125 for the unshuffled spacings), and the formula fails on the shuffled spacings ($\chi^2_3 = 1{,}026$). Shifting all zeros by $\pm10^{-6}$ changes $R_3$ by at most 0.04 errors ($\pm10^{-3}$: 1.3); unfolding with the exact $\theta(t)$ instead of the midpoint leaves $R_3$ unchanged; the quadrature of the model (10 height quantiles × Gauss $3 \times 3$) and the primes up to $10^5$ also in $Q$ and $S$ change it by at most 0.14 and 0.03 errors; $K = 20, 40, 80$ give for H3e10 (error from the scatter of the pieces) $\chi^2_3$ = 1.18, 1.13, 1.11 against the expectations 1.12, 1.05, 1.03. The independent recomputation (Section 3 (g)) reproduces the histograms, the errors with the model at the height of the piece and every amplitude to rounding ($10^{-14}$), every $\chi^2$ within $7\cdot10^{-6}$ (relative; its CUE coefficient comes from a finite $N$), and the formula from mpmath without tables within $4\cdot10^{-6}$ at $\ln(t/2\pi) \ge 6$, where 98 % of the zeros of the lowest block lie; below, the linear interpolation of the tables is amplified by $1/L^3$, up to $1.6\cdot10^{-4}$ at $L = 1$, a height used by two zeros. Shuffled with another seed, the spacings of S1 give $\chi^2 = 1.08$ for $D$.

![Figure 6. Left: deviation (data − formula)/error of every bin of $R_3$ in the highest block H3e10 ($t = 3.06\cdot10^{10}$; bins with $a, b \ge 0.2$). Right: distribution of the deviations of the 1,568 bins of L1e5 and H3e10 against the $t$ distribution with 39 degrees of freedom expected for an error from 40 pieces.](../figures/06_bin_deviations_en.png){width=95%}


# 5. Addendum to study III: quadratic L-functions at $|D| = 10^{12}$

Study III [3] measured the low zeros of 795,381 quadratic Dirichlet L-functions up to $|D| \approx 10^{10}$. Following the same suggestion of S. J. Miller, we computed with the script of study III (PARI/GP) all fundamental discriminants with $10^{12} < |D| \le 10^{12} + 15{,}000$: 9,123 L-functions (4,563 even, 4,560 odd), 39,150 zeros up to $x = 4.3$ (11.8 h on one computer, 750 blocks with SHA-256 checksums). Stopple [20] computed the low zeros of 1,003,089 odd characters (negative fundamental discriminants) with $10^{12} \le |d| \le 10^{12} + 3.3\cdot10^6$, a range that contains the odd characters of our three bands, and found that the lower-order terms of the ratios conjecture fit his data much better than those of the explicit formula; even characters at this conductor we have not found in the literature.

The ratios conjecture describes the one-level density ($\chi^2$ per bin 0.996; even 1.047, odd 0.948), the Katz–Sarnak limit USp does not ($\chi^2$ 4.48). The density in $[0, \frac12]$ is $0.5158 \pm 0.0106$, against 0.5054 from the ratios conjecture and 0.4105 for USp — 9.9 standard errors above the limit at $|D| = 10^{12}$. Near the central point the ratios conjecture fits here as well, unlike the quadratic twists of an elliptic curve, where the formula with lower-order terms describes the one-level density away from the central point but not the repulsion of the zeros from it [22]; that repulsion was later modelled by an excised orthogonal ensemble, which mimics the discreteness of the central values [23]. The mean first zero $0.7124 \pm 0.0031$ agrees with USp$(2N_{\mathrm{eff}})$ (0.7150, $N_{\mathrm{eff}} = 5.43$) and is far from the limit (0.7827). The explicit formula needs the dual term of the ratios conjecture also here (slope $0.988 \pm 0.012$ against 0.611 without it). The shift of the first zero of even characters, $-0.0034 \pm 0.0043$, does not decide between the laws $1/l$, $1/l^2$ and a free power fitted in study III (predictions 0.0035, 0.0017, 0.0022).

**Slow convergence.** S. J. Miller warned us that convergence in such families can be very slow and pointed us to the study of Watkins [21]. Zagier and Kramarz found in 1987 (as reported in [21]) that 23.3 % of the even cubic twists of $x^3 + y^3 = 1$ with $m \le 7\cdot10^4$ have a vanishing central value and suggested that a positive proportion does; Watkins extended the data to $m \le 10^7$ and found 22.9 %, 20.5 % and 17.7 % for $m \le 10^5$, $10^6$ and $10^7$, a count growing like $x^{0.935}$ against the expected $x^{5/6}(\log x)^C$, and concluded that the proportion more likely tends to zero. The data were right; the extrapolation from them was most likely not.

The quadratic family shows the same kind of slowness (Table 5, Figure 7). The density in $[0, \frac12]$ exceeds the USp value 0.4105 by 0.36 at $|D| \approx 2\cdot10^3$ and still by 0.094 at $10^{12}$ (three bands combined); the excess multiplied by $l = \ln(|D|/\pi)$ stays at 2.5–2.6 from $|D| \approx 5\cdot10^4$ to $10^{10}$ and is $2.48 \pm 0.16$ at $10^{12}$, so it falls only as $1/l$. Three descriptions of the seven points from $5\cdot10^4$ to $10^{12}$ (B10 … E100 and the three bands at $10^{12}$) can be compared. One parameter, $c/l$, fits well ($c = 2.57 \pm 0.02$, $\chi^2 = 4.0$ for 6 degrees of freedom, $p = 0.67$); a power of $|D|$ with two parameters fits worse (exponent $-0.062 \pm 0.002$, $\chi^2 = 15.8$ for 5, $p = 0.007$); and a free limit, $a + c/l$, gives $a = 0.002 \pm 0.005$. The data thus prefer the shape $1/l$ that the ratios conjecture predicts, but they do not fix the limit, and beyond them the descriptions part by orders of magnitude: at $|D| = 10^{120}$, $c/l$ gives an excess of 0.0094 and the power law $1.6\cdot10^{-8}$. Taken alone, the data do not show where the density goes. The ratios conjecture gives the value at every conductor without a free parameter, and its lower-order terms vanish as $l \to \infty$; evaluated with Gauss–Legendre quadrature (64 nodes, checked by adaptive quadrature), they bring the excess below 0.05 only near $|D| \approx 3.6\cdot10^{21}$ and below 0.01 near $2.5\cdot10^{101}$. The approach to USp is therefore an inference through the ratios conjecture, not an observation. The first zero gives a second example of the trap that Watkins describes. The shift of the mean first zero of even characters followed $0.094/l$ in the families up to $|D| \approx 10^8$ (study III, Section 4.5); the enlargement at $10^{10}$ preferred a faster fall (Section 4.8 of study III; a free power gives $l^{-1.6}$), and G15 ($-0.0034 \pm 0.0043$, against 0.0035, 0.0017 and 0.0022 for the laws $1/l$, $1/l^2$ and the free power) does not decide.

| family | $|D|$ | L-functions | density in $[0, \frac12]$ | ratios conjecture | (data − USp)$\cdot l$ | (ratios − USp)$\cdot l$ |
|------|------:|------:|------:|------:|------:|------:|
| B10 | $5\cdot10^4$ | 54,700 | $0.6718 \pm 0.0050$ | 0.6729 | $2.52 \pm 0.05$ | 2.53 |
| C100 | $1.5\cdot10^6$ | 607,941 | $0.6082 \pm 0.0014$ | 0.6084 | $2.58 \pm 0.02$ | 2.58 |
| D100 | $10^8$ | 60,802 | $0.5591 \pm 0.0043$ | 0.5595 | $2.57 \pm 0.07$ | 2.57 |
| E100 | $10^{10}$ | 121,557 | $0.5291 \pm 0.0030$ | 0.5267 | $2.59 \pm 0.06$ | 2.54 |
| G15 | $10^{12}$ | 9,123 | $0.5158 \pm 0.0106$ | 0.5054 | $2.79 \pm 0.28$ | 2.51 |
| G15b | $10^{12}$ | 9,118 | $0.4964 \pm 0.0104$ | 0.5054 | $2.27 \pm 0.28$ | 2.51 |
| G15c | $10^{12}$ | 9,114 | $0.5003 \pm 0.0105$ | 0.5054 | $2.38 \pm 0.28$ | 2.51 |
| G15–G15c combined | $10^{12}$ | 27,355 | $0.5040 \pm 0.0061$ | 0.5054 | $2.48 \pm 0.16$ | 2.51 |

*Table 5. One-level density of quadratic Dirichlet L-functions in $[0, \frac12]$ over eight decades of the conductor (families of study III and the three bands G15, G15b, G15c at $10^{12}$); USp gives 0.4105. The excess over the limit falls as $1/l$, $l = \ln(|D|/\pi)$, as the ratios conjecture predicts.*

![Figure 7. Left: excess of the density in $[0, \frac12]$ over USp times $l = \ln(|D|/\pi)$ for the families of study III and the three bands at $|D| = 10^{12}$ (diamond: combined), with the ratios conjecture (line) and the fit $c/l$ (dashed). Right: the excess itself up to $|D| = 10^{120}$: the ratios conjecture without a free parameter, $c/l$ and a power of $|D|$ fitted to the data; circles mark where the ratios conjecture brings the excess below 0.05 and 0.01.](../figures/07_slow_convergence_en.png){width=100%}


Study III found that families of consecutive discriminants have a smaller low-frequency power of the density than random samples. G15 is one narrow band, and the band statistics turned out to scatter strongly between bands of this width: three control stretches of the same width at $10^{10}$ give 0.80, 0.69 and 0.91 for the ratio of one stretch to random samples, where G15 gives 0.40. Two further bands of the same width, G15b ($10^{12} + 10^5 < |D| \le 10^{12} + 1.15\cdot10^5$, 9,118 L-functions) and G15c ($10^{12} + 2\cdot10^5 < |D| \le 10^{12} + 2.15\cdot10^5$, 9,114), were computed in the same way (750 blocks each with SHA-256 checksums) and processed by the same scripts (Figure 8). The ratios conjecture describes all three bands ($\chi^2$ per bin 0.996, 0.869, 1.211; USp 4.48, 3.96, 4.52). Combined by a weighted mean, the density in $[0, \frac12]$ is $0.5040 \pm 0.0061$, 15.4 standard errors above USp and 0.2 from the ratios conjecture, and the bands agree among themselves ($\chi^2 = 1.90$ for 2 degrees of freedom). The combined mean first zero, $0.7139 \pm 0.0018$, agrees with USp$(2N_{\mathrm{eff}})$ (0.7150); the combined slope of the explicit formula with the dual term is $0.993 \pm 0.007$; the combined shift of the even characters, $-0.0020 \pm 0.0025$, still does not decide between the three laws (predictions 0.0035, 0.0017, 0.0022). For the band statistics we wrote down a rule before G15c was computed: the hypothesis that bands at $10^{12}$ have a smaller ratio of one stretch to random samples than all controls at $10^{10}$ predicts the mean of the three bands below the smallest control (0.69); it is rejected if the mean lies above it, supported if all three bands lie below it, and undecided otherwise. The three bands give 0.40, 0.81, 0.68 (mean 0.63; Welch's $t = -1.27$ against the controls, $p = 0.29$): undecided. Three bands do not decide whether the band statistics change with the height. The window test with the dual term gives $\chi^2$ 49.5, 19.1, 33.3 (30 coefficients) against 18.6, 15.2, 30.4 for the controls; only G15 lies above what 30 coefficients allow ($p = 0.014$), and one control reaches 30.4.

![Figure 8. The three bands at $|D| = 10^{12}$ (G15, G15b, G15c) against the three control stretches of the same width at $10^{10}$ (E15a, E15b, E15c). Top: density in $[0, \frac12]$ with the ratios conjecture and USp, mean first zero with USp$(2N_{\mathrm{eff}})$, and the ratio of one stretch to random samples. Bottom: $z$ of the classes $D \bmod 840$ against random samples, the share of the reduction explained by the primes, and $\chi^2$ of the window test with the dual term. Dotted: means of the bands and of the controls; Welch's $t$ and $p$ in the titles.](../figures/08_three_bands_en.png){width=100%}


# 6. Statements and their status

| statement | status |
|--------------------------------------------------|------------------------------------------|
| The three-point correlation of the unfolded zeros at finite height is given by the formula of Conrey and Snaith from the ratios conjecture, including all lower-order terms. | Known theory [10, 11], conjectural; consistent here bin by bin in 14 blocks up to $3.06\cdot10^{10}$ ($\chi^2_3$ 0.905–1.145, mean 1.038, expected 1.054). |
| The differences reported by Conrey and Snaith for the first $10^5$ zeros show an error of the formula. | Not supported here: on the same zeros the formula agrees bin by bin ($\chi^2_3 = 1.050$), and differences of the reported size are the counting noise of the triples (standard deviation 0.027 against 0.0265 ± 0.0007); their bins are not known to us. |
| The larger differences of Conrey and Snaith came from evaluating the formula at one height. | Rejected (our first guess): their formula is an integral over heights. |
| GUE describes $R_3$ at these heights. | Rejected: $\chi^2_3$ 6.4–117; at $3.06\cdot10^{10}$ still about 31 errors below GUE in the bin $(0.5\text{–}0.6)^2$. |
| The deviation from GUE is carried by the pair correlation alone, or by the formula without arithmetic factors. | Rejected: $\chi^2_3$ 2.7–37 and 5.0–90. |
| Unitary matrices of an effective size describe the deviation. | Rejected: to first order with $N_{\mathrm{eff}}$ $\chi^2_3$ 27–446; fitted freely, $R_3$ needs a larger size than $R_2$ and the fit stays poor. |
| The amplitude of the deviation from GUE equals 1 at all heights. | Consistent here: $\lambda_3$ 0.993–1.010 (0.2–0.8 %), $\chi^2 = 13.9/14$, slope $(5.4 \pm 3.4)\cdot10^{-4}$ per unit of $\ln(t/2\pi)$. |
| The amplitude changes with height. | Not detected: slope $(5.4 \pm 3.4)\cdot10^{-4}$ per unit of $\ln(t/2\pi)$, 1.6 errors. |
| The fall of the amplitude inside the lowest blocks is a failure of the formula. | Rejected: it is predicted piece by piece by the formula averaged over the heights of each piece. |
| The errors from 40 pieces underestimate the scatter. | Rejected: the deviations of bins follow the $t$ distribution with 39 degrees of freedom (0.37 % with $|z| > 3$ against 0.47 %); for the amplitude, replicas of CUE(100) give $\lambda/\sigma$ with standard deviation 0.98, and 2,236 sub-blocks scatter 0.30 against 0.35. |
| The height inside a block distorted the $\chi^2$ per bin. | Rejected: rms error changes by at most 0.2 %. |
| Sets of $10^4$ zeros near $10^{20}$–$10^{21}$ distinguish the formula from GUE. | Rejected: errors 0.86–1.07 per set; the two sets give $-0.02 \pm 0.67$ together, all three sets $0.62 \pm 0.38$, consistent with both. |
| Quadratic L-functions at $|D| = 10^{12}$ follow the ratios conjecture and are not yet at the USp limit. | Consistent in all three bands: $\chi^2$ per bin 0.87–1.21; combined 15.4 errors above USp near the central point. |
| The excess of the one-level density over USp near the central point falls faster than $1/l$. | Rejected: $l$ times the excess stays at 2.5–2.6 from $5\cdot10^4$ to $10^{10}$ and is $2.48 \pm 0.16$ at $10^{12}$ (Table 5; slope against $l$ $0.0024 \pm 0.0055$, rule fixed before the last band). The approach to USp follows only through the ratios conjecture, which predicts the excess below 0.01 near $|D| \approx 2.5\cdot10^{101}$. |
| The shift of the first zero of even characters falls as $0.094/l$. | Not decided: preferred up to $10^8$ (study III), a faster fall preferred with the data at $10^{10}$, three bands at $10^{12}$ do not separate the laws. |
| The band statistics of consecutive discriminants change with the conductor between $10^{10}$ and $10^{12}$. | Not decided: stretch/random 0.40, 0.81, 0.68 for the three bands at $10^{12}$ against 0.80, 0.69, 0.91 for the controls at $10^{10}$ (rule fixed before the last band was computed). |

# 7. What is new and what is confirmed

**Confirmed.** The formula of Conrey and Snaith [10] and its agreement with the zeros; the lower-order terms of Bogomolny and Keating [9]; the slow approach to GUE; the need for arithmetic beyond an effective matrix size [16, 18].

**New, as far as we have been able to find** (we searched arXiv, the journals and the citations of [10] and [13] on 7–8 November 2026): a bin-by-bin test of the three-point formula with an error for every bin on $1.5\cdot10^8$ zeros at 14 heights up to $3.06\cdot10^{10}$, with null models of the pipeline; the differences reported in [10] shown to have the size of the counting noise of the triples; the amplitude of the deviation from GUE at 0.2–0.8 % in each of the 14 blocks, with its prediction inside blocks; the separation of the arithmetic three-point part from a product of pair corrections and from an effective matrix size; and, for quadratic L-functions, even characters at $|D| = 10^{12}$.

# 8. Limitations and future scope

**What is a theorem, a conjecture, an observation.** The GUE limit of $R_3$ for test functions of restricted support is a theorem [7, 8]; the formula with lower-order terms rests on the ratios conjecture, which is a conjecture; everything measured here is a numerical observation. Agreement supports the conjecture, it does not prove it, and nothing here bears on the Riemann hypothesis: the zeros were computed on the critical line, and what is tested is their statistics.

**Scope of the data.** The heights are those of the available tables: continuous up to $3\cdot10^7$, then four windows of about $2\cdot10^7$ zeros. Between $3\cdot10^{10}$ and $10^{20}$ this study has only the $10^4$ zeros near $2.7\cdot10^{11}$, and above $10^{20}$ only $2 \times 10^4$ zeros. The formula has an error of order $T^{-1/2+\varepsilon}$, far below our errors; at the level of $10^{-4}$, the truncation of the prime products and the tables would have to be refined.

**Slow convergence.** Statements at $|D| = 10^{12}$ or $t = 3\cdot10^{10}$ are statements at a finite conductor or height, not about the limit. Families can mislead at small conductor for a long time, as Watkins showed for cubic twists [21]; for the quadratic family the excess of the density over USp falls only as $1/l$ (Section 5, Table 5), and the limit is reached only through the ratios conjecture.

**Statistics.** The error of a bin comes from 40 pieces; the expected $\chi^2$ is therefore 1.054, and bins are correlated through shared triples, so the scatter of $\chi^2$ between blocks (0.066) is somewhat larger than for independent bins (0.055). Only triples with both gaps below 3 mean spacings are counted. In nearly empty bins the error from pieces fails, and we mask them; a Poisson likelihood would use them.

**Next steps.** (i) Larger sets far up: with the data of [18, 19] near $10^{23}$ the amplitude would be measured with an error of about 0.003. (ii) The four-point correlation [9] and the joint distribution of two consecutive spacings at finite height. (iii) For quadratic L-functions: more bands at $|D| = 10^{12}$ and a band at $10^{14}$. (iv) Applications of the three-point machinery to spectra with few levels (e.g. nuclear resonances, through the ratio of consecutive spacings).

# 9. Reproducibility

The repository (GitHub, archived on Zenodo) contains the code, the data and the figures of this paper. The zeros are not included: Platt's files are downloaded from the LMFDB with their md5 sums [15] and Odlyzko's tables from his pages [17]; the md5 sum of every file is checked before use. `zeta_zeros.py` reads the binary files (keeping the start of every block and the offsets of the zeros apart, so that differences of heights stay exact at $t \approx 3\cdot10^{10}$), unfolds and cuts the 14 blocks; `correlations.py` builds the histograms of $R_2$ and $R_3$ in 40 pieces and the random-matrix references; `ratios_formula.py` tabulates $\zeta$, $\zeta'/\zeta$, $(\zeta'/\zeta)'$, $A$ and $B$ on the line $1 + iy$ with mpmath and evaluates the formula; `compute_models.py` computes the models of every block and the grid of 56 heights (17 minutes on two cores); `compute_study_iv.py` computes every number of Sections 4.1–4.5 (about 25 minutes) and writes `reference_study_iv.json` together with the arrays from which the CSV exports (histograms, errors and models of every block, amplitudes of the pieces, the grid of models) and the figures are made; `data_checks.py` repeats the checks of Section 2. All random seeds are fixed, so the published numbers are reproduced exactly. The independent check `independent_check.py` shares no code with these scripts: it reads the files itself, builds its own histograms, evaluates GUE and the first-order CUE correction by its own quadrature (the latter as $N^2(R_N - R_\infty)$ at $N = 2000$), evaluates the formula at control bins directly with mpmath in the form (2.1)–(2.7) of [10], with the factor $P(x, y)$ of (2.7) instead of our $S(x, y)$, and repeats the null models with other seeds. Section 5 uses the generator of the zeros and the scripts of study III [3] (PARI/GP); the repository contains the four files of study III changed for the families G15, G15b and G15c, the intermediate outputs of these scripts and the SHA-256 checksum of every block of zeros. Text, figures and data are under CC BY 4.0, the code under MIT. The laboratory used for development is not part of the publication.

# References

1. Hess, D. (2026). Berry saturation and the statistics of the Riemann zeros from $t \approx 10^4$ to $10^{21}$: a reproducible numerical study. Zenodo. <https://doi.org/10.5281/zenodo.23815328>
2. Hess, D. (2026). Beyond the diagonal approximation: pair correlation, form factor and resonances of the Riemann zeros from $t \approx 10^4$ to $10^{21}$ — a reproducible numerical study II. Zenodo. <https://doi.org/10.5281/zenodo.23928017>
3. Hess, D. (2026). Low-lying zeros of 795,381 quadratic Dirichlet L-functions up to conductor $10^{10}$: the ratios conjecture, consecutive families and the first zero — a reproducible numerical study III. Zenodo. <https://doi.org/10.5281/zenodo.23991165>
4. Berry, M. V. (1988). Semiclassical formula for the number variance of the Riemann zeros. Nonlinearity 1, 399–407.
5. Bogomolny, E. B. & Keating, J. P. (1996). Gutzwiller's trace formula and spectral statistics: beyond the diagonal approximation. Phys. Rev. Lett. 77, 1472–1475.
6. Bogomolny, E. B. & Keating, J. P. (1996). Random matrix theory and the Riemann zeros II: n-point correlations. Nonlinearity 9, 911–935.
7. Hejhal, D. A. (1994). On the triple correlation of zeros of the zeta function. Int. Math. Res. Not. 1994(7), 293–302.
8. Rudnick, Z. & Sarnak, P. (1996). Zeros of principal L-functions and random matrix theory. Duke Math. J. 81, 269–322.
9. Bogomolny, E. B. & Keating, J. P. (1995). Random matrix theory and the Riemann zeros I: three- and four-point correlations. Nonlinearity 8, 1115–1131.
10. Conrey, J. B. & Snaith, N. C. (2008). Triple correlation of the Riemann zeros. J. Théor. Nombres Bordeaux 20, 61–106.
11. Conrey, J. B., Farmer, D. W. & Zirnbauer, M. R. (2008). Autocorrelation of ratios of L-functions. Commun. Number Theory Phys. 2, 593–636.
12. Conrey, J. B. & Snaith, N. C. (2007). Applications of the L-functions ratios conjectures. Proc. London Math. Soc. 94, 594–646.
13. Bogomolny, E. B. & Keating, J. P. (2013). A method for calculating spectral statistics based on random-matrix universality with an application to the three-point correlations of the Riemann zeros. J. Phys. A 46, 305203.
14. Platt, D. J. (2017). Isolating some non-trivial zeros of zeta. Math. Comp. 86, 2449–2467. <https://doi.org/10.1090/mcom/3198>
15. The LMFDB Collaboration. The L-functions and modular forms database, zeros of the Riemann zeta function. <https://www.lmfdb.org>
16. Bogomolny, E., Bohigas, O., Leboeuf, P. & Monastra, A. G. (2006). On the spacing distribution of the Riemann zeros: corrections to the asymptotic result. J. Phys. A 39, 10743–10754.
17. Odlyzko, A. M. Tables of zeros of the Riemann zeta function. <https://www-users.cse.umn.edu/~odlyzko/zeta_tables/>
18. Forrester, P. J. & Mays, A. (2015). Finite size corrections in random matrix theory and Odlyzko's data set for the Riemann zeros. Proc. R. Soc. A 471, 20150436. <https://doi.org/10.1098/rspa.2015.0436>
19. Bornemann, F., Forrester, P. J. & Mays, A. (2017). Finite size effects for spacing distributions in random matrix theory: circular ensembles and Riemann zeros. Stud. Appl. Math. 138, 401–437. <https://doi.org/10.1111/sapm.12160> (arXiv:1608.04638)
20. Stopple, J. (2009). The quadratic character experiment. Exp. Math. 18, 193–200. (arXiv:0802.4255)
21. Watkins, M. (2004). Rank distribution in a family of cubic twists. arXiv:math/0412427.
22. Huynh, D. K., Keating, J. P. & Snaith, N. C. (2009). Lower order terms for the one-level density of elliptic curve L-functions. J. Number Theory 129, 2883–2902. <https://doi.org/10.1016/j.jnt.2008.12.008>
23. Dueñez, E., Huynh, D. K., Keating, J. P., Miller, S. J. & Snaith, N. C. (2012). A random matrix model for elliptic curve L-functions of finite conductor. J. Phys. A 45, 115207. <https://doi.org/10.1088/1751-8113/45/11/115207>
