# Tříbodová korelace nul Riemannovy funkce zeta v konečné výšce: domněnka poměrů koš po koši na 147 milionech nul do t = 3·10¹⁰ (studie IV)

**David Hess (MAYAN ALFA V1008 LAB v.34.72) · 2026** · [English version](README.md)

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.24094318.svg)](https://doi.org/10.5281/zenodo.24094318) Zenodo: [doi:10.5281/zenodo.24094318](https://doi.org/10.5281/zenodo.24094318)

Čtvrtá studie řady, po studiích *Berry saturation and the statistics of the Riemann zeros* (Zenodo, [doi:10.5281/zenodo.23815328](https://doi.org/10.5281/zenodo.23815328)), *Beyond the diagonal approximation* (studie II, [doi:10.5281/zenodo.23928017](https://doi.org/10.5281/zenodo.23928017)) a *Low-lying zeros of 795,381 quadratic Dirichlet L-functions* (studie III, [doi:10.5281/zenodo.23991165](https://doi.org/10.5281/zenodo.23991165)). Měří tříbodovou korelaci R₃ nul Riemannovy funkce zeta v konečné výšce a srovnává ji se vzorcem Conreyho a Snaithové (2008) z domněnky poměrů, který obsahuje všechny členy nižšího řádu až po absolutní člen. Nuly spočítal D. Platt a šíří je LMFDB (147,2 milionu nul ve 14 blocích od první nuly do t = 3,06·10¹⁰); k nim tři sady po 10⁴ nulách A. Odlyzka u t = 10¹¹, 10²⁰ a 10²¹. Chyba každého koše je ze 40 souvislých kusů bloku.

| výsledek | hodnota |
|---|---|
| χ² na koš R₃ proti vzorci, 14 bloků (po 784 koších, a, b ≥ 0,2; očekávání (K−1)/(K−3) = 1,054) | 0,905 – 1,145, průměr 1,038, rozptyl 0,066; žádný blok za 3σ |
| totéž pro R₂ (60 košů) | 0,662 – 1,290, průměr 1,016 |
| R₃ proti GUE / vzorci bez aritmetických faktorů / jen párům / CUE(N_ef) v 1. řádu / CUE s volným N | 6,4 – 117 / 5,0 – 90 / 2,7 – 37 / 27 – 446 / 4,7 – 84 |
| R₃ v koši (0,5–0,6)² u t = 3,06·10¹⁰: data / vzorec / GUE | 0,2533 ± 0,0009 / 0,2519 / 0,2806 |
| amplituda λ₃ odchylky od GUE, R_data − R_GUE = λ (R_poměry − R_GUE), 14 bloků | 0,993 – 1,010, chyby 0,2 – 0,8 %; χ² = 13,9 pro 14 proti λ = 1; sklon (5,4 ± 3,4)·10⁻⁴ na jednotku ln(t/2π) |
| pokles amplitudy přes nejnižší blok (L1e5): pozorovaný / předpovězený vzorcem ve výšce každého kusu | −0,391 ± 0,053 / −0,401; χ² = 7,1 pro 14 sklonů |
| prvních 10⁵ nul (srovnání Conreyho a Snaithové 2008): směrodatná odchylka data − vzorec / čistý šum / jimi uvedená hodnota | 0,0271 / 0,0265 ± 0,0007 / asi 0,03; χ² na koš 1,050 |
| Odlyzkovy sady po 10⁴ nulách u 2,7·10¹¹, 1,4·10²⁰, 1,4·10²¹, spolu | λ₃ = 0,62 ± 0,38 (slučitelné se vzorcem i s GUE) |
| nulové modely pipeline: CUE(10), 10⁷ vlastních hodnot, χ² proti přesnému CUE₁₀ / GUE | 1,04 / 7,27 |
| 200 replik CUE(100) po 10⁴ hodnotách: λ₃/σ průměr a směrodatná odchylka | −0,02; 0,98 |
| nezávislý přepočet (bez sdíleného kódu): histogramy, chyby a amplitudy / všechna χ² / vzorec z mpmath bez tabulek při ln(t/2π) ≥ 6 | na zaokrouhlení (10⁻¹⁴) / do 7·10⁻⁶ / do 4·10⁻⁶ (`data/independent_check.json`) |
| oddíl 5 (dodatek ke studii III): kvadratické L-funkce u \|D\| = 10¹², tři pásma G15, G15b, G15c (27 355 L-funkcí); χ² na koš proti domněnce poměrů / proti USp | 0,996, 0,869, 1,211 / 4,48, 3,96, 4,52 |
| spojená hustota v [0, ½]: data / domněnka poměrů / USp (Katzova–Sarnakova limita) | 0,5040 ± 0,0061 / 0,5054 / 0,4105 (15,4 chyby nad USp) |
| přebytek nad USp krát l = ln(\|D\|/π) od B10 (5·10⁴) po tři pásma (10¹²) | 2,52, 2,58, 2,57, 2,59, 2,48 ± 0,16; fit 7 bodů: c/l p = 0,67, mocnina \|D\| p = 0,007, a + c/l: a = 0,002 ± 0,005; domněnka poměrů: přebytek pod 0,05 až u \|D\| ≈ 3,6·10²¹, pod 0,01 u 2,5·10¹⁰¹ |
| pásmové statistiky, úsek/náhodné (pravidlo zapsané před G15c) | pásma 0,40, 0,81, 0,68 proti kontrolám u 10¹⁰ 0,80, 0,69, 0,91: nerozhodnuto |

Vzorec Conreyho a Snaithové i limita GUE jsou známá teorie; přínosem je test koš po koši s chybou každého koše na 1,5·10⁸ nulách ve 14 výškách, ukázání, že rozdíly, které Conrey a Snaithová uvedli, mají velikost šumu počtu trojic, amplituda odchylky od GUE v každém ze 14 bloků s její předpovědí uvnitř bloků, oddělení aritmetické tříbodové části od součinu párových oprav a od efektivní velikosti matice a nulové modely celé pipeline. Shoda je numerickou podporou domněnky, ne jejím důkazem; o Riemannově hypotéze zde nic neplyne.

## Obsah

| složka | co |
|---|---|
| `paper/` | článek `three_point_correlation_en` (PDF, DOCX, Markdown, LaTeX `main.tex`); česká verze v `paper/cs/` |
| `figures/` | obrázky anglicky (`_en`) a česky (`_cz`), PNG 300 dpi a PDF |
| `data/` | `reference_study_iv.json` (všechna čísla oddílů 4.1–4.5), `arrays_study_iv.npz` (histogramy, chyby a modely každého bloku, amplitudy kusů, mřížka modelů ve výšce), exporty CSV (`blocks.csv`, `bins_r3.csv`, `bins_r2.csv`, `pieces.csv`, `odlyzko.csv`, `cs2008_bins.csv`, `grid_models.csv`), `data_checks.json` (oddíl 2), `independent_check.json`; oddíl 5: `reference_quadratic.json` (tabulka 5, fity, křivka domněnky poměrů, tři pásma a jejich kontroly), `reference_quadratic_table5.csv`, `three_bands.json`, `slow_convergence.json`, `quadratic_check.json` a `quadratic/` (výstupy skriptů studie III pro G15, G15b, G15c a kontrolní úseky E15a–E15c, jackknife okenního testu a manifesty s kontrolním součtem SHA-256 každého bloku nul) |
| `code/` | `zeta_zeros.py` (čtení Plattových souborů, rozvinutí, 14 bloků, Odlyzkovy tabulky), `correlations.py` (histogramy R₂ a R₃ ve 40 kusech, reference GUE a CUE, náhodné unitární matice), `ratios_formula.py` (tabulky na přímce 1 + iy, R₂ a R₃ z domněnky poměrů, mřížka výšek), `compute_models.py`, `compute_study_iv.py`, `data_checks.py`, `independent_check.py`, `export_csv.py`, `figures_iv.py`, `stream_histograms.py` a `big_set_amplitude.py` (histogramy čtené po částech a amplituda pro velmi velkou sadu nul v jedné výšce, např. 10⁹ nul u nuly číslo 10²³, kdyby byla k dispozici); oddíl 5 v `code/quadratic/`: `cut_controls.py` (kontrolní úseky E15a, E15b, E15c z E100), `band_windows_jackknife.py`, `three_bands.py`, `slow_convergence.py` (domněnka poměrů u centrálního bodu přes konduktor), `collect_quadratic.py`, `quadratic_check.py` (nezávislá kontrola), `figures_quadratic.py` a v `study_iii_update/` čtyři skripty studie III s rodinami G15, G15b, G15c |

## Data

Nuly nejsou součástí repozitáře. Plattovy nuly: <https://beta.lmfdb.org/data/riemann-zeta-zeros/> (soubory `zeros_<N>.dat` a `md5.txt`); studie používá 30 souborů `zeros_14.dat` … `zeros_27746000.dat` (souvisle nuly 1 … 68 276 991) a trojice u t ≈ 10⁸, 10⁹, 10¹⁰ a 3,06·10¹⁰ uvedené v `code/zeta_zeros.py` (asi 1,8 GB). Odlyzkovy tabulky `zeros3`, `zeros4`, `zeros5` (a `zeros6` pro kontrolu): <https://www-users.cse.umn.edu/~odlyzko/zeta_tables/>. Prvních 10⁵ nul LMFDB v textu používá jen `data_checks.py`.

## Reprodukce

Příkazy a časy jsou v anglickém [README](README.md#reproducing). Oddíl 5 se dá bez nul přepočítat z mezivýsledků v `data/quadratic/` za několik minut (příkazy tamtéž). Náhodná semínka jsou pevná, takže zveřejněná čísla vyjdou přesně.

Výpočty byly vyvinuty a ověřeny v chráněné laboratoři MAYAN ALFA V1008 (interní, verze 34.72; mayanalfa.com), která není součástí repozitáře.

## Licence

Text, obrázky a data: CC BY 4.0 (`LICENSE-CC-BY-4.0.md`). Kód: MIT (`LICENSE`). Shrnutí česky: `LICENCE.cs.md`.

## Citace

Viz `CITATION.cff`. Tato studie: Hess, D. (2026), Three-point correlation of the Riemann zeros at finite height: the ratios conjecture bin by bin on 147 million zeros up to t = 3·10¹⁰, the amplitude of the deviation from GUE, and smaller sets further up — a reproducible numerical study IV, Zenodo, [doi:10.5281/zenodo.24094318](https://doi.org/10.5281/zenodo.24094318). Předchozí studie: Hess, D. (2026), Berry saturation and the statistics of the Riemann zeros from t ≈ 10⁴ to 10²¹: a reproducible numerical study, Zenodo, doi:10.5281/zenodo.23815328; Hess, D. (2026), Beyond the diagonal approximation: pair correlation, form factor and resonances of the Riemann zeros from t ≈ 10⁴ to 10²¹ — a reproducible numerical study II, Zenodo, doi:10.5281/zenodo.23928017; Hess, D. (2026), Low-lying zeros of 795,381 quadratic Dirichlet L-functions up to conductor 10¹⁰: the ratios conjecture, consecutive families and the first zero — a reproducible numerical study III, Zenodo, doi:10.5281/zenodo.23991165.
