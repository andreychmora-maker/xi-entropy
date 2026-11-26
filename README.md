# Asymptotics of ξ-Entropy and Optimal Parameter Selection for Reed–Solomon Masked Cryptosystems

[![Magma](https://img.shields.io/badge/Magma-V2.29%2B-blue.svg)](http://magma.maths.usyd.edu.au/magma/)
[![License: MIT](https://img.shields.io/badge/License-MIT-green.svg)](LICENSE)

Self-contained **Magma V2.29** verification and Monte Carlo simulation suite (`xi_entropy_verification.m`) accompanying the technical note:

> **Andrey L. Chmora**, *"Asymptotics of ξ-Entropy and Optimal Parameter Selection for Reed–Solomon Masked Cryptosystems"*, 2026.

This repository formalizes and empirically validates the asymptotic and finite-length behavior of ephemeral fixed-weight blind noise <b>ξ</b><sub><i>i</i></sub> ~ U(ℋ<sub><i>n</i>, <i>w</i><sub>ξ</sub></sub>) across Reed–Solomon (RS) masked cryptographic primitives, including the monolithic **A-Scheme**, **Stealth ABBE**, and **Fuzzy Covert Time-Lock Puzzles (FCTLP)**.

---

## Core Analytical Results

For a family of linear codes 𝒞<sub><i>n</i></sub>[<i>n</i>, <i>k</i><sub><i>n</i></sub>, <i>d</i><sub><i>n</i></sub>]<sub><i>q</i><sub><i>n</i></sub></sub> with asymptotic rate <i>R</i>, relative distance δ, unique decoding radius ℵ<sub><i>n</i></sub> = ⌊(<i>d</i><sub><i>n</i></sub> − 1)/2⌋, and noise budget share β = <i>w</i><sub>ξ</sub> / ℵ<sub><i>n</i></sub> ∈ (0, 1] (τ = βδ/2), the suite verifies four structural laws:

1. **Static ξ-Entropy & Syndrome Utilization:**
   <p align="center"><i>H</i><sub>ξ</sub>(<i>n</i>, <i>q</i><sub><i>n</i></sub>, <i>w</i><sub>ξ</sub>) = log<sub>2</sub> C(<i>n</i>, <i>w</i><sub>ξ</sub>) + <i>w</i><sub>ξ</sub> log<sub>2</sub>(<i>q</i><sub><i>n</i></sub> − 1) = Θ(<i>n</i> log <i>q</i><sub><i>n</i></sub>), &nbsp;&nbsp;&nbsp; lim<sub><i>n</i>→∞</sub> <i>H</i><sub>ξ</sub> / ((<i>n</i> − <i>k</i><sub><i>n</i></sub>) log<sub>2</sub> <i>q</i><sub><i>n</i></sub>) = β / 2 ≤ 1/2.</p>

2. **Operational Pairwise Collision Security & Alphabet Saturation:**
   Against cross-session linear addition attacks (<b>c</b>′<sub>1</sub> ⊕ <b>c</b>′<sub>2</sub>) with decoder radius θ = lim<sub><i>n</i>→∞</sub> <i>t</i><sub><i>n</i></sub> / <i>n</i>, pairwise unlinkability grows linearly in <i>n</i> but saturates at a **field-independent** hypergeometric support-overlap ceiling as <i>q</i><sub><i>n</i></sub> → ∞ (with optimal exact symbol cancellations <i>c</i><sup>*</sup> = 𝒪(1)):
   <p align="center">lim<sub><i>n</i>→∞</sub> (1/<i>n</i>) <i>H</i><sub>ξ</sub><sup>pair</sup>(<i>n</i>, <i>q</i><sub><i>n</i></sub>, <i>w</i><sub>ξ</sub>, <i>t</i><sub><i>n</i></sub>) = <i>I</i><sub>hyp</sub>(2τ − θ; τ) ≤ <i>h</i><sub>2</sub>(τ).</p>

3. **Guruswami–Sudan Phase Transition & The Low-Rate Trap:**
   In MDS/RS codes (θ<sub>GS</sub> = 1 − √<i>R</i>), positive linear pairwise security requires β > β<sub>crit</sub><sup>GS</sup>(<i>R</i>) = 2(1 − <i>R</i><sup>1/4</sup>) / (1 − <i>R</i>), yielding a universal minimum code rate:
   <p align="center"><i>R</i> > <i>R</i><sub>min</sub> = (<i>u</i><sup>*</sup>)<sup>4</sup> ≈ 0.08737803, &nbsp;&nbsp;&nbsp; where &nbsp; <i>u</i><sup>3</sup> + <i>u</i><sup>2</sup> + <i>u</i> = 1.</p>
   Reducing the code rate from <i>R</i> = 0.75 to <i>R</i> = 0.50 nearly doubles <i>H</i><sub>ξ</sub> yet cuts Guruswami–Sudan (GS) pairwise security roughly in half, confirming that <i>R</i> = 0.75 lies near the global GS optimum <i>R</i><sup>*</sup><sub>GS</sub>(β).

4. **Over-Threshold Scaling (<i>W</i> > <i>t</i><sub><i>n</i></sub>):**
   For super-radius time-lock masks (as in FCTLP), symbol-cancellation grinding requires Δ<sub><i>t</i></sub> = <i>W</i> − <i>t</i><sub><i>n</i></sub> exact cancellations in 𝔽<sub><i>q</i><sub><i>n</i></sub></sub><sup>*</sup>, restoring full Θ(<i>n</i> log <i>q</i><sub><i>n</i></sub>) alphabet scaling.

---

## Repository Contents

* **`xi_entropy_verification.m`** — Self-contained Magma V2.29 verification script organized into six modules:
  * **Part 1:** High-precision root solver for <i>u</i><sup>*</sup> and the universal GS cutoff <i>R</i><sub>min</sub> = (<i>u</i><sup>*</sup>)<sup>4</sup>, plus critical thresholds β<sub>crit</sub><sup>BM</sup>(<i>R</i>) and β<sub>crit</sub><sup>GS</sup>(<i>R</i>).
  * **Part 2:** Golden-section optimization of <i>R</i><sup>*</sup><sub>BM</sub>(β) and <i>R</i><sup>*</sup><sub>GS</sub>(β) across β ∈ [0.65, 0.90] (Table 1).
  * **Part 3:** Alphabet-size decoupling analysis at fixed block length <i>n</i> = 255 across <i>m</i> = log<sub>2</sub> <i>q</i> ∈ {8, …, 32} (Fig. 2).
  * **Part 4:** Exact finite-length double-sum evaluation vs. asymptotic predictors <i>n</i> · <i>I</i><sub>hyp</sub>(2τ − θ; τ) and saddle-point cancellations <i>c</i><sup>*</sup> across 𝔽<sub>2<sup><i>m</i></sup></sub> for <i>m</i> ∈ {8, …, 12} (Table 2 & Fig. 3).
  * **Part 5:** Live algebraic masking and Berlekamp–Massey decoding trials over `ReedSolomonCode(GF(256), 64)`.
  * **Part 6:** Empirical Monte Carlo simulation (30,000 trials) demonstrating the low-rate trap between <i>R</i> = 0.75 (`[255, 192, 64]`) and <i>R</i> = 0.50 (`[255, 128, 128]`).

---

## Summary of Verified Benchmarks

### Selected Finite-Length Profiles (<i>R</i> = 0.75, 𝒞<sub>RS</sub>[2<sup><i>m</i></sup> − 1, <i>k</i>, <i>d</i>]<sub>2<sup><i>m</i></sup></sub>)

| Profile | Code [<i>n</i>, <i>k</i>, <i>d</i>]<sub>2<sup><i>m</i></sup></sub> | (<i>w</i><sub>ξ</sub>, <i>w</i><sub>attr</sub>) | <i>H</i><sub>ξ</sub> (bits) | λ<sub>coll</sub><sup>BM</sup> (asymp.) | λ<sub>coll</sub><sup>GS</sup> (asymp.) | Base-3 LUT (8B) |
| :--- | :---: | :---: | :---: | :---: | :---: | :---: |
| **Compact Hybrid** | [255, 192, 64]<sub>2<sup>8</sup></sub> | (23, 8) | 291.8 | 39.27 (39.82) | 25.72 (24.63) | 51.3 KB (L1/L2) |
| **Tight Hybrid** | [255, 192, 64]<sub>2<sup>8</sup></sub> | (25, 6) | 314.3 | 54.62 (58.02) | 39.30 (39.63) | 51.3 KB (L1/L2) |
| **Low-Rate Trap (<i>m</i> = 8)** | [255, 128, 128]<sub>2<sup>8</sup></sub> | (47, 16) | 547.6 | 50.52 (50.94) | 13.60 (11.50) | 51.3 KB (L1/L2) |
| **Standalone 128b** | [1023, 768, 256]<sub>2<sup>10</sup></sub> | (95, 32) | 1401.5 | 166.78 (167.83) | 116.71 (115.38) | 461.3 KB (L2/L3) |
| **Standalone GS-128b** | [1023, 768, 256]<sub>2<sup>10</sup></sub> | (100, 27) | 1467.7 | 208.40 (212.35) | 153.33 (153.38) | 461.3 KB (L2/L3) |
| **High-Margin** | [2047, 1536, 512]<sub>2<sup>11</sup></sub> | (191, 64) | 3011.7 | 337.47 (338.50) | 239.71 (238.11) | 1.35 MB (L3/GS) |

### Live Algebraic & Monte Carlo Verification (𝔽<sub>256</sub>)

| Experiment | Configuration | Trials | Empirical Result | Theoretical Benchmark |
| :--- | :--- | :---: | :---: | :---: |
| **In-Cohort Blind BM Decode** | [255, 192, 64]<sub>256</sub>, <i>d</i><sub>H</sub> = 8, <i>w</i><sub>ξ</sub> = 23 | 55,000 | **55,000 / 55,000** exact | 100% deterministic (wt ≤ 31) |
| **Cross-Session Frame Sum** | <b>c</b>′<sub>1</sub> ⊕ <b>c</b>′<sub>2</sub>, BM radius ℵ = 31 | 54,999 | **0 / 54,999** valid decodes | Pr ≈ 2<sup>−39.27</sup> |
| **Monte Carlo (<i>R</i> = 0.75)** | [255, 192, 64]<sub>256</sub>, <i>w</i><sub>ξ</sub> = 23, <i>J</i> = 34 | 30,000 | μ̂ = 43.91, min = 38, **0** GS hits | μ = 43.92, 𝔼[GS] ≈ 5.4 × 10<sup>−4</sup> |
| **Monte Carlo (<i>R</i> = 0.50)** | [255, 128, 128]<sub>256</sub>, <i>w</i><sub>ξ</sub> = 47, <i>J</i> = 75 | 30,000 | μ̂ = 85.30, min = 74, **2** GS hits | μ = 85.30, 𝔼[GS] ≈ 2.42 |

---

## Reproduction

### Local Magma V2.29+
```bash
magma xi_entropy_verification.m
```
Full execution with `N_LIVE := 55000` completes in **`53.270 s`**.

### Magma Online Calculator
Paste the contents of `xi_entropy_verification.m` directly into the [Magma Online Calculator](http://magma.maths.usyd.edu.au/calc/) and click **Submit** (fits within the 60-second sandbox limit; optionally set `N_LIVE := 10000;` in Part 5 for a ~15-second run).

---

## Related Works

* **The A-Scheme:** A. L. Chmora, *"The A-Scheme for Information-Theoretic Covert Transport"*, Zenodo (2026). DOI: [`10.5281/zenodo.23574958`](https://doi.org/10.5281/zenodo.23574958)
* **Stealth ABBE:** A. L. Chmora, *"Stealth ABBE: Carry-Free Homomorphic Packing and Blind Reed–Solomon Decoding"*, Zenodo (2026). DOI: [`10.5281/zenodo.23677113`](https://doi.org/10.5281/zenodo.23677113)
* **Fuzzy Covert Time-Lock Puzzles:** A. L. Chmora, *"Fuzzy Covert Time-Lock Puzzles via Sequential Noise Stripping"*, Zenodo (2026). DOI: [`10.5281/zenodo.24382646`](https://doi.org/10.5281/zenodo.24382646)

## Citation

```bibtex
@misc{chmora2026xientropy,
  author       = {Chmora, Andrey L.},
  title        = {Asymptotics of $\xi$-Entropy and Optimal Parameter Selection for Reed--Solomon Masked Cryptosystems},
  year         = {2026},
  publisher    = {Zenodo},
  doi          = {10.5281/zenodo.XXXXXXX},
  url          = {[https://doi.org/10.5281/zenodo.XXXXXXX](https://doi.org/10.5281/zenodo.XXXXXXX)}
}
```
