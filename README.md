# HLS18 Local Law — Numerical Implementation

Numerical solution of the degree-4 polynomial equation arising from the local spectral law for sparse sample covariance matrices, as described in:

> Reinhard Klette, *Concise Computer Vision: An Introduction into Theory and Algorithms*, Springer, London, 2014.

> **This project is no longer maintained and is provided as-is, without warranty of any kind.**

---

## Background

Given an N × M sparse random matrix X with iid entries of variance 1/N and fourth cumulant κ, the empirical spectral distribution of X^T X converges locally to a measure whose Stieltjes transform m(z) satisfies a degree-4 polynomial P₄(m, z) = 0. The challenge is selecting the physically admissible root consistently across all z = E + iη, for small η near the real line.

---

## Files

| File | Description |
|------|-------------|
| `code.R` | Core solver: polynomial roots, greedy continuity reordering, Stieltjes-admissible root selection, and user-facing density function |
| `missing_functions.R` | Simulation, empirical Stieltjes transform, Marchenko-Pastur comparison, and all plotting functions |

---

## Dependencies

- R (≥ 3.5)
- No packages required for `code.R`
- `missing_functions.R` uses base graphics only

---

## Usage

```r
source("code.R")

# Spectral density at a given resolution
N     <- 2048
M     <- round(N / 2.5)    # aspect ratio d = 2.5
p     <- N^(-0.8)           # sparsity ~ 0.002
Kappa <- (1 - p) / p        # fourth cumulant for Bernoulli(p)/sqrt(p) entries
eta   <- N^(-0.6)           # intermediate resolution

E   <- seq(0, 6, length.out = 500)
rho <- d.HLS18(E, N, M, Kappa, eta)
plot(E, rho, type = 'l')
```

```r
source("missing_functions.R")

# Full simulation comparison (reproduces fig:sim.Progresion)
set.seed(42)
plot.progression(N = 2048)
```

---

## Key Functions

### `code.R`

**`d.HLS18(E, N, M, Kappa, eta)`** — main interface. Returns the local spectral density at energies `E` and resolution `eta`.

**`m.HLS18(z, N, M, Kappa)`** — the Stieltjes-admissible root m(z). Tracks the correct root from a stable high-η regime down to the target η via exponentially spaced steps.

**`continuity(A, B)`** — greedy reordering of root set `B` to minimize total displacement from `A`. Enforces path-continuity of algebraic curves across steps.

**`m.4.HLS18(z, N, M, Kappa)`** — all four roots of P₄, unordered.

### `missing_functions.R`

**`simulate.eigenvalues(N, M, p)`** — simulates eigenvalues of X^T X for a sparse Bernoulli matrix.

**`empirical.density(eigenvalues, E, eta)`** — empirical spectral density at resolution η via the inverse Stieltjes formula.

**`marchenko.pastur(x, d)`** — Marchenko-Pastur density for aspect ratio d (baseline comparison).

**`plot.progression(N, ...)`** — reproduces the multi-resolution comparison figure: histogram, empirical Stieltjes, HLS18 prediction, and Marchenko-Pastur.

**`plot.roots(N, M, Kappa, eta)`** — visualizes the four roots in the complex plane across an energy grid.

**`plot.roots.continuity(N, M, Kappa, eta)`** — shows the effect of `continuity()` on root tracking.

---

## Parameters

| Symbol | Argument | Meaning |
|--------|----------|---------|
| N, M | `N`, `M` | Matrix dimensions; d = N/M is the aspect ratio |
| κ | `Kappa` | Fourth cumulant of normalized entries; for Bernoulli(p)/√p: κ = (1−p)/p |
| η | `eta` | Imaginary resolution; smaller → finer scale, more instability |
| p | `p` | Sparsity (entry nonzero probability); in the paper p = N^{−0.8} |

---

## License

GNU General Public License v2 or later.  
Copyright 2019 Diego de Jesus Caudillo Amador — CIMAT.
