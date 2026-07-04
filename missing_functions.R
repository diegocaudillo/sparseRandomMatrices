# =============================================================================
# Simulation and Comparison Functions — HLS18 Local Law
# =============================================================================
# Author : Diego Caudillo  |  dcaudillo@cimat.mx
# =============================================================================
#
# This file complements code.R with the functions needed to:
#   (1)  Simulate sparse random matrices and compute their eigenvalues
#   (2)  Compute the empirical Stieltjes transform from eigenvalues
#   (3)  Evaluate the Marchenko-Pastur law (baseline comparison)
#   (4)  Reproduce the figures in the paper:
#          - fig:raices          (four roots of P4 across the complex plane)
#          - fig:raicesContinuidad (continuity tracking of roots)
#          - fig:sim.Progresion  (density comparison at multiple resolutions)
# =============================================================================
#
# DEPENDENCIES
#   library(ggplot2)   -- for plotting functions
#   source("code.R")   -- for m.4.HLS18, m.HLS18, d.HLS18
# =============================================================================


# =============================================================================
# SECTION 1 — SIMULATION
# =============================================================================

# -----------------------------------------------------------------------------
# kappa.from.p
#   Computes the fourth cumulant Kappa of the normalized entry
#   X_{ij} = B_{ij} / sqrt(p), where B_{ij} ~ Bernoulli(p).
#
#   For this model Var(X_{ij}) = 1, and the excess kurtosis equals (1-p)/p.
#   Since HLS18 uses the cumulant of the entry scaled by 1/sqrt(N), the
#   Kappa passed to d.HLS18 should be this value (independent of N).
#
#   Parameters
#     p : sparsity probability (e.g. N^{-0.8})
#
#   Returns : scalar fourth cumulant
# -----------------------------------------------------------------------------
kappa.from.p <- function(p) {
  return( (1 - p) / p )
}


# -----------------------------------------------------------------------------
# simulate.sparse.matrix
#   Generates an N x M sparse random matrix with iid entries.
#   Each entry is X_{ij} = B_{ij} / sqrt(N * p), where B_{ij} ~ Bernoulli(p),
#   so that Var(X_{ij}) = 1/N.
#
#   Parameters
#     N : number of rows
#     M : number of columns
#     p : entry sparsity (probability of being nonzero)
#
#   Returns : N x M numeric matrix
# -----------------------------------------------------------------------------
simulate.sparse.matrix <- function(N, M, p) {
  entries <- rbinom(N * M, size = 1, prob = p) / sqrt(N * p)
  return( matrix(entries, nrow = N, ncol = M) )
}


# -----------------------------------------------------------------------------
# simulate.eigenvalues
#   Generates eigenvalues of X^T X for a sparse random matrix X.
#   Returns all M eigenvalues (including any zeros from rank deficiency).
#
#   Parameters
#     N : number of rows
#     M : number of columns
#     p : entry sparsity
#
#   Returns : numeric vector of length M with eigenvalues of X^T X
# -----------------------------------------------------------------------------
simulate.eigenvalues <- function(N, M, p) {
  X   <- simulate.sparse.matrix(N, M, p)
  XtX <- crossprod(X)           # M x M, faster than t(X) %*% X
  ev  <- eigen(XtX, symmetric = TRUE, only.values = TRUE)$values
  return( sort(ev) )
}


# =============================================================================
# SECTION 2 — EMPIRICAL STIELTJES TRANSFORM
# =============================================================================

# -----------------------------------------------------------------------------
# empirical.stieltjes
#   Computes the empirical Stieltjes transform of a spectral measure:
#     m_N(z) = (1/M) * sum_i  1 / (lambda_i - z)
#
#   For z = E + i*eta with eta > 0 this is always well-defined.
#   Im(m_N(z)) / pi gives the empirical density at resolution eta.
#
#   Parameters
#     eigenvalues : numeric vector of eigenvalues of X^T X
#     z           : single complex number (evaluation point)
#
#   Returns : complex scalar
# -----------------------------------------------------------------------------
empirical.stieltjes <- function(eigenvalues, z) {
  return( mean( 1 / (eigenvalues - z) ) )
}


# -----------------------------------------------------------------------------
# empirical.density
#   Evaluates the empirical spectral density at resolution eta over a
#   grid of real energies, using the inverse Stieltjes formula:
#     rho_N(E) = Im( m_N(E + i*eta) ) / pi
#
#   Parameters
#     eigenvalues : numeric vector of eigenvalues
#     E           : real vector of energy grid points
#     eta         : imaginary resolution
#
#   Returns : numeric vector of density values (same length as E)
# -----------------------------------------------------------------------------
empirical.density <- function(eigenvalues, E, eta) {
  z    <- complex(real = E, imaginary = eta)
  dens <- sapply(z, function(zz) Im(empirical.stieltjes(eigenvalues, zz)) / pi)
  return( pmax(dens, 0) )
}


# =============================================================================
# SECTION 3 — MARCHENKO-PASTUR LAW
# =============================================================================

# -----------------------------------------------------------------------------
# marchenko.pastur
#   Computes the Marchenko-Pastur density for aspect ratio d = N/M.
#
#   Support : [lambda_-, lambda_+]  where lambda_+/- = (1 +/- sqrt(d))^2
#   Density : rho(x) = sqrt((lambda_+ - x)(x - lambda_-)) / (2 * pi * d * x)
#   Atom    : point mass of size (1 - 1/d) at 0, when d > 1
#
#   Parameters
#     x : numeric vector of evaluation points (>= 0)
#     d : aspect ratio N/M
#
#   Returns : numeric vector of density values (atom at 0 not included)
# -----------------------------------------------------------------------------
marchenko.pastur <- function(x, d) {
  lambda.plus  <- (1 + sqrt(d))^2
  lambda.minus <- (1 - sqrt(d))^2
  inside       <- (x >= lambda.minus) & (x <= lambda.plus) & (x > 0)
  dens         <- numeric(length(x))
  dens[inside] <- sqrt( (lambda.plus - x[inside]) * (x[inside] - lambda.minus) ) /
                  ( 2 * pi * d * x[inside] )
  return( dens )
}


# =============================================================================
# SECTION 4 — PLOTTING FUNCTIONS
# =============================================================================

# -----------------------------------------------------------------------------
# plot.roots
#   Reproduces fig:raices.
#   Plots the four roots of P4 in the complex plane for z = E + i*eta,
#   E ranging over a real interval, at a fixed resolution eta.
#
#   Parameters
#     N, M  : matrix dimensions
#     Kappa : fourth cumulant
#     eta   : imaginary resolution (controls how close to the real line)
#     E.range : numeric vector c(E.min, E.max)
#     n.pts   : number of grid points (default 200)
# -----------------------------------------------------------------------------
plot.roots <- function(N, M, Kappa, eta, E.range = c(-5, 5), n.pts = 200) {
  E  <- seq(E.range[1], E.range[2], length.out = n.pts)
  z  <- complex(real = E, imaginary = eta)
  rs <- m.4.HLS18(z, N, M, Kappa)   # 4 x n.pts matrix

  cols <- c("tomato", "steelblue", "seagreen", "darkorange")
  xlim <- range(Re(rs))
  ylim <- range(Im(rs))

  plot(NULL, xlim = xlim, ylim = ylim,
       xlab = "Re(m)", ylab = "Im(m)",
       main = bquote(paste("Roots of ", P[4], ",  ", eta == .(eta))))
  for (k in 1:4)
    points(Re(rs[k,]), Im(rs[k,]), col = cols[k], pch = 16, cex = 0.4)
  legend("topright", legend = paste("Root", 1:4), col = cols,
         pch = 16, bty = "n")
}


# -----------------------------------------------------------------------------
# plot.roots.continuity
#   Reproduces fig:raicesContinuidad.
#   Plots the imaginary parts of the four roots as a function of E
#   (i) without reordering and (ii) after applying continuity().
#
#   Parameters
#     N, M  : matrix dimensions
#     Kappa : fourth cumulant
#     eta   : target imaginary resolution
#     E0    : fixed real part (single value)
#     high  : starting imaginary value for the descent
#     K     : number of exponential steps
# -----------------------------------------------------------------------------
plot.roots.continuity <- function(N, M, Kappa, eta, E0 = 1.0, high = 1e2, K = 40) {
  etas <- exp( seq(log(eta), log(high), length.out = K) )[K:1]
  ZZ   <- complex(real = E0, imaginary = etas)
  raw  <- m.4.HLS18(ZZ, N, M, Kappa)   # 4 x K, unordered
  ord  <- raw
  for (i in 2:K) ord[,i] <- continuity(ord[,i-1], ord[,i])

  par(mfrow = c(1, 2))
  cols <- c("tomato","steelblue","seagreen","darkorange")

  matplot(etas, t(Im(raw)), type = "l", lty = 1, col = cols,
          xlab = expression(eta), ylab = "Im(root)",
          main = "Without reordering", log = "x")

  matplot(etas, t(Im(ord)), type = "l", lty = 1, col = cols,
          xlab = expression(eta), ylab = "Im(root)",
          main = "After continuity()", log = "x")
  par(mfrow = c(1, 1))
}


# -----------------------------------------------------------------------------
# plot.progression
#   Reproduces fig:sim.Progresion.
#   Compares, for each resolution eta in a list, three curves:
#     - Histogram of simulated eigenvalues
#     - Empirical inverse Stieltjes density
#     - HLS18 theoretical density
#   and optionally the Marchenko-Pastur density as a baseline.
#
#   Parameters
#     N, M  : matrix dimensions
#     p     : sparsity
#     etas  : numeric vector of resolutions to show (default: N^c for c in seq)
#     E     : real energy grid
#     n.sim : number of simulation replicas to average (default 1)
# -----------------------------------------------------------------------------
plot.progression <- function(N, M = round(N / 2.5), p = N^{-0.8},
                              etas  = N^seq(-0.3, -0.9, length.out = 8),
                              E     = seq(0, (1 + sqrt(N/M))^2 * 1.3, length.out = 400),
                              n.sim = 1) {

  Kappa <- kappa.from.p(p)
  d     <- N / M

  # Simulate eigenvalues (average over replicas if n.sim > 1)
  cat("Simulating eigenvalues...\n")
  ev.list <- lapply(seq_len(n.sim), function(i) simulate.eigenvalues(N, M, p))
  ev      <- unlist(ev.list)

  # Theoretical HLS18 densities (computed once per eta)
  n.eta <- length(etas)
  nc    <- ceiling(sqrt(n.eta))
  nr    <- ceiling(n.eta / nc)
  par(mfrow = c(nr, nc), mar = c(3, 3, 2, 1))

  for (i in seq_along(etas)) {
    eta <- etas[i]

    d.emp  <- empirical.density(ev, E, eta)
    d.theo <- d.HLS18(E, N, M, Kappa, eta)
    d.mp   <- marchenko.pastur(E, d)

    ymax <- max(d.emp, d.theo, d.mp) * 1.1

    # Histogram of a single replica
    hist(ev.list[[1]], breaks = 60, freq = FALSE,
         xlim = range(E), ylim = c(0, ymax),
         col  = "grey85", border = "grey70",
         main = bquote(eta == N^{.(round(log(eta)/log(N), 2))}),
         xlab = "E", ylab = "density")

    lines(E, d.emp,  col = "steelblue",  lwd = 1.8)   # empirical Stieltjes
    lines(E, d.theo, col = "tomato",     lwd = 1.8)   # HLS18
    lines(E, d.mp,   col = "seagreen",   lwd = 1.2, lty = 2)  # Marchenko-Pastur

    if (i == 1)
      legend("topright", bty = "n", cex = 0.7,
             legend = c("Empirical Stieltjes", "HLS18", "Marchenko-Pastur"),
             col    = c("steelblue", "tomato", "seagreen"),
             lty    = c(1, 1, 2), lwd = c(1.8, 1.8, 1.2))
  }
  par(mfrow = c(1, 1))
}


# =============================================================================
# SECTION 5 — QUICK DEMO
# =============================================================================
# Uncomment to reproduce fig:sim.Progresion with the paper parameters.
#
# source("code.R")
# set.seed(42)
# N <- 2048
# plot.progression(N)
# =============================================================================
