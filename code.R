# =============================================================================
# Numerical Solution of the HLS18 Functional Equation
# Local Law for Sparse Sample Covariance Matrices
# =============================================================================
# Author : Diego Caudillo  |  dcaudillo@cimat.mx
# Date   : April 2019
# =============================================================================
#
# CONTEXT
# -------
# Given an N x M sparse random matrix X with iid entries of variance 1/N
# and fourth cumulant Kappa, the empirical spectral distribution of X^T X
# converges (in an appropriate local sense) to a measure whose Stieltjes
# transform m(z) satisfies the degree-4 polynomial equation P4(m,z) = 0.
#
# This file implements:
#   (1)  m.4.HLS18  -- all four roots of P4, unordered
#   (2)  continuity -- greedy reordering to enforce path-continuity of roots
#   (3)  m.HLS18    -- the unique Stieltjes-admissible root
#   (4)  d.HLS18    -- the local spectral density (user-facing interface)
#
# USER INTERFACE
# --------------
#   d.HLS18(E, N, M, Kappa, eta)
#
#   Arguments:
#     E     -- numeric vector of real energies (evaluation points)
#     N, M  -- matrix dimensions  (d = N/M is the aspect ratio)
#     Kappa -- fourth cumulant of the (normalized) matrix entries
#             (for Bernoulli(p)/sqrt(p) entries: Kappa = (1-p)/p)
#     eta   -- imaginary resolution (eta > 0); smaller = finer scale
#
#   Returns:
#     Numeric vector of the same length as E with the density values.
#
# EXAMPLE
#   E   <- seq(-0.5, 5, length.out = 500)
#   rho <- d.HLS18(E, N=2048, M=2048/2.5, Kappa=500, eta=2048^{-0.6})
#   plot(E, rho, type='l')
# =============================================================================


# -----------------------------------------------------------------------------
# (1)  m.4.HLS18
#      Returns all four roots of P4(m, z) = 0 for each z in the input vector.
#      Roots are computed via R's polyroot() and returned with NO guaranteed
#      ordering -- use m.HLS18 for the admissible root.
#
#      Parameters
#        z     : complex vector (typically z = E + i*eta)
#        N, M  : matrix dimensions
#        Kappa : fourth cumulant of normalized entries
#
#      Returns : 4 x length(z) complex matrix (columns = roots per z)
# -----------------------------------------------------------------------------
m.4.HLS18 <- function(z,N,M,Kappa){
  d <- N/M
  D <- 1-1/d
  
  c.1 <-       D
  c.2 <-   N*D*D*Kappa
  c.3 <- 2*N  *D*Kappa
  c.4 <-   N    *Kappa
  
  r <- sapply( z , function(Z){
    t.r <- polyroot( c(
      1,
      Z + c.1,
      Z + c.2,
      Z * c.3,
      Z*Z * c.4) )
    return( t.r)
  }
  )
  
  return(r)
}


# -----------------------------------------------------------------------------
# (2)  continuity
#      Greedy root-reordering to enforce path-continuity between two
#      consecutive root sets A (previous step) and B (current step).
#
#      For each position j, considers all swaps (j, k) with k > j and
#      applies the swap if it reduces the total displacement:
#        |B[j]-A[j]| + |B[k]-A[k]|  -->  |B[k]-A[j]| + |B[j]-A[k]|
#
#      The theoretical backing for this greedy approach is the equivalence
#      between combinatorial and Euclidean distance for algebraic curves
#      (see Tao, "Topics in Random Matrix Theory", 2012).
#
#      Parameters
#        A : complex vector of roots at the previous step
#        B : complex vector of roots at the current step (to be reordered)
#
#      Returns : reordered version of B
# -----------------------------------------------------------------------------
continuity <- function(A , B){
  if( length(A) != length(B) ) return(NULL)
  for( j in 1:(length(A)-1) ){
    for(l in 1:length(B) )
    for(k in (j+1):length(B) ){
      if(  
        (Mod(B[j]-A[j]) + Mod(B[k]-A[k])) # Status Quo
        > 
          (Mod(B[k]-A[j]) + Mod(B[j]-A[k])) # Change
      ){
        a <- B[k]
        B[k] <- B[j]
        B[j] <- a
      }
    }
  }
  return(B)
}


# -----------------------------------------------------------------------------
# (3)  m.HLS18
#      Selects the unique Stieltjes-admissible root of P4 for each z.
#
#      Strategy:
#        - Far from the real line (Im(z) = high), only one root satisfies
#          z * m(z) ~ -1, which is the Stieltjes condition at infinity.
#        - Starting from that stable regime, approach Im(z) = Im(z_target)
#          via exponentially spaced steps, applying continuity() at each step
#          to track the correct root.
#
#      Parameters
#        z    : complex vector of evaluation points
#        N, M : matrix dimensions
#        Kappa: fourth cumulant
#        high : imaginary part considered "far from real line" (default 1e2)
#        Step : multiplicative step size for the exponential ladder (default 2)
#
#      Returns : complex vector of the admissible m(z) values
# -----------------------------------------------------------------------------
m.HLS18 <- function(z,N,M,Kappa,high=1e2,Step=2){
  r <- sapply( z , function(Z){
    
    K <- max(ceiling( (log(high) - log(Im(Z)))/log(Step) ),1)
    etas <- exp( seq(log(Im(Z)),log(high),length.out=K) )[K:1]
    ZZ <- complex( real=rep(Re(Z),K) , imaginary=etas )
    m <- m.4.HLS18(ZZ,N,M,Kappa)
    for( i in 2:K )  m[,i] <- continuity( m[,i-1] , m[,i] )
    m.fin <- continuity(m[,K],m.4.HLS18(Z,N,M,Kappa))
    R.K <- m.fin[which.min( Mod(m[,1]*ZZ[1] +1) )]
    
    return( R.K )
  }
  )
  
  return(r)
}


# -----------------------------------------------------------------------------
# (4)  dirac.cauchy
#      Cauchy (Poisson) kernel approximating a Dirac delta at 0.
#      Used to remove the atom at zero contributed by the rank deficiency
#      when N > M (i.e. d = N/M > 1 implies a point mass of size 1 - 1/d).
#
#      Parameters
#        x   : real vector
#        eta : half-width of the Cauchy kernel (= imaginary resolution)
#
#      Returns : (eta/pi) / (x^2 + eta^2)
# -----------------------------------------------------------------------------
dirac.cauchy <- function(x,eta)  return( (eta/pi)/( x*x + eta*eta ) )


# -----------------------------------------------------------------------------
# (5)  d.HLS18  [USER INTERFACE]
#      Spectral density of the HLS18 local law at resolution eta.
#
#      Computed as the inverse Stieltjes transform:
#        rho(E) = (N/M) * Im( m(E + i*eta) ) / pi
#      minus the atom at zero (approximated by a scaled Cauchy kernel).
#
#      Parameters
#        E     : real vector of energy values
#        N, M  : matrix dimensions
#        Kappa : fourth cumulant of normalized entries
#        eta   : imaginary resolution (e.g. N^{-0.6} for intermediate scale)
#
#      Returns : non-negative numeric vector of density values (same length as E)
# -----------------------------------------------------------------------------
d.HLS18 <- function(E,N,M,Kappa,eta){
  z <- complex(real = E,imaginary = eta)
  m <- m.HLS18(z,N,M,Kappa)
  dens <- pmax( (N/M)*Im(m)/pi , 0)
  atom <- dirac.cauchy(E,eta)*(1-1/d)*d
  dens <- pmax(dens - atom,0)     # HLS18 si considera el atomo en cero
  return( dens )
}
