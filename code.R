# -----------------------------------------------------------------------------
# Funciones para Resolver la Ecuación Funcional en HLS18
# Autor: Diego Caudillo - Abril 2019 dcaudillo@cimat.mx
# -----------------------------------------------------------------------------


# Como usarse:
# La mayoria de las funciones son de uso interno.
# La interfaz de usuario es
#      d.HLS18(E,N,M,Kappa,eta){
# donde E es un arreglo de reales, M y N son las dimensiones de la matriz
# Kappa es el cuarto cumulante de las entradas normalizadas
# y eta>0 es la resolucion deseada
#      La funcion regresa la densidad de HLS18 en un arreglo del mismo 
# tamano que E.


# Las 4 raices del polinomio por lote sin orden alguno
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

# Algoritmo greedy para asegurar estabilidad numerica.
# Intercambio del orden de las raices considerando longitud de linea
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

# Escoge la mejor raiz usando continuidad y cercania a (-z) lejos de R
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

# Densidad Atomo en Cero [Cauchy Escalada]
dirac.cauchy <- function(x,eta)  return( (eta/pi)/( x*x + eta*eta ) )

# Curva de densidad para un nivel eta
d.HLS18 <- function(E,N,M,Kappa,eta){
  z <- complex(real = E,imaginary = eta)
  m <- m.HLS18(z,N,M,Kappa)
  dens <- pmax( (N/M)*Im(m)/pi , 0)
  atom <- dirac.cauchy(E,eta)*(1-1/d)*d
  dens <- pmax(dens - atom,0)     # HLS18 si considera el atomo en cero
  return( dens )
}

