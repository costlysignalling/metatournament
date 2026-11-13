# Population reference for a focal actor's first talent-informed response.
# X1 is uniform on [0,B], X2=B-X1, independent of the lognormal talents.
# E1=E2=B/2; each M is calculated separately, never set equal by assumption.
# Requires only base R. Linear spoils returns use an analytic expectation.
random_allocation_reference <- function(beta1,beta2,phi1,phi2,B=50,
    meanlog=0,sdlog=sqrt(log((1+sqrt(2))/2)),nquad=64L) {
  pars <- c(beta1,beta2,phi1,phi2,B,sdlog)
  if(any(!is.finite(pars)) || any(pars<=0) || !is.finite(meanlog))
    stop("Exponents, budget and log-scale SD must be finite and positive.")
  nquad <- as.integer(nquad)
  if(length(nquad)!=1L || is.na(nquad) || nquad<2L) stop("nquad must be at least 2.")
  J <- matrix(0,nquad,nquad)
  J[cbind(1:(nquad-1),2:nquad)] <- sqrt((1:(nquad-1))/2)
  J <- J+t(J)
  gh <- eigen(J,symmetric=TRUE)
  talent <- exp(meanlog+sqrt(2)*sdlog*gh$values)
  weights <- gh$vectors[1,]^2
  moment <- function(beta,phi) {
    if(phi==1) return(exp(meanlog+sdlog^2/2)*
      (expm1((beta+1)*log1p(B))/(B*(beta+1))-1))
    integrand <- function(u) {
      training <- expm1(beta*log1p(B*u))
      z <- expm1(phi*log1p(outer(talent,training)))
      ans <- as.vector(crossprod(weights,z))
      if(any(!is.finite(ans))) stop("Population expectation exceeds floating-point range.")
      ans
    }
    integrate(integrand,0,1,rel.tol=1e-9,subdivisions=500L)$value
  }
  M1 <- moment(beta1,phi1); M2 <- moment(beta2,phi2)
  if(any(!is.finite(c(M1,M2))) || min(M1,M2)<=0)
    stop("Invalid population mean spoils weights.")
  list(E1=B/2,E2=B/2,M1=M1,M2=M2,s=M2/M1)
}
