# Exact one-dimensional population integrals when phi1=phi2=1.
# Talents are independent and identically lognormal. No investment grid or
# talent quadrature is used: integrate() evaluates smooth finite integrals.
linear_population_means <- function(beta1,beta2,scaler_ratio,B=50,
    meanlog=0,sdlog=sqrt(log((1+sqrt(2))/2)),rel.tol=1e-9) {
  if(any(!is.finite(c(beta1,beta2,scaler_ratio,B,meanlog,sdlog))) ||
     min(beta1,beta2,scaler_ratio,B,sdlog)<=0 || max(beta1,beta2)>=1)
    stop("Require 0<beta1,beta2<1 and positive finite ratio, B and sdlog.")
  cutoff <- function(x) log(beta2)-log(beta1)-log(scaler_ratio)+
    (beta2-1)*log1p(B-x)+(1-beta1)*log1p(x)
  sigma <- sqrt(2)*sdlog
  mu <- exp(meanlog+sdlog^2/2)
  quad <- function(f) integrate(f,0,B,rel.tol=rel.tol,
                                abs.tol=rel.tol*min(1,B),subdivisions=1000L)$value
  E1 <- quad(function(x) pnorm(cutoff(x)/sigma,lower.tail=FALSE))
  E2 <- quad(function(x) pnorm(cutoff(x)/sigma))
  M1 <- mu*quad(function(x) beta1*exp((beta1-1)*log1p(x))*
    pnorm((cutoff(x)-sdlog^2)/sigma,lower.tail=FALSE))
  M2 <- mu*quad(function(x) beta2*exp((beta2-1)*log1p(B-x))*
    pnorm((cutoff(x)+sdlog^2)/sigma))
  list(E1=E1,E2=E2,M1=M1,M2=M2,snew=(E1/M1)/(E2/M2))
}

linear_population_equilibrium <- function(beta1,beta2,B=50,meanlog=0,
    sdlog=sqrt(log((1+sqrt(2))/2)),rel.tol=1e-9) {
  residual <- function(logratio) {
    ratio <- exp(logratio)
    st <- linear_population_means(beta1,beta2,ratio,B,meanlog,sdlog,rel.tol)
    log(st$snew)-logratio
  }
  # The end limits have opposite signs for 0<beta1,beta2<1. This bracket
  # covers the supplied parameters; fail clearly rather than hide underflow.
  # Centre on the ratio that makes the midpoint stationary for equal talents.
  # This avoids normal-tail underflow for very unequal training exponents.
  centre <- log(beta2/beta1)+(beta2-beta1)*log1p(B/2)
  bracket <- centre+c(-8,8)
  ff <- vapply(bracket,residual,1)
  if(any(!is.finite(ff)) || prod(sign(ff))>0)
    stop("Could not bracket a linear-population equilibrium; rescale parameters.")
  ratio <- exp(uniroot(residual,bracket,tol=1e-10)$root)
  st <- linear_population_means(beta1,beta2,ratio,B,meanlog,sdlog,rel.tol)
  c(s=ratio,E1=st$E1,E2=st$E2,M1=st$M1,M2=st$M2)
}
