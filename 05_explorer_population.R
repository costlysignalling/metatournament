# Deterministic population moments for the interactive asymmetric explorer.
# Rcpp caches compilation; subsequent calls reuse the loaded functions.
source("05_energy_allocation_quadrature.R")
explorer_population <- local({
  engine <- new.env(parent=globalenv())
  function(beta1,beta2,phi1,phi2,B,fixed_ratio,responsive=TRUE,
           nquad=64L,nx=801L) {
    if(!exists("explorer_moments_cpp",engine,inherits=FALSE)) {
      code <- readLines("05_energy_allocation_population_phispace.cpp",warn=FALSE)
      wrapper <- '
// [[Rcpp::export]]
NumericVector explorer_moments_cpp(double b1,double b2,double p1,double p2,
 double rho,NumericVector aa,NumericVector aw,NumericVector nodes,
 NumericVector weights,double B,double sd,int nx) {
 ConditionalCell cell{b1,b2,p1,p2,B,sd,rho,nx};
 ConditionalState z=cell.state(aa,aw,nodes,weights);
 return NumericVector::create(_["E1"]=z.e,_["E2"]=B-z.e,
 _["M1"]=z.m1,_["M2"]=z.m2,_["s"]=rho);
}'
      Rcpp::sourceCpp(code=paste(c(code,wrapper),collapse="\n"),
        env=engine,cacheDir="05_explorer_cpp_cache")
    }
    sd <- sqrt(log((1+sqrt(2))/2))
    # Integrate the lower-phi activity outside the switch-split integral;
    # for equal phi use the smaller beta. Reverse labels on return.
    swap <- phi1>phi2 || (phi1==phi2 && beta1>beta2)
    bb <- if(swap)c(beta2,beta1) else c(beta1,beta2)
    pp <- if(swap)c(phi2,phi1) else c(phi1,phi2)
    ratio <- if(swap)1/fixed_ratio else fixed_ratio
    gh <- energy_gh(as.integer(nquad),sd); gl <- energy_gl(64L)
    moments <- function(rho) {
      z <- as.list(engine$explorer_moments_cpp(bb[1],bb[2],pp[1],pp[2],
        rho,gh$a,gh$w,gl$x,gl$w,B,sd,as.integer(nx)))
      if(swap) z <- list(E1=z$E2,E2=z$E1,M1=z$M2,M2=z$M1,s=1/z$s)
      z
    }
    first <- moments(ratio)
    eq <- NULL
    if(responsive) {
      z <- engine$energy_population_cpp(bb[1],bb[2],pp[1],pp[2],ratio,ratio,
        gh$a,gh$w,gl$x,gl$w,B,sd,as.integer(nx),1L)
      if(z$status!=1L) stop("Continuous-population equilibrium did not converge; increase nquad/nxpop.")
      eq <- moments(z$ratio_equilibrium)
      eq$residual <- log((eq$E1/eq$M1)/(eq$E2/eq$M2)/eq$s)
      if(!is.finite(eq$residual) || abs(eq$residual)>2e-6)
        stop("Population equilibrium failed its independent moment residual check.")
    }
    list(first=first,responsive=eq)
  }
})
