# General asymmetric tournaments: first response and responsive population.
# Run from this folder. Only Rcpp (and Rtools on Windows) is needed.
# The C++ helper contains the allocation calculation; population feedback,
# validation, input controls and plotting are kept here in readable R.
#
# Quick start in R:
# source('06B_predict_many_tournaments.R')
# example <- many_example()      # prints allocations and draws two pie charts
#
# IMPORTANT: A responsive result is a numerically checked self-consistent
# solution, not a guarantee of equilibrium existence, uniqueness, or convergence
# of repeated simultaneous best responses. Increase the population sample to
# check integration accuracy. Nonconcave actors have reported utility bounds.

many_compile <- function() {
  if(!exists('many_allocations_cpp',mode='function')) {
    if(!requireNamespace('Rcpp',quietly=TRUE)) stop('Install Rcpp first.')
    # Compile once per R session in R's local temporary folder. Reusing a DLL
    # from a synchronized project cache can fail Windows application control.
    # Rebuild from the supplied source; no security settings are changed.
    Rcpp::sourceCpp('engine/06B_many_tournaments_core.cpp',env=environment(many_compile),
      cacheDir=file.path(tempdir(),'06B_many_cpp_cache'),rebuild=TRUE)
  }
}

many_allocate <- function(alpha,beta,phi,B=50,scaler=rep(1,length(beta)),
                          tolerance=1e-4,max_grid=2048L) {
  many_compile()
  if(is.null(dim(alpha))) alpha <- matrix(alpha,nrow=1)
  alpha <- as.matrix(alpha); T <- ncol(alpha)
  if(length(beta)!=T || length(phi)!=T || length(scaler)!=T ||
     any(!is.finite(c(alpha,beta,phi,B,scaler))) || any(alpha<0) ||
     any(beta<=0 | beta>=1) || any(phi<=0) || any(scaler<=0) || any(B<0) ||
     !(length(B)%in%c(1L,nrow(alpha))) || !is.finite(tolerance) || tolerance<=0 ||
     length(max_grid)!=1L || max_grid<128 || max_grid!=as.integer(max_grid))
    stop('Use nonnegative talents/budgets, 0<beta<1, phi>0, positive scalers, and matching dimensions.')
  # A common multiplier changes utility units, but not the allocation.
  scaler <- scaler/exp(mean(log(scaler)))
  ans <- many_allocations_cpp(alpha,rep(B,length.out=nrow(alpha)),beta,phi,
    scaler,128L,as.integer(max_grid),tolerance)
  if(any(!is.finite(c(ans$allocation,ans$utility,ans$utility_gap_bound))))
    stop('The calculation exceeds floating-point range; rescale inputs or reduce extreme parameters.')
  colnames(ans$allocation) <- paste0('tournament_',seq_len(T))
  ans$relative_gap_bound <- ans$utility_gap_bound/(1+abs(ans$utility))
  ans$within_tolerance <- ans$relative_gap_bound<=tolerance
  ans$scaler <- scaler
  ans
}

many_weights <- function(n,weights=NULL) {
  if(is.null(weights)) weights <- rep(1,n)
  if(length(weights)!=n || any(!is.finite(weights)) || any(weights<0) || sum(weights)<=0)
    stop('Population weights must be nonnegative and have positive total.')
  weights/sum(weights)
}

many_means <- function(alpha,X,beta,phi,weights=NULL) {
  alpha <- as.matrix(alpha); X <- as.matrix(X)
  if(!identical(dim(alpha),dim(X))) stop('Talents and allocations must have the same dimensions.')
  weights <- many_weights(nrow(alpha),weights)
  W <- X
  for(t in seq_along(beta)) W[,t] <- expm1(phi[t]*log1p(alpha[,t]*expm1(beta[t]*log1p(X[,t]))))
  E <- colSums(X*weights); M <- colSums(W*weights)
  list(E=E,M=M,scaler=E/M)
}

# A deterministic low-discrepancy population. Full JOINT distributions matter:
# Sigma is the covariance matrix of log talents, not of original-scale talents.
# For other distributions, supply your own alpha matrix and optional weights.
many_halton <- function(n,d,skip=31L) {
  primes <- c(2,3,5,7,11,13,17,19,23,29,31,37,41,43,47,53,59,61,67,71)
  if(d>length(primes)) stop('Supply a population matrix for more than 20 dimensions.')
  ans <- matrix(0,n,d)
  for(j in seq_len(d)) {
    index <- seq_len(n)+skip; factor <- 1/primes[j]
    while(any(index>0)) {
      ans[,j] <- ans[,j]+factor*(index%%primes[j]);index <- index%/%primes[j]
      factor <- factor/primes[j]
    }
  }
  ans
}
many_lognormal_population <- function(n=4096,T=3,meanlog=rep(0,T),
  Sigma=diag(log((1+sqrt(2))/2),T)) {
  if(length(meanlog)!=T || !all(dim(Sigma)==c(T,T)))stop('Check meanlog and Sigma dimensions.')
  Z <- qnorm(many_halton(n,T))%*%chol(Sigma)
  exp(sweep(Z,2,meanlog,'+'))
}

# Fixed reference: independent talents and a UNIFORM SIMPLEX allocation,
# Dirichlet(1,...,1), independent of talents. Each marginal allocation is
# B*Beta(1,T-1), not Uniform(0,B) except when T=2. E_t=B/T; M_t differs.
# Supply many_means(alpha,reference_X,...) instead for a measured reference,
# correlated talent-budget data, or a different initial allocation rule.
many_random_reference <- function(beta,phi,B=50,meanlog=rep(0,length(beta)),
  sdlog=rep(sqrt(log((1+sqrt(2))/2)),length(beta)),nquad=64L) {
  T <- length(beta)
  if(length(phi)!=T || length(meanlog)!=T || length(sdlog)!=T || B<=0 ||
     any(beta<=0 | beta>=1) || any(phi<=0) || any(sdlog<0))stop('Invalid reference parameters.')
  J <- matrix(0,nquad,nquad)
  J[cbind(seq_len(nquad-1),2:nquad)] <- sqrt(seq_len(nquad-1))
  gh <- eigen(J+t(J),symmetric=TRUE)
  weights <- gh$vectors[1,]^2
  M <- numeric(T)
  for(t in seq_len(T)) {
    talent <- exp(meanlog[t]+sdlog[t]*gh$values)
    f <- function(u) {
      x <- if(T==1) rep(B,length(u)) else B*(-expm1(log1p(-u)/(T-1)))
      as.vector(crossprod(weights,expm1(phi[t]*log1p(outer(talent,expm1(beta[t]*log1p(x)))))))
    }
    M[t] <- if(T==1)f(.5) else integrate(f,0,1,rel.tol=1e-9,subdivisions=500)$value
  }
  E <- rep(B/T,T)
  list(E=E,M=M,scaler=E/M)
}

many_responsive <- function(alpha,beta,phi,B=50,start=NULL,weights=NULL,
  residual_tolerance=2e-7,max_iterations=35L,allocation_tolerance=1e-4,
  max_grid=2048L) {
  alpha <- as.matrix(alpha);T <- ncol(alpha);weights <- many_weights(nrow(alpha),weights)
  if(is.null(start))start <- rep(1,T)
  if(length(start)!=T || any(start<=0) || any(!is.finite(start)))stop('Invalid starting scalers.')
  z <- log(start[-T]/start[T])
  evaluate <- function(z) {
    scaler <- exp(c(z,0))
    a <- many_allocate(alpha,beta,phi,B,scaler,allocation_tolerance,max_grid)
    if(!all(a$within_tolerance))stop('Some global allocation bounds exceed tolerance; increase max_grid or explicitly relax allocation_tolerance.')
    m <- many_means(alpha,a$allocation,beta,phi,weights)
    if(any(!is.finite(m$scaler)) || any(m$E<=0) || any(m$M<=0))
      stop('An inactive tournament has E=M=0; an entry analysis is needed before claiming equilibrium.')
    target <- log(m$scaler[-T]/m$scaler[T])
    list(residual=target-z,means=m,allocation=a)
  }
  state <- evaluate(z);history <- numeric()
  for(iter in seq_len(max_iterations)) {
    err <- if(T==1)0 else max(abs(state$residual));history <- c(history,err)
    if(err<residual_tolerance)break
    # Small differences can be dominated by strategy switches in a finite
    # population. Keep broad differences until close to a root, then resolve
    # the smooth local response precisely. This does not relax convergence.
    step <- if(err>max(1e-5,50*residual_tolerance)) .02 else 1e-6
    J <- sapply(seq_len(T-1L),function(j) {
      zz <- z;zz[j] <- zz[j]+step
      (evaluate(zz)$residual-state$residual)/step
    })
    J <- matrix(J,T-1,T-1)
    direction <- tryCatch(as.vector(solve(J,-state$residual)),error=function(e)state$residual)
    if(max(abs(direction))>2)direction <- direction*2/max(abs(direction))
    improved <- FALSE
    # A damped consistency update is a fallback when the estimated Newton
    # direction is unhelpful near a change of optimal strategy.
    for(proposal in list(direction,state$residual)) {
      if(max(abs(proposal))>2)proposal <- proposal*2/max(abs(proposal))
      for(damping in 2^-(0:10)) {
        trial <- evaluate(z+damping*proposal)
        if(sum(trial$residual^2)<sum(state$residual^2)) {
          z <- z+damping*proposal;state <- trial;improved <- TRUE;break
        }
      }
      if(improved)break
    }
    if(!improved)break
  }
  residual <- if(T==1)0 else max(abs(state$residual))
  result <- c(state$means,list(allocation=state$allocation$allocation,
    response_scaler=exp(c(z,0)),residual=residual,converged=residual<residual_tolerance,
    max_allocation_gap=max(state$allocation$relative_gap_bound),history=history,
    n=nrow(alpha),weights=weights))
  if(!result$converged)warning('Equilibrium residual did not reach tolerance; do not label this an equilibrium.')
  result
}

many_pie <- function(result,title='Predicted allocation',colours=NULL) {
  x <- as.numeric(result$allocation[1,]);T <- length(x)
  if(is.null(colours))colours <- if(T==3)c('#EE0066','#0088FF','#40A56A') else hcl.colors(T,'Dark 3')
  if(sum(x)==0){plot.new();title(main=paste(title,'(zero budget)'));return(invisible(result))}
  # Compact three-line labels avoid clipping in a two-panel plotting window.
  labels <- ifelse(x>0,paste0('T',seq_len(T),'\n',formatC(x,digits=2,format='f'),
    '\n',formatC(100*x/sum(x),digits=1,format='f'),'%'),'')
  pie(x,labels=labels,col=colours,border='white',main=title,clockwise=TRUE,
    radius=.72,cex=.8)
  invisible(result)
}

many_example <- function(beta=c(.6,.4,.4), phi=c(.66,1,1.2), B=50,
                         actor=rep(1,length(beta)), n_population=4096L,
                         draw=TRUE) {
  # Change the arguments to explore other tournaments and focal actors.
  # Baseline talent marginals have median = 1 and SD = 0.5.
  population <- many_lognormal_population(n_population,T=length(beta))
  fixed <- many_random_reference(beta,phi,B)
  responsive <- many_responsive(population,beta,phi,B,start=fixed$scaler)
  if(!responsive$converged)stop(sprintf(
    'Responsive solve did not converge (residual %.3g). Try increasing n_population; inspect many_responsive() for diagnostics.',
    responsive$residual))
  first <- many_allocate(actor,beta,phi,B,fixed$scaler)
  equilibrium <- many_allocate(actor,beta,phi,B,responsive$scaler)
  table <- rbind(first$allocation,equilibrium$allocation)
  rownames(table) <- c('Fixed equal E_t reference','Responsive E_t')
  print(round(table,3))
  if(draw){old <- par(mfrow=c(1,2));on.exit(par(old));many_pie(first,expression(Fixed~equal~E[t]~reference));many_pie(equilibrium,expression(Responsive~E[t]))}
  invisible(list(beta=beta,phi=phi,B=B,population=population,fixed=fixed,
    responsive=responsive,actor=actor,first=first,equilibrium=equilibrium))
}

if(sys.nframe()==0L)many_example()


many_example(beta=c(.6,.8,.4,0.8),phi=c(.66,1,1.2,1.4),B=50,actor=c(1,1,1,1.5))
