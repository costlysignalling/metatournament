source('06B_predict_many_tournaments.R')
dir.create('06B_many_results',showWarnings=FALSE)
beta <- c(.9,.6,.6);phi <- c(.8,1.5,1);B <- 50
fixed <- many_random_reference(beta,phi,B)
records <- list()
for(n in c(4096L,16384L)) {
  started <- proc.time()[3]
  A <- many_lognormal_population(n,3)
  first <- many_allocate(A,beta,phi,B,fixed$scaler)
  stopifnot(all(first$within_tolerance))
  initial <- many_means(A,first$allocation,beta,phi)
  eq <- many_responsive(A,beta,phi,B,start=fixed$scaler)
  stopifnot(eq$converged)
  records[[as.character(n)]] <- list(first=initial,responsive=eq)
  cat('n=',n,'seconds=',proc.time()[3]-started,'first E=',initial$E,
      'responsive E=',eq$E,'residual=',eq$residual,'\n')
}
A <- many_lognormal_population(16384L,3)
eq <- records[['16384']]$responsive
saveRDS(list(beta=beta,phi=phi,B=B,fixed=fixed,checks=records),
        '06B_many_results/three_tournament_population.rds')
write.csv(data.frame(tournament=1:3,beta,phi,reference_E=fixed$E,reference_M=fixed$M,
  first_E=records[['16384']]$first$E,responsive_E=eq$E,responsive_M=eq$M,
  responsive_scaler=eq$scaler),'06B_many_results/three_tournament_population.csv',row.names=FALSE)

# Pair restriction: recover the full optimum's 1:2 split conditional on its
# total budget for those two tournaments, without using alpha_3.
take <- seq_len(1024)
pairB <- rowSums(eq$allocation[take,1:2,drop=FALSE])
# Use the exact scalers at which the saved allocation was optimized.
# Population-implied scalers can still differ by the equilibrium residual.
pair <- many_allocate(A[take,1:2],beta[1:2],phi[1:2],pairB,eq$response_scaler[1:2])
pair_error <- max(abs(pair$allocation-eq$allocation[take,1:2]))
stopifnot(pair_error<2e-6)
cat('Conditional pair recovery max absolute error:',pair_error,'\n')

# Analytic symmetric phi=1 solution, T=2,...,10.
linear_error <- 0
for(T in 2:10) {
  a <- exp(seq(-.7,.7,length.out=T));b <- .6
  z <- many_allocate(a,rep(b,T),rep(1,T),B)
  v <- a^(1/(1-b));o <- order(v,decreasing=TRUE)
  k <- max(which((B+seq_len(T))*v[o]>cumsum(v[o])))
  expected <- numeric(T);expected[o[seq_len(k)]] <- (B+k)*v[o[seq_len(k)]]/sum(v[o[seq_len(k)]])-1
  linear_error <- max(linear_error,max(abs(z$allocation-expected)))
}
stopifnot(linear_error<2e-6)

# Nonconcave DP bounds must contain independently finer-grid objective values.
nonconcave_checks <- list()
for(a in list(c(.457,.457,.457),c(2,.3,1.5),c(.8,1.3,.4))) {
  small <- many_allocate(a,rep(.5,3),rep(3,3),1,tolerance=1,max_grid=128L)
  fine <- many_allocate(a,rep(.5,3),rep(3,3),1,tolerance=1e-7,max_grid=1024L)
  stopifnot(fine$utility>=small$utility-1e-10,
    fine$utility<=small$utility+small$utility_gap_bound+1e-10)
  nonconcave_checks[[length(nonconcave_checks)+1]] <- list(alpha=a,coarse=small,fine=fine)
}

# Three-tournament parameter checks: first-response trends and feedback.
# Screening is labelled as such; bounds and convergence are retained per case.
scenarios <- list(
  linear_difficulty=list(beta=c(.3,.5,.7),phi=c(1,1,1)),
  equal_difficulty=list(beta=rep(.4,3),phi=c(.6,1,1.4)),
  extreme_returns=list(beta=rep(.4,3),phi=c(.6,1,4)),
  mixed_example=list(beta=c(.6,.4,.4),phi=c(.66,1,1.2)),
  harder_extreme=list(beta=c(.3,.5,.7),phi=c(2,.8,.8)),
  easier_extreme=list(beta=c(.3,.5,.7),phi=c(.8,.8,2)))
screen <- list();P <- many_lognormal_population(2048,3)
for(name in names(scenarios)) {
  s <- scenarios[[name]];ref <- many_random_reference(s$beta,s$phi,B)
  first <- many_allocate(P,s$beta,s$phi,B,ref$scaler,tolerance=1e-3,max_grid=1024L)
  m <- many_means(P,first$allocation,s$beta,s$phi)
  fit <- tryCatch(many_responsive(P,s$beta,s$phi,B,start=ref$scaler,
    allocation_tolerance=1e-3,max_grid=1024L),error=function(e)list(converged=FALSE,error=conditionMessage(e)))
  screen[[name]] <- list(parameters=s,first=m,first_gap=max(first$relative_gap_bound),responsive=fit)
  cat(name,': first E=',m$E,'equilibrium=',if(isTRUE(fit$converged))fit$E else 'not certified','\n')
}
saveRDS(list(pair_error=pair_error,linear_error=linear_error,nonconcave=nonconcave_checks,
  screening=screen),'06B_many_results/verification.rds')
write.csv(do.call(rbind,lapply(names(screen),function(name){
  x <- screen[[name]];data.frame(example=name,tournament=1:3,
    beta=x$parameters$beta,phi=x$parameters$phi,first_E=x$first$E,
    responsive_E=if(isTRUE(x$responsive$converged))x$responsive$E else NA_real_,
    first_global_gap_bound=x$first_gap,
    responsive_residual=if(is.null(x$responsive$residual))NA_real_ else x$responsive$residual)
})),'06B_many_results/parameter_checks.csv',row.names=FALSE)
cat('All completed verification checks passed.\n')
