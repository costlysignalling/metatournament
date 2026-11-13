source('06C_observed_pair_prediction.R')
z <- readRDS('06B_many_results/three_tournament_population.rds')
n <- 16384L;A <- many_lognormal_population(n,3)
full <- z$checks[[as.character(n)]]$responsive
b <- rowSums(full$allocation[,1:2])
conditional <- predict_observed_pair(A[,1:2],b,z$beta[1:2],z$phi[1:2])
error <- max(abs(conditional$allocation-full$allocation[,1:2]))
ratio_error <- abs(log((conditional$scaler[1]/conditional$scaler[2])/(full$scaler[1]/full$scaler[2])))
stopifnot(conditional$converged,error<1e-4,ratio_error<1e-5)
cat('Conditional equilibrium without target splits: max allocation error=',error,
    'log-ratio error=',ratio_error,'residual=',conditional$residual,'\n')
saveRDS(list(error=error,ratio_error=ratio_error,residual=conditional$residual),
  '06B_many_results/conditional_equilibrium_check.rds')
