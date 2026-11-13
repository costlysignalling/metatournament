source('06B_predict_many_tournaments.R')
path <- '06B_many_results/three_tournament_population.rds'
z <- readRDS(path)
A <- many_lognormal_population(65536L,3)
first <- many_allocate(A,z$beta,z$phi,z$B,z$fixed$scaler)
stopifnot(all(first$within_tolerance))
eq <- many_responsive(A,z$beta,z$phi,z$B,start=z$checks[['16384']]$responsive$scaler)
stopifnot(eq$converged)
z$checks[['65536']] <- list(first=many_means(A,first$allocation,z$beta,z$phi),responsive=eq)
saveRDS(z,path)
write.csv(do.call(rbind,lapply(names(z$checks),function(n){
  v <- z$checks[[n]]
  data.frame(population_nodes=as.integer(n),tournament=1:3,first_E=v$first$E,
    responsive_E=v$responsive$E,consistency_residual=v$responsive$residual)
})), '06B_many_results/population_convergence.csv',row.names=FALSE)
write.csv(data.frame(tournament=1:3,beta=z$beta,phi=z$phi,reference_E=z$fixed$E,
  reference_M=z$fixed$M,first_E=z$checks[['65536']]$first$E,responsive_E=eq$E,
  responsive_M=eq$M,responsive_scaler=eq$scaler),
  '06B_many_results/three_tournament_population.csv',row.names=FALSE)
print(z$checks[['65536']]$first$E);print(eq[c('E','scaler','residual')])
cat('Change from 16384 nodes:',eq$E-z$checks[['16384']]$responsive$E,'\n')
