# Share the validated continuous-population equilibria used in Script 06.
# The HDR limits only coexistence tests; these ratios use the full population.
file <- '06_coexistence_beta_inputs.rds'
if(!file.exists(file)) stop('Prepare the Script 06 beta-plane inputs first.')
coex <- readRDS(file)
stopifnot(coex$B>0,identical(coex$panel_vals,c(.66,1,1.5)))
ratios <- vector('list',9L);dim(ratios)<-c(3L,3L)
for(i in 1:3)for(j in 1:3) {
 ratios[[i,j]] <- coex$panels[[i,j]]$rho
 stopifnot(all(is.finite(ratios[[i,j]]) & ratios[[i,j]]>0))
}
output <- list(B=coex$B,betas=coex$axis,phis=coex$panel_vals,ratios=ratios,
 talent=list(meanlog=0,sd_original=.5,median=1),
 source=list(file=file,method='deterministic switch-split continuous-population equilibrium',
             md5=unname(tools::md5sum(file))))
saveRDS(output,'05_equal_talent_geometry_inputs.rds')
message('Equal-talent beta inputs synchronized with Script 06.')
