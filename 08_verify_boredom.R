# Read the supplied precision checks, or recalculate independent comparisons.
source('08_boredom_model.R')
verify_boredom<-function(recalculate=FALSE){
 root<-'08_boredom_results'
 if(!recalculate)return(lapply(setNames(c('qmc_example_checks.rds','allocation_precision_checks.rds','parameter_precision_checks.rds'),c('integration','allocation','parameters')),function(n)readRDS(file.path(root,n))))
 examples<-readRDS(file.path(root,'examples_c_0.5_N_30.rds'));integration<-allocation<-list()
 actor<-boredom_qmc_density(256)$talent
 for(k in seq_along(examples$prep)){
  high<-examples$prep[[k]]
  low<-boredom_population(high$beta,high$phi,density=boredom_qmc_density(8192),grid_half=120L,tolerance=2e-6,start=high$responsive$moment)
  integration[[k]]<-list(low=low,high=high)
  m<-high$responsive$moment;bg<-29*c(m[1],50-m[1],m[2:3]);core<-boredom_compile()
  q1<-core$boredom_map_cpp(actor,bg,high$beta[1],high$beta[2],high$phi[1],high$phi[2],50,.5,80L)
  q2<-core$boredom_map_cpp(actor,bg,high$beta[1],high$beta[2],high$phi[1],high$phi[2],50,.5,400L)
  allocation[[k]]<-max(abs(q1-q2))/50
 }
 d<-readRDS(file.path(root,'energy_c_0.5_N_30.rds'));parameters<-list()
 for(i in c(1,3))for(j in 1:3)for(xy in list(c(3,19),c(19,3),c(11,11),c(17,13))){
  # Choose the same physical beta positions when changing the base grid size.
  ref<-seq(.01,.99,length.out=21)[xy]
  x<-which.min(abs(d$beta-ref[1]));y<-which.min(abs(d$beta-ref[2]));p<-d$panels[[i,j]]
  fit<-boredom_population(d$beta[c(x,y)],c(p$phi1,p$phi2),density=boredom_qmc_density(8192),grid_half=160L,tolerance=1e-5,maxit=80L)
  if(!fit$responsive$converged)stop('Independent precision calculation unresolved.')
  low<-c(p$first$s1[x,y],p$first$s0[x,y],p$responsive$s1[x,y],p$responsive$s0[x,y]);high<-c(fit$first$mapped[c(1,4)],fit$responsive$mapped[c(1,4)])/50
  parameters[[length(parameters)+1L]]<-list(phi=c(p$phi1,p$phi2),beta=d$beta[c(x,y)],difference=high-low)
 }
 result<-list(integration=integration,allocation=allocation,parameters=parameters)
 saveRDS(integration,file.path(root,'qmc_example_checks.rds'));saveRDS(allocation,file.path(root,'allocation_precision_checks.rds'));saveRDS(parameters,file.path(root,'parameter_precision_checks.rds'))
 result
}
