# Rebuild the expected-density examples and the beta-space figure data.
# Run from the scripts folder. Redraws do not require compilation.
source('08_boredom_model.R')
source('08_prepare_boredom_density_grid.R')
source('08_refine_boredom_density_grid.R')
source('08_refine_boredom_boundaries.R')
source('08_prepare_more_boredom.R')
prepare_boredom_figures<-function(B=50,compulsory=.5,N=30,nexample=32768L,nbeta=21L,nparameter=2048L){
 dir.create('08_boredom_results',showWarnings=FALSE)
 density<-boredom_qmc_density(nexample);boredom_compile()
 tag<-paste0('c_',format(compulsory,trim=TRUE),'_N_',N)
 examples<-list(list(beta=c(.8,.8),phi=c(1,1)),list(beta=c(.8,.8),phi=c(1.25,1.25)),
  list(beta=c(.8,.6),phi=c(1.5,1)),list(beta=c(.6,.4),phi=c(.66,1)))
 prep<-lapply(examples,function(sc){
  p<-boredom_population(sc$beta,sc$phi,B,compulsory,N,density,160L,2e-6)
  if(!p$responsive$converged)stop('Example did not converge; increase integration resolution.')
  p
 })
 alpha<-seq(0,2.5,length.out=251)
 for(k in seq_along(prep)){
  prep[[k]]$maps<-lapply(c('first','responsive'),function(mode)boredom_map(prep[[k]],mode,alpha))
  prep[[k]]$boundaries<-lapply(c('first','responsive'),function(mode)boredom_boundaries(prep[[k]],mode,alpha))
 }
 saveRDS(list(prep=prep,alpha=alpha,B=B,compulsory=compulsory,N=N),file.path('08_boredom_results',paste0('examples_',tag,'.rds')))
 prepare_boredom_density_grid(B=B,compulsory=compulsory,N=N,nbeta=nbeta,n=nparameter)
 refine_boredom_density_grid(compulsory,N)
 if(B==50 && compulsory==.5 && N==30)refine_boredom_boundaries()
 prepare_more_boredom(B,compulsory,N)
 invisible(prep)
}
if(sys.nframe()==0L)prepare_boredom_figures()
