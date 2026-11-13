# Keep S33 unchanged; add four comparisons at the same B, c and N.
source('08_boredom_model.R')
prepare_more_boredom<-function(B=50,compulsory=.5,N=30,n=32768L){
 scenarios<-list(list(beta=c(.8,.8),phi=c(1.5,1.5)),list(beta=c(.8,.8),phi=c(2,2)),
  list(beta=c(.8,.8),phi=c(1.5,2)),list(beta=c(.8,.8),phi=c(2,2.25)))
 density<-boredom_qmc_density(n);alpha<-seq(0,2.5,length.out=251);prep<-list()
 for(k in seq_along(scenarios)){
  sc<-scenarios[[k]];p<-boredom_population(sc$beta,sc$phi,B,compulsory,N,density,160L,1e-5,maxit=60L)
  if(!p$responsive$converged)p<-boredom_population(sc$beta,sc$phi,B,compulsory,N,boredom_qmc_density(32768),200L,1e-5,start=p$responsive$moment,maxit=70L)
  if(!p$responsive$converged)stop('New example unresolved; no result published.')
  p$maps<-lapply(c('first','responsive'),function(mode)boredom_map(p,mode,alpha))
  p$boundaries<-lapply(c('first','responsive'),function(mode)boredom_boundaries(p,mode,alpha))
  prep[[k]]<-p;cat('More-boredom example',k,'complete; residual',p$responsive$residual,'\n');flush.console()
 }
 tag<-paste0('c_',format(compulsory,trim=TRUE),'_N_',N)
 result<-list(prep=prep,alpha=alpha,B=B,compulsory=compulsory,N=N)
 saveRDS(result,file.path('08_boredom_results',paste0('more_boredom_',tag,'.rds')));invisible(result)
}
if(sys.nframe()==0L)prepare_more_boredom()
