# Refine the near-corner transitions, where small integration errors can bend
# the A/B boundary noticeably. This changes numerical integration, not contours.
source('08_boredom_model.R')
refine_boredom_density_grid<-function(compulsory=.5,N=30,n=8192L){
 tag<-paste0('c_',format(compulsory,trim=TRUE),'_N_',N)
 file<-file.path('08_boredom_results',paste0('energy_',tag,'.rds'));d<-readRDS(file)
 den<-boredom_qmc_density(n);beta<-d$beta;B<-d$B;checks<-list()
 for(j in 1:3)for(i in 1:3){
  xy<-expand.grid(x=seq_along(beta),y=seq_along(beta))
  xy<-xy[(i==3 & beta[xy$x]>=.6 & beta[xy$y]>=.8)|(j==3 & beta[xy$x]<=.2 & beta[xy$y]>=.8),,drop=FALSE]
  p<-d$panels[[i,j]]
  for(k in seq_len(nrow(xy))){
   x<-xy$x[k];y<-xy$y[k]
   fit<-boredom_population(beta[c(x,y)],c(p$phi1,p$phi2),B,compulsory,N,den,140L,1e-5,maxit=80L)
   if(!fit$responsive$converged)fit<-boredom_population(beta[c(x,y)],c(p$phi1,p$phi2),B,compulsory,N,boredom_qmc_density(32768),160L,1e-5,start=fit$responsive$moment,maxit=100L)
   if(!fit$responsive$converged)stop('Refined density cell unresolved; no result saved.')
   for(mode in c('first','responsive')){m<-fit[[mode]]$mapped;p[[mode]]$s0[x,y]<-m[4]/B;p[[mode]]$s1[x,y]<-m[1]/B;p[[mode]]$s2[x,y]<-1-(m[4]+m[1])/B}
   p$residual[x,y]<-fit$responsive$residual;p$order[x,y]<-fit$order
   checks[[length(checks)+1L]]<-c(i=i,j=j,x=x,y=y,n=fit$order,residual=fit$responsive$residual)
  }
  d$panels[[i,j]]<-p
 }
 saveRDS(d,file);saveRDS(checks,file.path('08_boredom_results','corner_precision_checks.rds'))
 invisible(d)
}
if(sys.nframe()==0L)refine_boredom_density_grid()
