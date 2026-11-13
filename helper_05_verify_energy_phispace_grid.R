d<-readRDS('05_energy_allocation_phispace_refined_inputs.rds')
stopifnot(identical(d$panel_betas,c(.1,.5,.9)),length(d$phis)==61L,
 abs(min(d$phis)-.1)<1e-12,abs(max(d$phis)-3.1)<1e-12)
n<-length(d$phis)
for(i in 1:3)for(j in 1:3){
 p<-d$panels[[i,j]];q<-d$panels[[j,i]]
 stopifnot(all(is.finite(p$first)),all(is.finite(p$equilibrium)),all(p$ratio_equilibrium>0),
 min(p$first,p$equilibrium)>=0,max(p$first,p$equilibrium)<=d$B,
 max(abs(p$first+t(q$first)-d$B))<1e-10,
 max(abs(p$equilibrium+t(q$equilibrium)-d$B))<1e-10)
 hit<-matrix(FALSE,n-1,n-1);d0<-p$first-d$B/2;d1<-p$equilibrium-d$B/2
 for(f in list(d0,d1,d1-d0,d1+d0)){
  z<-list(f[-n,-n],f[-1,-n],f[-n,-1],f[-1,-1]);lo<-Reduce(pmin,z);hi<-Reduce(pmax,z)
  hit<-hit|(lo<=1e-8 & hi>= -1e-8 & hi-lo>1e-8)
 }
 corners<-p$direct[-n,-n]&p$direct[-1,-n]&p$direct[-n,-1]&p$direct[-1,-1]
 stopifnot(all(corners[hit]))
 cat('beta1=',p$beta1,'beta2=',p$beta2,':',sum(hit),'boundary cells verified\n')
}
cat('All phi-space domain, exchange-symmetry and boundary-coverage checks passed.\n')
