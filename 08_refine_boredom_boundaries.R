# Directly solve equality of first-response and responsive productive imbalance.
# Only geometry is refined here; colours retain the integrated population means.
source('08_boredom_model.R')
refine_boredom_boundaries<-function(n=32768L){
 d<-readRDS('08_boredom_results/energy_c_0.5_N_30.rds');beta<-d$beta;den<-boredom_qmc_density(n);out<-list()
 bilinear<-function(z,x,y){i<-pmax(1L,pmin(length(beta)-1L,findInterval(x,beta)));j<-pmax(1L,pmin(length(beta)-1L,findInterval(y,beta)));tx<-(x-beta[i])/(beta[i+1]-beta[i]);ty<-(y-beta[j])/(beta[j+1]-beta[j]);(1-tx)*(1-ty)*z[cbind(i,j)]+tx*(1-ty)*z[cbind(i+1,j)]+(1-tx)*ty*z[cbind(i,j+1)]+tx*ty*z[cbind(i+1,j+1)]}
 for(spec in list(list(i=1,j=3,along='y',at=seq(.75,.99,length.out=9),range=c(.01,.45)),list(i=3,j=2,along='x',at=seq(.6,.99,length.out=12),range=c(.8,.99)))){
  p<-d$panels[[spec$i,spec$j]];field<-(p$responsive$s1-p$responsive$s2)-(p$first$s1-p$first$s2);warm<-NULL;records<-list()
  for(a in spec$at){
   coords<-function(v)if(spec$along=='y')c(v,a)else c(a,v)
   old<-function(v){xy<-coords(v);bilinear(field,xy[1],xy[2])}
   ab<-spec$range;fa<-old(ab[1]);fb<-old(ab[2])
   guess<-if(fa*fb<=0)uniroot(old,ab,tol=1e-7)$root else mean(ab)
   value<-function(v){
    fit<-boredom_population(coords(v),c(p$phi1,p$phi2),density=den,grid_half=60L,tolerance=1e-5,start=warm,maxit=70L)
    if(!fit$responsive$converged)fit<-boredom_population(coords(v),c(p$phi1,p$phi2),density=boredom_qmc_density(131072),grid_half=80L,tolerance=1e-5,start=fit$responsive$moment,maxit=120L)
    if(!fit$responsive$converged)stop('Boundary population did not converge.')
    warm<<-fit$responsive$moment
    (2*(fit$responsive$mapped[1]-fit$first$mapped[1])+fit$responsive$mapped[4]-fit$first$mapped[4])/50
   }
   # A bracketed solve is stable even when the transition is shallow.
   lo<-max(ab[1],guess-.035);hi<-min(ab[2],guess+.035)
   flo<-value(lo);fhi<-value(hi)
   if(flo*fhi>0){lo<-ab[1];hi<-ab[2];flo<-value(lo);fhi<-value(hi)}
   if(flo*fhi>0){
    x<-if(flo>0)ab[1]-.01 else ab[2]+.01;fx<-NA_real_
   }else{
   root<-uniroot(value,c(lo,hi),f.lower=flo,f.upper=fhi,tol=.0005,maxiter=30)
   x<-root$root;fx<-root$f.root
   }
   records[[length(records)+1L]]<-c(at=a,root=x,residual=fx,original=guess)
   cat('Boundary',spec$i,spec$j,'at',a,'root',x,'residual',fx,'\n');flush.console()
  }
  rr<-do.call(rbind,records)
  inside<-is.finite(rr[,'residual'])
  for(k in which(diff(inside)!=0)){
   offrow<-if(!inside[k])k else k+1L
   edge<-if(rr[offrow,'root']<spec$range[1])spec$range[1]else spec$range[2]
   edge_value<-function(anchor){
    pair<-if(spec$along=='y')c(edge,anchor)else c(anchor,edge)
    fit<-boredom_population(pair,c(p$phi1,p$phi2),density=den,grid_half=80L,tolerance=1e-5,start=warm,maxit=80L)
    if(!fit$responsive$converged)stop('Panel-edge population unresolved.')
    warm<<-fit$responsive$moment
    (2*(fit$responsive$mapped[1]-fit$first$mapped[1])+fit$responsive$mapped[4]-fit$first$mapped[4])/50
   }
   aa<-rr[k:(k+1),'at'];ee<-vapply(aa,edge_value,0.)
   if(prod(ee)<=0){cross<-uniroot(edge_value,aa,f.lower=ee[1],f.upper=ee[2],tol=.0005);rr<-rbind(rr,c(at=cross$root,root=edge,residual=cross$f.root,original=NA_real_))}
  }
  spec$records<-rr[order(rr[,'at']),,drop=FALSE];spec$n<-n;out[[paste(spec$i,spec$j,sep='_')]]<-spec
  saveRDS(out,'08_boredom_results/boundary_roots_c_0.5_N_30.rds')
 }
 saveRDS(out,'08_boredom_results/boundary_roots_c_0.5_N_30.rds');invisible(out)
}
if(sys.nframe()==0L)refine_boredom_boundaries()
