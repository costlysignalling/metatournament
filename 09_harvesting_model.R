# Tournament 1 has fixed performance. Tournament 2 is fixed or trained.
# The population integrals condition on talent 2 and split talent 1 at the
# analytical specialist boundaries; no Monte Carlo samples or pixel counts.
source('08_boredom_model.R') # Gaussian and Legendre nodes; lognormal convention
harvest_core<-local({core<-NULL;function(){
 if(is.null(core)){core<<-new.env();Rcpp::sourceCpp('09_harvesting_core.cpp',env=core)}
 core
}})
harvest_engine<-function(beta2=.4,phi1=.66,phi2=1.5,B=50,both_fixed=FALSE,order=48L,inner=16L){
 gh<-boredom_gauss(order);gl<-boredom_gauss(inner,FALSE)
 harvest_core()$harvesting_engine_cpp(beta2,phi1,phi2,B,exp(sqrt(2)*boredom_slog()*gh$x),gh$w,gl$x,2*gl$w,both_fixed)
}
harvest_population<-function(beta2=.4,phi1=.66,phi2=1.5,B=50,both_fixed=FALSE,order=48L,inner=16L,start=NULL){
 core<-harvest_core();en<-harvest_engine(beta2,phi1,phi2,B,both_fixed,order,inner)
 ref<-core$harvesting_reference_cpp(en);fixed<-ref[4]/ref[3]
 moment<-function(s)core$harvesting_moments_cpp(en,exp(s))
 f<-function(s){v<-moment(s);if(any(!is.finite(v))||any(v<=0))stop('Population moments are numerically empty.');log(v[1])+log(v[4])-log(v[2])-log(v[3])-s}
 centre<-if(is.null(start))log(fixed)else log(start)
 safe<-function(s)tryCatch(f(s),error=function(e)NA_real_)
 fc<-safe(centre)
 if(!is.finite(fc)){centre<-log(fixed);fc<-f(centre)}
 if(abs(fc)<1e-9)root<-centre else{
  found<-FALSE
  for(step in c(.15,.4,1,2,4,8)){
   l<-centre-step;h<-centre+step;fl<-safe(l);fh<-safe(h)
   if(is.finite(fl)&&fl*fc<=0){br<-c(l,centre);found<-TRUE;break}
   if(is.finite(fh)&&fc*fh<=0){br<-c(centre,h);found<-TRUE;break}
  }
  if(!found)stop('Could not bracket a responsive harvesting equilibrium.')
  root<-uniroot(f,br,tol=2e-9)$root
 }
 residual<-abs(f(root));stopifnot(residual<1e-7)
 list(beta=c(NA,beta2),phi=c(phi1,phi2),B=B,both_fixed=both_fixed,
  fixed_ratio=fixed,responsive_ratio=exp(root),reference=ref,
  first=core$harvesting_moments_cpp(en,fixed),responsive=moment(root),residual=residual,order=order,inner=inner)
}
harvest_alloc<-function(p,mode,actor){
 en<-harvest_engine(p$beta[2],p$phi[1],p$phi[2],p$B,p$both_fixed,order=16L,inner=8L)
 harvest_core()$harvesting_map_cpp(en,if(mode=='first')p$fixed_ratio else p$responsive_ratio,as.matrix(actor))
}
harvest_boundaries<-function(p,mode,a2){
 en<-harvest_engine(p$beta[2],p$phi[1],p$phi[2],p$B,p$both_fixed,order=16L,inner=8L)
 harvest_core()$harvesting_boundaries_cpp(en,if(mode=='first')p$fixed_ratio else p$responsive_ratio,a2)
}
prepare_harvest_examples<-function(){
 dir.create('09_harvesting_results',showWarnings=FALSE)
 sets<-list(
  fixed=list(list(beta2=1,phi1=.66,phi2=.66,both_fixed=TRUE),list(beta2=1,phi1=1.5,phi2=1.5,both_fixed=TRUE),list(beta2=1,phi1=.66,phi2=1.5,both_fixed=TRUE),list(beta2=1,phi1=.66,phi2=1.5)),
  beta04=lapply(list(c(.66,1.5),c(1.5,1.5),c(.66,.66),c(1.5,.66)),function(phi)list(beta2=.4,phi1=phi[1],phi2=phi[2])),
  beta08=lapply(list(c(.66,1.5),c(1.5,1.5),c(.66,.66),c(1.5,.66)),function(phi)list(beta2=.8,phi1=phi[1],phi2=phi[2])))
 out<-lapply(sets,function(set)lapply(set,function(sc){message('Harvesting example ',paste(unlist(sc),collapse=','));do.call(harvest_population,sc)}))
 saveRDS(out,'09_harvesting_results/examples.rds');invisible(out)
}
