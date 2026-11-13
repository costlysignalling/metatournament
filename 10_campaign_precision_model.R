# Accurate recipient expectations first; smooth their numerical representation
# in campaign pressure, then integrate individual optimal campaign decisions.
# No smoothing is applied across tournament parameters or to the final maps.
source('10_campaign_efficiency_model.R')
source('10_campaign_linear_exact_integration.R')
campaign_precision_compile<-local({core<-NULL;function(){if(is.null(core)){
 core<<-new.env();Rcpp::sourceCpp('10_campaign_precision_core.cpp',env=core,
  cacheDir=getOption('metatournament.precision.cpp.cache',file.path(tempdir(),'10_campaign_precision_cpp')))};core}})
.campaign_precision_environment<-new.env(parent=environment(convincing_adaptive_moments))
.campaign_precision_environment$convincing_compile<-campaign_precision_compile
.campaign_precision_environment$convincing_adaptive_moments<-convincing_adaptive_moments
environment(.campaign_precision_environment$convincing_adaptive_moments)<-.campaign_precision_environment
campaign_precision_baseline<-convincing_adaptive_baseline
environment(campaign_precision_baseline)<-.campaign_precision_environment
campaign_precision_knots<-function(refine=FALSE){
 x<-sort(unique(c(seq(-2.5,2.5,by=.25),c(-10,-7.5,-5,-3.5,3.5,5,7.5,10))))
 if(refine)x<-sort(c(x,(head(x,-1)+tail(x,-1))/2));x
}
campaign_precision_profile<-function(beta,phi,N=30,B=50,refine=FALSE,tolerance=1e-5){
 core<-campaign_precision_compile();den<-convincing_density(2048L)
 linear<-N>1e30&&all(phi==1)
 if(linear){
  balance<-function(z){m<-campaign_linear_recipient_moments(beta,B,log_ratio=z)
   z-log(m[1]/m[2]*m[3]/m[4])}
  ratio<-if(beta[1]==beta[2])0 else uniroot(balance,c(-5,5),tol=1e-9)$root
  moment<-campaign_linear_recipient_moments(beta,B,log_ratio=ratio)[1:3]
 }else{
  ref<-campaign_plane_reference(beta,phi,B)
  q<-core$convincing_receivers_cpp(den$talent,c(B/2,1,ref['Lratio']),0,beta[1],beta[2],phi[1],phi[2],B,1e60,32L)[,1]
  initial<-c(sum(den$weights*q),sum(den$weights*boredom_w(q,den$talent[,1],beta[1],phi[1])),
   sum(den$weights*boredom_w(B-q,den$talent[,2],beta[2],phi[2])))
  eq<-core$convincing_reference_density_cpp(den$talent,den$weights,initial,
   beta[1],beta[2],phi[1],phi[2],B,.5,N,40L,250L,2e-5,TRUE)
  if(N>1e30){
   evaluate<-function(z){a<-core$convincing_adaptive_moments_cpp(c(B/2,1,exp(z)),0,
    beta[1],beta[2],phi[1],phi[2],B,N,boredom_slog(),80L,tolerance)$moments[1,]
    list(residual=z-log(a[1]/a[2]*a[3]/a[4]),moment=a[1:3])}
   if(all(beta==beta[1])&&all(phi==phi[1])){ans<-evaluate(0);moment<-ans$moment;moment[1]<-B/2;moment[2:3]<-mean(moment[2:3])
   }else{
    z<-log(eq$moment[1]/eq$moment[2]*eq$moment[3]/(B-eq$moment[1]));ans<-evaluate(z)
    for(it in 1:8){if(abs(ans$residual)<3e-6)break
     slope<-(evaluate(z+.01)$residual-ans$residual)/.01
     z<-z-max(-.5,min(.5,ans$residual/slope));ans<-evaluate(z)}
    if(abs(ans$residual)>3e-5)stop('Accurate global baseline did not converge: ',ans$residual)
    moment<-ans$moment
   }
  }else{
   p<-list(beta=beta,phi=phi,B=B,N=N,baseline=list(responsive=list(moment=eq$moment)))
   moment<-campaign_precision_baseline(p,tolerance=2e-6)$responsive$moment
  }
 }
 pressure<-campaign_precision_knots(refine)
 if(linear){moments<-t(vapply(pressure,function(v)campaign_linear_recipient_moments(beta,B,v,ratio),numeric(4)))
  error<-matrix(0,length(pressure),4)
 }else{
  symmetric<-all(beta==beta[1])&&all(phi==phi[1])
  wanted<-if(symmetric)which(pressure>=0)else seq_along(pressure)
  quad<-core$convincing_adaptive_moments_cpp(moment,pressure[wanted],beta[1],beta[2],phi[1],phi[2],B,N,boredom_slog(),80L,tolerance)
  moments<-error<-matrix(NA_real_,length(pressure),4)
  moments[wanted,]<-quad$moments;error[wanted,]<-quad$error_estimate
  if(symmetric){negative<-which(pressure<0);reflected<-match(-pressure[negative],pressure)
   moments[negative,]<-moments[reflected,c(4,3,2,1)];error[negative,]<-error[reflected,c(4,3,2,1)]}
 }
 stopifnot(all(is.finite(moments)),all(moments>0))
 list(beta=beta,phi=phi,N=N,B=B,pressure=pressure,moments=moments,error=error,baseline=moment,
  method='adaptive recipient integration; cubic pressure interpolation',version=1L,tolerance=tolerance)
}
campaign_precision_evaluate<-function(profile,efficiencies=c(.05,.1,.2),n=2048L,npressure=161L,nx=32L){
 p<-profile
 if(p$N>1e30&&all(p$phi==1)){
  values<-lapply(efficiencies,function(e)campaign_linear_global_mean(p$beta,p$B,efficiency=e,npressure=npressure,tolerance=1e-5)$result[1:5])
  names(values)<-sprintf('%.2f',efficiencies);return(values)
 }
 den<-convincing_density(n);core<-campaign_precision_compile();B<-p$B
 values<-cbind(log(p$moments[,1]/p$moments[,4]),log(p$moments[,2]),log(p$moments[,3]))
 # Hyman-filtered cubic interpolation keeps the exact monotone directions.
 expected<-c(1,1,-1);fits<-vector('list',3)
 for(k in 1:3){y<-expected[k]*values[,k];correction<-max(abs(isoreg(y)$yf-y))
  if(correction>2e-3)stop('Recipient moments are not sufficiently monotone: ',correction)
  fits[[k]]<-splinefun(p$pressure,expected[k]*cummax(isoreg(y)$yf),method='hyman')
 }
 out<-list()
 for(e in efficiencies){effort<-seq(-B,B,length.out=npressure);delta<-e*effort
  eta<-fits[[1]](delta);moments<-cbind(B*plogis(eta),exp(fits[[2]](delta)),exp(fits[[3]](delta)),B*plogis(-eta))
  o<-core$convincing_focal_cpp(den$talent,moments,effort,p$beta[1],p$beta[2],p$phi[1],p$phi[2],B,p$N,nx,TRUE,FALSE)
  out[[sprintf('%.2f',e)]]<-c(campaign=sum(den$weights*abs(o[,3]))/B,
   campaign1=sum(den$weights*pmax(o[,3],0))/B,campaign2=sum(den$weights*pmax(-o[,3],0))/B,
   training1=sum(den$weights*o[,1])/B,training2=sum(den$weights*o[,2])/B)
 }
 out
}
