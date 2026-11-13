# Denser talent integration at a modest set of pressure anchors, followed by
# monotone cubic interpolation. Adaptive quadrature supplies the rare tails.
source('10_campaign_precision_model.R')
campaign_dense_profile<-function(beta,phi,N=30,B=50,n=32768L,refine=FALSE,tolerance=1e-5){
 core<-campaign_precision_compile();den<-convincing_density(n)
 pressure<-campaign_precision_knots(refine);ref<-campaign_plane_reference(beta,phi,B)
 linear<-N>1e30&&all(phi==1)
 if(linear){
  balance<-function(z){m<-campaign_linear_recipient_moments(beta,B,log_ratio=z)
   z-log(m[1]/m[2]*m[3]/m[4])}
  ratio<-if(beta[1]==beta[2])0 else uniroot(balance,c(-5,5),tol=1e-9)$root
  moment<-campaign_linear_recipient_moments(beta,B,log_ratio=ratio)[1:3]
  moments<-t(vapply(pressure,function(v)campaign_linear_recipient_moments(beta,B,v,ratio),numeric(4)))
  error<-matrix(0,length(pressure),4);baseline_residual<-0
 }else{
  q0<-core$convincing_receivers_cpp(den$talent,c(B/2,1,ref['Lratio']),0,
   beta[1],beta[2],phi[1],phi[2],B,1e60,40L)[,1]
  initial<-c(sum(den$weights*q0),sum(den$weights*boredom_w(q0,den$talent[,1],beta[1],phi[1])),
   sum(den$weights*boredom_w(B-q0,den$talent[,2],beta[2],phi[2])))
  eq<-core$convincing_reference_density_cpp(den$talent,den$weights,initial,
   beta[1],beta[2],phi[1],phi[2],B,.5,N,40L,80L,2e-5,TRUE)
  if(!eq$converged){
   p<-list(beta=beta,phi=phi,B=B,N=N,baseline=list(responsive=list(moment=eq$moment)))
   base<-campaign_precision_baseline(p,tolerance=2e-6)
   moment<-base$responsive$moment;baseline_residual<-base$responsive$residual
  }else{moment<-eq$moment;baseline_residual<-eq$residual}
  q<-core$convincing_receivers_cpp(den$talent,moment,pressure,
   beta[1],beta[2],phi[1],phi[2],B,N,80L)
  moments<-t(vapply(seq_along(pressure),function(k){x<-q[,k];c(sum(den$weights*x),
   sum(den$weights*boredom_w(x,den$talent[,1],beta[1],phi[1])),
   sum(den$weights*boredom_w(B-x,den$talent[,2],beta[2],phi[2])),sum(den$weights*(B-x)))},numeric(4)))
  error<-matrix(NA_real_,length(pressure),4)
  tails<-which(abs(pressure)>2.5|apply(moments<=0,1,any)|pmin(moments[,1],moments[,4])<B*1e-3)
  if(length(tails)){
   quad<-core$convincing_adaptive_moments_cpp(moment,pressure[tails],
    beta[1],beta[2],phi[1],phi[2],B,N,boredom_slog(),80L,tolerance)
   moments[tails,]<-quad$moments;error[tails,]<-quad$error_estimate
  }
 }
 if(all(beta==beta[1])&&all(phi==phi[1])){
  reflected<-match(-pressure,pressure)
  moments<-(moments+moments[reflected,c(4,3,2,1)])/2
 }
 stopifnot(all(is.finite(moments)),all(moments>0),baseline_residual<3e-5)
 list(beta=beta,phi=phi,N=N,B=B,pressure=pressure,moments=moments,error=error,baseline=moment,
  baseline_residual=baseline_residual,method='dense stratified recipient integration with adaptive rare tails; monotone cubic pressure interpolation',
  version=2L,tolerance=tolerance,recipient_order=n)
}
