# Run from the scripts folder. All model computations are in R/C++.
# One potential persuader per group; recipient decisions use the original
# non-convincing equilibrium. This is one response round, not a simultaneous
# equilibrium in training and persuasion. No existing manuscript is modified.
source('08_finite_pool_model.R')
convincing_compile<-local({env<-NULL;function(){
 if(is.null(env)){env<<-new.env();Rcpp::sourceCpp('10_convincing_core.cpp',env=env,
  cacheDir=file.path(tempdir(),'10_convincing_cpp'))};env
}})
convincing_density<-function(n=8192L){
 # Integrate every tail separately, rather than leaving rare competitors out
 # of a finite quantile sample. Paired cells retain exchange symmetry.
 edge<-c(0,pnorm(c(-6,-4,-2,2,4,6)),1);mass<-diff(edge)
 radical<-function(i,base){x<-numeric(length(i));f<-1/base;while(any(i>0)){x<-x+f*(i%%base);i<-i%/%base;f<-f/base};x}
 talent<-list();weight<-list();k<-0L
 for(i in seq_along(mass))for(j in i:length(mass)){
  count<-2L*ceiling(max(n/128,n*mass[i]*mass[j])/2)
  m<-if(i==j)count/2 else count;idx<-seq_len(m)+113L
  u<-cbind(edge[i]+mass[i]*radical(idx,2),edge[j]+mass[j]*radical(idx,3))
  u[]<-pmax(1e-15,pmin(1-1e-15,u));a<-matrix(exp(boredom_slog()*qnorm(u)),ncol=2)
  k<-k+1L;talent[[k]]<-rbind(a,a[,2:1]);weight[[k]]<-rep(mass[i]*mass[j]/count,2*m)
 }
 den<-list(talent=do.call(rbind,talent),weights=unlist(weight),order=n)
 stopifnot(abs(sum(den$weights)-1)<1e-12);den
}
convincing_baseline<-function(beta,phi,B,N,density){
 ref<-boredom_reference(beta,phi,B,compulsory=0);core<-finite_pool_compile()
 first<-core$finite_pool_density_cpp(density$talent,density$weights,ref,beta[1],beta[2],phi[1],phi[2],B,.5,N,160L,1L,2e-6,FALSE)
 eq<-core$finite_pool_density_cpp(density$talent,density$weights,first$mapped[1:3],beta[1],beta[2],phi[1],phi[2],B,.5,N,160L,200L,2e-6,TRUE)
 if(!eq$converged)stop('Non-convincing reference did not converge: ',eq$residual)
 if(identical(beta[1],beta[2])&&identical(phi[1],phi[2])){
  eq$moment[1]<-B/2;eq$moment[2:3]<-mean(eq$moment[2:3])
 }
 list(beta=beta,phi=phi,B=B,N=N,reference=ref,first=first,responsive=eq,order=density$order)
}
convincing_pchip_function<-function(x,y){
 h<-diff(x);d<-diff(y)/h;n<-length(x);m<-numeric(n)
 for(j in 2:(n-1))if(d[j-1]*d[j]>0){w1<-2*h[j]+h[j-1];w2<-h[j]+2*h[j-1];m[j]<-(w1+w2)/(w1/d[j-1]+w2/d[j])}
 endpoint<-function(a,b,h0,h1){v<-((2*h0+h1)*a-h0*b)/(h0+h1);if(v*a<=0)0 else if(a*b<0&&abs(v)>3*abs(a))3*a else v}
 m[1]<-endpoint(d[1],d[2],h[1],h[2]);m[n]<-endpoint(d[n-1],d[n-2],h[n-1],h[n-2])
 function(out){j<-pmax(1,pmin(n-1,findInterval(out,x)));u<-(out-x[j])/h[j]
  (2*u^3-3*u^2+1)*y[j]+(u^3-2*u^2+u)*h[j]*m[j]+(-2*u^3+3*u^2)*y[j+1]+(u^3-u^2)*h[j]*m[j+1]}
}
convincing_background_function<-function(p){
 e<-convincing_pchip_function(p$effort,log(p$moments[,1])-log(p$moments[,4]))
 m1<-convincing_pchip_function(p$effort,log(p$moments[,2]));m2<-convincing_pchip_function(p$effort,log(p$moments[,3]));B<-p$B
 function(effort){eta<-e(effort);cbind(E1=B*plogis(eta),E2=B*plogis(-eta),M1=exp(m1(effort)),M2=exp(m2(effort)))}
}
convincing_background<-function(p,effort)convincing_background_function(p)(effort)
convincing_adaptive_moments<-function(p,pressure,moment=p$baseline$responsive$moment,
 tolerance=3e-6,nx=80L){
 ans<-convincing_compile()$convincing_adaptive_moments_cpp(moment,pressure,
  p$beta[1],p$beta[2],p$phi[1],p$phi[2],p$B,p$N,boredom_slog(),nx,tolerance)
 unit<-c(p$B,pmax(1,moment[2:3]),p$B)
 ans$relative_error_estimate<-max(sweep(ans$error_estimate,2,unit,'/'))
 if(ans$relative_error_estimate>5*tolerance)
  stop('Adaptive recipient integration exceeds its global error budget: ',ans$relative_error_estimate)
 ans
}
convincing_adaptive_baseline<-function(p,tolerance=1e-6){
 encode<-function(m)c(qlogis(m[1]/p$B),log(m[2:3]))
 decode<-function(z)c(p$B*plogis(z[1]),exp(z[2:3]))
 evaluate<-function(z){
  a<-convincing_adaptive_moments(p,0,decode(z),tolerance=tolerance)
  mapped<-a$moments[1,1:3]
  list(residual=encode(mapped)-z,mapped=mapped)
 }
 z<-encode(p$baseline$responsive$moment);a<-evaluate(z)
 for(it in 1:15){
  if(max(abs(a$residual))<2e-6)break
  jac<-sapply(1:3,function(j){h<-rep(0,3);h[j]<-1e-3
   (evaluate(z+h)$residual-a$residual)/1e-3})
  step<-tryCatch(as.vector(solve(jac,-a$residual)),error=function(e)a$residual)
  scale<-1;accepted<-FALSE
  while(scale>=1/32){
   trial<-z+scale*step;b<-evaluate(trial)
   if(max(abs(b$residual))<max(abs(a$residual))){z<-trial;a<-b;accepted<-TRUE;break}
   scale<-scale/2
  }
  if(!accepted)stop('Adaptive non-convincing reference did not converge.')
 }
 if(max(abs(a$residual))>=2e-6)stop('Adaptive baseline residual remains too large.')
 p$baseline$responsive$moment<-decode(z)
 p$baseline$responsive$mapped<-a$mapped
 p$baseline$responsive$residual<-max(abs(a$residual))
 p$baseline$recipient_integration<-'adaptive';p$baseline
}
convincing_refine_recipients<-function(p,npressure=161L,tolerance=3e-6,nx=80L){
 stopifnot(npressure>=5,npressure%%2==1)
 if(!identical(p$baseline$recipient_integration,'adaptive'))
  p$baseline<-convincing_adaptive_baseline(p,tolerance=min(tolerance,1e-6))
 p$npressure<-npressure;p$effort<-seq(-p$B,p$B,length.out=npressure)
 p$pressure<-p$efficiency*p$effort
 ans<-convincing_adaptive_moments(p,p$pressure,tolerance=tolerance,nx=nx)
 p$moments<-ans$moments;p$adaptive_integration<-ans
 p$recipient_integration<-'adaptive'
 p$zero_error<-max(abs(p$moments[(npressure+1)/2,1:3]-p$baseline$responsive$moment)/
  pmax(1,p$baseline$responsive$moment))
 p$receiver_population<-convincing_compile()$convincing_receivers_cpp(p$density$talent,
  p$baseline$responsive$moment,p$pressure,p$beta[1],p$beta[2],p$phi[1],p$phi[2],p$B,p$N,p$nx)
 convincing_refresh(p)
}
convincing_prepare<-function(beta=c(.4,.4),phi=c(1,1),B=50,N=30,
 efficiency=.1,n=8192L,npressure=81L,nx=160L,baseline=NULL,aligned=FALSE,
 adaptive=!(identical(beta[1],beta[2])&&identical(phi[1],phi[2]))){
 stopifnot(length(beta)==2,length(phi)==2,all(is.finite(c(beta,phi,B,N,efficiency))),
  all(beta>0 & beta<1),all(phi>0),B>0,N>1,efficiency>0,n>=32,n%%2==0,
  npressure>=5,npressure%%2==1)
 den<-convincing_density(n)
 if(is.null(baseline))baseline<-convincing_baseline(beta,phi,B,N,den)
 stopifnot(isTRUE(all.equal(baseline$beta,beta)),isTRUE(all.equal(baseline$phi,phi)),
  baseline$B==B,baseline$N==N)
 effort<-seq(-B,B,length.out=npressure)
 pressure<-efficiency*effort;core<-convincing_compile()
 q<-core$convincing_receivers_cpp(den$talent,baseline$responsive$moment,pressure,
  beta[1],beta[2],phi[1],phi[2],B,N,nx)
 mom<-vapply(seq_along(pressure),function(k){x<-q[,k];c(sum(den$weights*x),
  sum(den$weights*boredom_w(x,den$talent[,1],beta[1],phi[1])),
  sum(den$weights*boredom_w(B-x,den$talent[,2],beta[2],phi[2])),sum(den$weights*(B-x)))},numeric(4))
 # Retain the independently evaluated zero-pressure response for residual checks.
 mid<-(npressure+1L)/2L;zero_error<-max(abs(mom[1:3,mid]-baseline$responsive$moment)/pmax(1,baseline$responsive$moment))
 if(identical(beta[1],beta[2])&&identical(phi[1],phi[2])){
  e<-(mom[1,]+rev(mom[4,]))/2;m<-(mom[2,]+rev(mom[3,]))/2
  mom[1,]<-e;mom[4,]<-rev(e);mom[2,]<-m;mom[3,]<-rev(m)
 }
 p<-list(beta=beta,phi=phi,B=B,N=N,efficiency=efficiency,n=n,npressure=npressure,
  nx=nx,aligned=aligned,energy_includes_convincing=TRUE,baseline=baseline,effort=effort,pressure=pressure,moments=t(mom),
  density=den,receiver_population=q,zero_error=zero_error)
 if(adaptive)convincing_refine_recipients(p,max(161L,npressure))else convincing_refresh(p)
}
convincing_refresh<-function(p){
 if(ncol(p$moments)!=4L)stop('This preparation uses the older integration. Recreate it with convincing_prepare() or convincing_example().')
 if(!isTRUE(p$energy_includes_convincing))stop('Recreate this preparation: convincing now contributes to tournament energy, and efficiency is per unit of time.')
 p$version<-if(identical(p$recipient_integration,'adaptive'))'adaptive-recipients-v9'else'energy-convincing-v8'
 # Anchor the zero-intervention node to the equilibrium being perturbed.
 # In particular, changing the integration order must not change which
 # tournament defines the actor's original preference.
 p$moments[(length(p$effort)+1L)/2L,]<-c(p$baseline$responsive$moment,p$B-p$baseline$responsive$moment[1])
 den<-p$density;B<-p$B;N<-p$N;q<-p$receiver_population;mid<-(length(p$effort)+1)/2
 p$actors<-convincing_optimize(p,den$talent)
 p$influence_weights<-convincing_effort_weights(p,p$actors[,3],den$weights)
 bg<-convincing_background(p,p$actors[,3])
 p$group_budget<-c(tournament1=sum(den$weights*((N-1)*bg[,1]+p$actors[,1])),
  tournament2=sum(den$weights*((N-1)*bg[,2]+p$actors[,2])),
  convincing=sum(den$weights*abs(p$actors[,3])))
 p$group_energy<-p$group_budget[1:2]+c(
  sum(den$weights*pmax(p$actors[,3],0)),sum(den$weights*pmax(-p$actors[,3],0)))
 stopifnot(abs(sum(p$group_budget)-N*B)<1e-7*N*B,
  abs(sum(p$group_energy)-N*B)<1e-7*N*B)
 q0<-q[,mid]
 p$summary<-c(mean_convincing=sum(den$weights*abs(p$actors[,3])),
  persuader_fraction=sum(den$weights*(abs(p$actors[,3])>1e-4*B)),
  group_convincing_share=sum(den$weights*abs(p$actors[,3]))/(N*B),
  mean_recipient_abs_shift=sum(den$weights*(abs(q-q0)%*%p$influence_weights)),
  mean_gain=sum(den$weights*p$actors[,7]))
 p
}
convincing_optimize<-function(p,actor,nx=40L){
 if(ncol(p$moments)!=4L)stop('This preparation uses the older integration. Recreate it with convincing_prepare() or convincing_example().')
 if(!isTRUE(p$energy_includes_convincing))stop('Recreate this preparation with convincing_prepare(): the energy accounting and efficiency convention have changed.')
 actor<-as.matrix(actor)
 stopifnot(is.numeric(actor),ncol(actor)==2L,all(is.finite(actor)),all(actor>=0),nx>=8)
 out<-convincing_compile()$convincing_focal_cpp(as.matrix(actor),p$moments,p$effort,
  p$beta[1],p$beta[2],p$phi[1],p$phi[2],p$B,p$N,nx,TRUE,isTRUE(p$aligned))
 colnames(out)<-c('x1','x2','signed_convincing','utility','baseline_x1','baseline_utility','gain',
  'other_direction_utility','direction_gap','tied_directions')
 out
}
convincing_effort_weights<-function(p,c,weights){
 # Linear interpolation of recipient responses, matching the moment interpolation.
 step<-diff(p$effort)[1];j<-pmax(1,pmin(length(p$effort)-1,findInterval(c,p$effort)))
 h<-pmax(0,pmin(1,(c-p$effort[j])/step));v<-numeric(length(p$effort))
 for(i in seq_along(c)){v[j[i]]<-v[j[i]]+weights[i]*(1-h[i]);v[j[i]+1]<-v[j[i]+1]+weights[i]*h[i]};v
}
convincing_map<-function(p,alpha=seq(0,2.5,length.out=181L)){
 actor<-as.matrix(expand.grid(alpha,alpha));out<-convincing_optimize(p,actor)
 keep<-which(p$influence_weights>1e-8)
 q<-convincing_compile()$convincing_receivers_cpp(actor,p$baseline$responsive$moment,
  p$pressure[keep],p$beta[1],p$beta[2],p$phi[1],p$phi[2],p$B,p$N,p$nx)
 q0<-out[,5]
 # Zero pressure preserves the chosen baseline strategy, even if the
 # baseline has several equally good optima and another search finds a
 # different member of that set. Switching a tie is not persuasion.
 zero<-which(p$pressure[keep]==0)
 if(length(zero))q[,zero]<-q0
 shift<-q-q0
 list(alpha=alpha,allocation=out,
  recipient_abs=as.vector(abs(shift)%*%p$influence_weights[keep]),
  recipient_signed=as.vector(shift%*%p$influence_weights[keep]))
}
convincing_example<-function(beta=c(.4,.4),phi=c(1,1),B=50,N=30,
 efficiency=.1,actor=c(1.5,.7),n=8192L,npressure=81L){
 p<-convincing_prepare(beta,phi,B,N,efficiency,n,npressure)
 print(convincing_optimize(p,matrix(actor,nrow=1)));print(p$summary);invisible(p)
}
