# Exact population integration when both tournaments offer only 0/B optima.
# beta*phi>=1 and k*phi*beta<=eta are sufficient: for alpha>=1 utility
# is convex; for alpha<1 its derivative starts below eta and can turn only
# from decreasing to increasing. In either case its maximum is at 0 or B.
# The winning full-budget tournament is determined by k_t*w_t(B,alpha_t).
# Independence leaves a one-dimensional normal-density integral for each mean.
optout_endpoint_state<-function(k,beta,phi,eta=.5,B=50,sdlog=sqrt(log((1+sqrt(2))/2))){
 if(any(beta*phi<1-1e-12)||any(k*phi*beta>eta+1e-10))stop('Endpoint condition is not certified.')
 g<-expm1(beta*log1p(B));E<-M<-numeric(2)
 for(t in 1:2){u<-3-t
  a0<-expm1(log1p(eta*B/k[t])/phi[t])/g[t];lo<-log(a0)/sdlog
  density<-function(z,weighted=FALSE){a<-exp(sdlog*z);w<-expm1(phi[t]*log1p(a*g[t]));limit<-expm1(log1p(k[t]/k[u]*w)/phi[u])/g[u]
   out<-dnorm(z)*pnorm(log(limit)/sdlog);if(weighted)out<-out*w;out}
  E[t]<-B*integrate(density,lo,max(12,lo+8,phi[t]*sdlog+12),rel.tol=2e-9,abs.tol=1e-12,subdivisions=300)$value
  M[t]<-integrate(function(z)density(z,TRUE),lo,max(12,lo+8,phi[t]*sdlog+12),rel.tol=2e-9,abs.tol=1e-12,subdivisions=300)$value
 }
 list(E0=B-sum(E),E1=E[1],E2=E[2],M=M,scaler=E/M)
}
optout_endpoint_equilibrium<-function(k,beta,phi,eta=.5,B=50,tol=1e-9){
 fun<-function(v){kk<-exp(v);if(any(kk*phi*beta>eta))return(rep(NA_real_,2));st<-optout_endpoint_state(kk,beta,phi,eta,B);log(st$scaler/kk)}
 v<-log(k);err<-Inf
 for(iter in 1:60){r<-fun(v);err<-max(abs(r));if(err<tol)break
  hh<-1e-4;J<-sapply(1:2,function(t){z<-v;z[t]<-z[t]+hh;(fun(z)-r)/hh})
  step<-tryCatch(solve(J,-r),error=function(e)r/2);step<-pmax(-1,pmin(1,step))
  found<-FALSE
  for(a in 0:12){vv<-v+step*2^-a;rr<-fun(vv);if(all(is.finite(rr))&&max(abs(rr))<err){v<-vv;found<-TRUE;break}}
  if(!found){v<-v+.1*r}
 }
 if(err>=tol)stop('Endpoint residual ',err)
 list(state=optout_endpoint_state(exp(v),beta,phi,eta,B),scaler=exp(v),residual=err)
}

# Replace certified endpoint-only cells with direct one-dimensional integration.
# This removes finite-population winner-count steps near beta1=beta2 when
# phi1=phi2=1.5. Other cells retain the general global-allocation solver.
optout_refine_endpoint_cells<-function(data){
 for(i in 1:3)for(j in i:3){
  p<-data$panels[[i,j]];phi<-c(p$phi1,p$phi2);n<-length(data$beta)
  p$endpoint_exact<-matrix(FALSE,n,n)
  p$endpoint_first_exact<-matrix(FALSE,n,n)
  for(a in 1:n)for(b in 1:n){if(i==j&&a>b)next
   beta<-data$beta[c(a,b)];k<-c(p$scaler1[a,b],p$scaler2[a,b])
   if(any(beta*phi<1)||any(k*phi*beta>data$eta))next
   eq<-tryCatch(optout_endpoint_equilibrium(k,beta,phi,data$eta,data$B),error=function(e)NULL)
   if(is.null(eq))next
   sc<-optout_fixed_scalers(beta,phi,data$B)$scaler
   first_valid<-all(sc*phi*beta<=data$eta)
   first<-if(first_valid)optout_endpoint_state(sc,beta,phi,data$eta,data$B)else NULL
   # The two conditions are certified separately: a fixed-reference failure
   # does not prevent exact integration of a responsive endpoint-only state.
   for(mode in c(if(first_valid)'first','responsive')){
    st<-if(mode=='first')first else eq$state
    p[[mode]]$s0[a,b]<-st$E0/data$B;p[[mode]]$s1[a,b]<-st$E1/data$B;p[[mode]]$s2[a,b]<-st$E2/data$B
    if(i==j&&a<b){p[[mode]]$s0[b,a]<-p[[mode]]$s0[a,b];p[[mode]]$s1[b,a]<-p[[mode]]$s2[a,b];p[[mode]]$s2[b,a]<-p[[mode]]$s1[a,b]}
   }
   p$endpoint_exact[a,b]<-TRUE;p$endpoint_first_exact[a,b]<-first_valid
   if(i==j&&a<b){p$endpoint_exact[b,a]<-TRUE;p$endpoint_first_exact[b,a]<-first_valid}
   p$scaler1[a,b]<-eq$scaler[1];p$scaler2[a,b]<-eq$scaler[2];p$residual[a,b]<-eq$residual;p$status[a,b]<-1L
   if(i==j&&a<b){p$scaler1[b,a]<-eq$scaler[2];p$scaler2[b,a]<-eq$scaler[1];p$residual[b,a]<-eq$residual;p$status[b,a]<-1L}
  }
  data$panels[[i,j]]<-p
  if(i<j){swap<-function(z)list(s0=t(z$s0),s1=t(z$s2),s2=t(z$s1));data$panels[[j,i]]<-list(phi1=p$phi2,phi2=p$phi1,first=swap(p$first),responsive=swap(p$responsive),scaler1=t(p$scaler2),scaler2=t(p$scaler1),residual=t(p$residual),status=t(p$status),endpoint_exact=t(p$endpoint_exact),endpoint_first_exact=t(p$endpoint_first_exact),refinement_population=if(!is.null(p$refinement_population))t(p$refinement_population)else NULL)}
 }
 data$signature$endpoint_integration<-TRUE
 data$signature$n_endpoint_cells<-sum(vapply(data$panels,function(p)sum(p$endpoint_exact),1))
 data$signature$n_endpoint_first_cells<-sum(vapply(data$panels,function(p)sum(p$endpoint_first_exact),1))
 data
}
