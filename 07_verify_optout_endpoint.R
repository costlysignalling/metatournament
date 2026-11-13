source('07_optout_model.R')
source('07_optout_endpoint_integration.R')
beta<-c(.95,.99);phi<-c(1.5,1.5);eta<-.5;B<-50
ref<-optout_fixed_scalers(beta,phi,B)$scaler
# A broad, deterministic talent grid tests entries, maxxers, and switching.
a<-as.matrix(expand.grid(exp(seq(-3,3,length.out=31)),exp(seq(-3,3,length.out=31))))
z<-optout_allocate(a,beta,phi,ref,eta,B,nscan=513L)
u1<-ref[1]*optout_w(B,a[,1],beta[1],phi[1]);u2<-ref[2]*optout_w(B,a[,2],beta[2],phi[2])
stopifnot(max(abs(z[,'utility']-pmax(eta*B,u1,u2)))<1e-7)
eq<-optout_endpoint_equilibrium(ref,beta,phi,eta,B)
stopifnot(eq$residual<1e-8,abs(eq$state$E0+eq$state$E1+eq$state$E2-B)<1e-10)
eq<-optout_endpoint_equilibrium(rep(ref[2],2),rep(.99,2),phi,.8,B)
a0<-expm1(log1p(.8*B/eq$scaler[1])/phi[1])/expm1(.99*log1p(B))
closed<-(1-pnorm(log(a0)/optout_slog())^2)/2
stopifnot(abs(eq$state$E1/B-closed)<1e-10)
cat('Global allocation, budget, continuous-scaler, and exact symmetric probability checks passed.\n')
