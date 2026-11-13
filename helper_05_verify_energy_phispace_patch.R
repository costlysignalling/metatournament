p<-readRDS('05_energy_allocation_phispace_population_patch.rds');z<-p$panel;a<-p$axis;B<-50;n<-length(a)
stopifnot(all(is.finite(z$first)),all(is.finite(z$equilibrium)),min(z$first,z$equilibrium)>=0,max(z$first,z$equilibrium)<=B,all(z$ratio_equilibrium>0))
hit<-matrix(FALSE,n-1,n-1);d0<-z$first-25;d1<-z$equilibrium-25
for(f in list(d0,d1,d1-d0,d1+d0)){
 v<-list(f[-n,-n],f[-1,-n],f[-n,-1],f[-1,-1]);lo<-Reduce(pmin,v);hi<-Reduce(pmax,v);hit<-hit|(lo<=1e-8&hi>= -1e-8&hi-lo>1e-8)
}
direct<-z$direct[-n,-n]&z$direct[-1,-n]&z$direct[-n,-1]&z$direct[-1,-1];stopifnot(all(direct[hit]))
c<-readRDS('06_coexistence_phi_inputs.rds');ids<-match(round(a,12),round(c$axis,12));q<-c$panels[[2,3]]
shared<-z$direct&q$direct[ids,ids]
error<-max(abs(z$equilibrium[shared]-50*q$E1share[ids,ids][shared]))
cat('Shared direct-node energy discrepancy:',error,'\n');stopifnot(error<1e-4)
cat('Companion energy figure: finite values, valid range, direct boundary coverage and cross-figure consistency passed.\n')
