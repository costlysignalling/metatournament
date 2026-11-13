# Increase population resolution only at cells that fail the requested residual.
# The separate beta2_values argument lets the C++ solver evaluate one point.
optout_repair_energy_panel<-function(p,beta,eta,B,engine,tol=.0001,nx=401L){
 bad<-which(p$status!=1,arr.ind=TRUE)
 if(!nrow(bad))return(p)
 if(is.null(p$refinement_population))p$refinement_population<-matrix(0L,length(beta),length(beta))
 for(row in seq_len(nrow(bad))){a<-bad[row,1];b<-bad[row,2];if(p$status[a,b]==1)next;bb<-beta[c(a,b)];pp<-c(p$phi1,p$phi2)
  ref<-optout_fixed_scalers(bb,pp,B)$scaler;k<-c(p$scaler1[a,b],p$scaler2[a,b]);ok<-FALSE
  for(np in c(8192L,32768L,65536L)){
   message('Refine cell: beta=',paste(bb,collapse=','),' phi=',paste(pp,collapse=','),' n=',np)
   z<-engine$optout_energy_grid_cpp(bb[1],pp[1],pp[2],optout_population(np),ref[1],ref[2],eta,B,nx,300L,tol,1L,matrix(k[1]),matrix(k[2]),bb[2])
   k<-c(z$scaler1[1],z$scaler2[1])
   if(z$status[1]==1){ok<-TRUE;break}
  }
  if(!ok)stop('Unresolved energy cell after population refinement: beta=',paste(bb,collapse=','))
  for(mode in c('first','responsive'))for(v in c('s0','s1','s2'))p[[mode]][[v]][a,b]<-z[[mode]][[v]][1]
  p$refinement_population[a,b]<-np
  p$scaler1[a,b]<-k[1];p$scaler2[a,b]<-k[2];p$residual[a,b]<-z$residual[1];p$status[a,b]<-1L
  if(p$phi1==p$phi2 && a!=b){p$refinement_population[b,a]<-np;for(mode in c('first','responsive')){p[[mode]]$s0[b,a]<-p[[mode]]$s0[a,b];p[[mode]]$s1[b,a]<-p[[mode]]$s2[a,b];p[[mode]]$s2[b,a]<-p[[mode]]$s1[a,b]};p$scaler1[b,a]<-k[2];p$scaler2[b,a]<-k[1];p$residual[b,a]<-z$residual[1];p$status[b,a]<-1L}
 }
 p
}
