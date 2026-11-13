# Read the saved, independent high-resolution checks and verify their signs.
# For a fresh calculation use recalculate=TRUE (requires Rcpp/Rtools).
verify_optout_energy_boundaries<-function(recalculate=FALSE){
 f<-'07_energy_boundary_validation.rds'
 if(recalculate){
  source('07_optout_model.R')
  env<-new.env(parent=globalenv());Rcpp::sourceCpp('07_optout_energy_core.cpp',env=env,
    cacheDir=file.path(tempdir(),'07_boundary_verify'),rebuild=TRUE)
  d<-readRDS('07_optout_energy_results/energy_eta_0.5.rds');corners<-bc<-list()
  calculate<-function(bb,pp,np,nx,k){ref<-optout_fixed_scalers(bb,pp)$scaler
   env$optout_energy_grid_cpp(bb[1],pp[1],pp[2],optout_population(np),ref[1],ref[2],.5,50,nx,400L,3e-5,1L,matrix(k[1]),matrix(k[2]),bb[2])}
  for(i in 1:3)for(np in c(16384L,65536L)){
   p<-d$panels[[i,3]];z<-calculate(c(.01,.99),c(d$phi[i],1.5),np,401L,c(p$scaler1[1,length(d$beta)],p$scaler2[1,length(d$beta)]))
   corners[[paste(i,np)]]<-list(phi=c(d$phi[i],1.5),np=np,result=z)
  }
  for(ix in list(c(45,31),c(45,32),c(45,33),c(59,26),c(59,27),c(59,28))){
   p<-d$panels[[3,3]];bb<-d$beta[ix];z<-calculate(bb,c(1.5,1.5),65536L,601L,c(p$scaler1[ix[1],ix[2]],p$scaler2[ix[1],ix[2]]))
   bc[[paste(ix,collapse='_')]]<-list(indices=ix,beta=bb,n_pop=65536L,nx=601L,result=z)
  }
  saveRDS(list(corners=corners,bc=bc),f)
 }
 checks<-readRDS(f)
 delta<-function(z,mode)z[[mode]]$s1[1]-z[[mode]]$s2[1]
 for(rec in c(checks$corners,checks$bc)){
  z<-rec$result;stopifnot(z$status[1]==1,z$residual[1]<=3e-5)
  for(mode in c('first','responsive'))stopifnot(abs(z[[mode]]$s0[1]+z[[mode]]$s1[1]+z[[mode]]$s2[1]-1)<1e-10)
 }
 for(rec in checks$corners){a<-delta(rec$result,'first');b<-delta(rec$result,'responsive');stopifnot(a*b>0,abs(b)<abs(a))}
 for(pair in list(c('45_31','45_33'),c('59_26','59_28'))){
  stopifnot(delta(checks$bc[[pair[1]]]$result,'responsive')>0,
            delta(checks$bc[[pair[2]]]$result,'responsive')<0)
 }
 cat('Corner-region B classifications, B/C brackets and resource accounting verified.\n')
 invisible(checks)
}
