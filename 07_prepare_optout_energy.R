# Rebuild population energy for the opt-out parameter-space figures.
# Source this file from the scripts folder, then call prepare_optout_energy().
# Data and numerical settings are saved separately from plotting.
source('07_optout_model.R')
source('07_optout_endpoint_integration.R')
source('07_repair_optout_energy.R')
prepare_optout_energy <- function(eta_values=.5,n_beta=61L,
    n_pop=4096L,nx=301L,tol=.0001,threads=6L,B=50) {
 if(any(eta_values<0|eta_values>1))stop('Use 0<=eta<=1 for active-branch comparisons.')
 core_md5<-unname(tools::md5sum('07_optout_energy_core.cpp'))
 if(!exists('.optout_energy_cpp',envir=.GlobalEnv,inherits=FALSE) ||
    !identical(attr(get('.optout_energy_cpp',envir=.GlobalEnv),'source_md5'),core_md5)) {
  env<-new.env(parent=globalenv())
  Rcpp::sourceCpp('07_optout_energy_core.cpp',env=env,
    cacheDir=file.path(tempdir(),'07_energy_cpp'),rebuild=TRUE)
  attr(env,'source_md5')<-core_md5
  assign('.optout_energy_cpp',env,envir=.GlobalEnv)
 }
 beta<-seq(.01,.99,length.out=n_beta);phi<-c(.66,1,1.5)
 population<-optout_population(n_pop)
 ref<-sapply(phi,function(p) vapply(beta,function(b) optout_fixed_scalers(b,p,B)$scaler[1],1))
 dir.create('07_optout_energy_results',showWarnings=FALSE)
 signature<-list(version=2L,beta=beta,phi=phi,n_pop=n_pop,nx=nx,tol=tol,B=B,
   core=unname(tools::md5sum('07_optout_energy_core.cpp')),
   endpoint=unname(tools::md5sum('07_optout_endpoint_integration.R')))
 output<-character()
 for(eta in eta_values) {
  panels<-vector('list',9);dim(panels)<-c(3,3)
  for(j in 1:3)for(i in 1:j) {
   f<-file.path('07_optout_energy_results',sprintf('eta_%s_phi_%s_%s_nb%d_np%d_nx%d.rds',
     optout_eta_tag(eta),phi[i],phi[j],n_beta,n_pop,nx))
   cached<-if(file.exists(f)) readRDS(f) else NULL
   if(!is.null(cached)&&identical(cached$signature,signature)) z<-cached$data else {
    message('eta=',eta,' phi1=',phi[i],' phi2=',phi[j])
    tm<-system.time(z<-.GlobalEnv$.optout_energy_cpp$optout_energy_grid_cpp(beta,phi[i],phi[j],
      population,ref[,i],ref[,j],eta,B,nx,180L,tol,threads))
    message('  ',round(tm[3],1),' seconds; ',sum(z$status==0),' unresolved cells')
    saveRDS(list(signature=signature,data=z),f)
   }
   panels[[i,j]]<-list(phi1=phi[i],phi2=phi[j],first=z$first,responsive=z$responsive,
     scaler1=z$scaler1,scaler2=z$scaler2,residual=z$residual,status=z$status)
   if(i<j) {
    swap<-function(s) list(s0=t(s$s0),s1=t(s$s2),s2=t(s$s1))
    panels[[j,i]]<-list(phi1=phi[j],phi2=phi[i],first=swap(z$first),responsive=swap(z$responsive),
      scaler1=t(z$scaler2),scaler2=t(z$scaler1),residual=t(z$residual),status=t(z$status))
   }
  }
  data<-list(beta=beta,phi=phi,eta=eta,B=B,panels=panels,signature=signature)
  f<-file.path('07_optout_energy_results',paste0('energy_eta_',optout_eta_tag(eta),'.rds'))
  data<-optout_refine_endpoint_cells(data)
  for(j in 1:3)for(i in 1:j){
   p<-optout_repair_energy_panel(data$panels[[i,j]],beta,eta,B,
     get('.optout_energy_cpp',envir=.GlobalEnv),tol,max(401L,nx))
   data$panels[[i,j]]<-p
   if(i<j){swap<-function(z)list(s0=t(z$s0),s1=t(z$s2),s2=t(z$s1))
    data$panels[[j,i]]<-list(phi1=p$phi2,phi2=p$phi1,first=swap(p$first),responsive=swap(p$responsive),
      scaler1=t(p$scaler2),scaler2=t(p$scaler1),residual=t(p$residual),status=t(p$status),endpoint_exact=t(p$endpoint_exact),refinement_population=if(!is.null(p$refinement_population))t(p$refinement_population)else NULL)}
  }
  saveRDS(data,f);output<-c(output,f)
 }
 invisible(output)
}
