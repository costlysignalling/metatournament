# Refine the population calculations near the high-phi B/C boundary and
# the small A/B corner regions. Run from the scripts folder.
source('07_optout_model.R')
refine_optout_energy_boundaries<-function(eta=.5,n_pop=32768L,nx=401L,tol=3e-5){
 if(eta<=0||eta>=1)stop('Use 0<eta<1 for active-branch boundary refinement.')
 source('07_optout_endpoint_integration.R')
 f<-file.path('07_optout_energy_results',paste0('energy_eta_',optout_eta_tag(eta),'.rds'))
 d<-readRDS(f);beta<-d$beta;n<-length(beta)
 if(n_pop<4096||n_pop%%2||nx<51||tol<=0)stop('Use an even population >=4096, nx>=51 and tol>0.')
 core_md5<-unname(tools::md5sum('07_optout_energy_core.cpp'))
 if(exists('.optout_energy_cpp',envir=.GlobalEnv,inherits=FALSE)&&
    identical(attr(get('.optout_energy_cpp',envir=.GlobalEnv),'source_md5'),core_md5)){
  env<-get('.optout_energy_cpp',envir=.GlobalEnv)
 }else{
  env<-new.env(parent=globalenv());Rcpp::sourceCpp('07_optout_energy_core.cpp',env=env,
    cacheDir=file.path(tempdir(),'07_boundary_cpp'),rebuild=TRUE)
  attr(env,'source_md5')<-core_md5;assign('.optout_energy_cpp',env,envir=.GlobalEnv)
 }
 cells<-data.frame(i=integer(),j=integer(),a=integer(),b=integer())
 add<-function(i,j,a,b){b<-b[b>=1&b<=n];if(length(b))cells<<-rbind(cells,data.frame(i=i,j=j,a=a,b=b))}
 # B/C is separated by zero responsive energy imbalance. Use a band around
 # each grid crossing, retaining the exact reflection in equal-phi panels.
 p<-d$panels[[3,3]];delta<-p$responsive$s1-p$responsive$s2
 for(a in seq_len(n)){
  bb<-which(beta<beta[a]-.025)
  if(length(bb)<2)next
  nonzero<-which(abs(delta[a,bb])>1e-8)
  if(length(nonzero)<2)next
  left<-bb[nonzero[-length(nonzero)]];right<-bb[nonzero[-1]]
  initial<-p$first$s1-p$first$s2
  hit<-which(delta[a,left]*delta[a,right]<0 & initial[a,left]*initial[a,right]>0)
  for(h in hit)add(3,3,a,seq(left[h]-2L,right[h]+2L))
 }
 # Close to the equal-beta diagonal, small imbalances can be dominated by
 # finite winner counts. Refine noncertified cells in the endpoint-prone range.
 for(a in which(beta*d$phi[3]>=1))for(b in seq.int(max(1L,a-3L),a)){
  if(b>=a||beta[b]*d$phi[3]<1)next
  if(!is.null(p$endpoint_exact)&&isTRUE(p$endpoint_exact[a,b]))next
  add(3,3,a,b)
 }
 # A/B caps: compare absolute imbalances before and after feedback.
 for(i in 1:3){p<-d$panels[[i,3]]
  change<-abs(p$responsive$s1-p$responsive$s2)-abs(p$first$s1-p$first$s2)
  for(a in which(beta<=.18)){
   bb<-which(beta>=.90);hit<-which(change[a,bb[-length(bb)]]*change[a,bb[-1]]<=0)
   for(h in hit)add(i,3,a,bb[h]+(-1:2))
  }
  add(i,3,1,n)
 }
 cells<-unique(cells);cells<-cells[order(cells$i,cells$j,cells$a,cells$b),]
 dir.create('07_boundary_cache',showWarnings=FALSE)
 records<-vector('list',nrow(cells))
 for(row in seq_len(nrow(cells))){
  v<-cells[row,];p<-d$panels[[v$i,v$j]];bb<-beta[c(v$a,v$b)];pp<-c(p$phi1,p$phi2)
  cache<-file.path('07_boundary_cache',paste0(paste(c(optout_eta_tag(eta),d$B,v$i,v$j,v$a,v$b,n_pop,nx,format(tol,scientific=FALSE)),collapse='_'),'.rds'))
  cache_signature<-list(eta=eta,B=d$B,beta=bb,phi=pp,n_pop=n_pop,nx=nx,tol=tol,
    core=unname(tools::md5sum('07_optout_energy_core.cpp')),
    model=unname(tools::md5sum('07_optout_model.R')))
  previous<-if(file.exists(cache))readRDS(cache)else NULL
  if(!is.null(previous)&&identical(attr(previous,'signature'),cache_signature)){
   z<-previous;np<-attr(z,'population')
  }else{
   ref<-optout_fixed_scalers(bb,pp,d$B)$scaler;k<-c(p$scaler1[v$a,v$b],p$scaler2[v$a,v$b])
   for(np in unique(c(as.integer(n_pop),max(n_pop,65536L),max(n_pop,131072L)))){
    z<-env$optout_energy_grid_cpp(bb[1],pp[1],pp[2],optout_population(np),ref[1],ref[2],eta,d$B,
      nx,350L,tol,1L,matrix(k[1]),matrix(k[2]),bb[2])
    k<-c(z$scaler1[1],z$scaler2[1]);if(z$status[1]==1)break
   }
   if(z$status[1]!=1)stop('Boundary cell failed tolerance: ',paste(v,collapse=','))
   attr(z,'population')<-np;attr(z,'signature')<-cache_signature
   temporary<-paste0(cache,'.tmp');saveRDS(z,temporary)
   if(file.exists(cache))file.remove(cache)
   if(!file.rename(temporary,cache))stop('Could not save boundary cache.')
  }
  for(mode in c('first','responsive'))for(s in c('s0','s1','s2'))p[[mode]][[s]][v$a,v$b]<-z[[mode]][[s]][1]
  p$scaler1[v$a,v$b]<-z$scaler1[1];p$scaler2[v$a,v$b]<-z$scaler2[1];p$residual[v$a,v$b]<-z$residual[1];p$status[v$a,v$b]<-1L
  if(is.null(p$refinement_population))p$refinement_population<-matrix(0L,n,n)
  p$refinement_population[v$a,v$b]<-np
  if(v$i==v$j&&v$a!=v$b){
   for(mode in c('first','responsive')){p[[mode]]$s0[v$b,v$a]<-p[[mode]]$s0[v$a,v$b];p[[mode]]$s1[v$b,v$a]<-p[[mode]]$s2[v$a,v$b];p[[mode]]$s2[v$b,v$a]<-p[[mode]]$s1[v$a,v$b]}
   p$scaler1[v$b,v$a]<-z$scaler2[1];p$scaler2[v$b,v$a]<-z$scaler1[1];p$residual[v$b,v$a]<-z$residual[1];p$status[v$b,v$a]<-1L;p$refinement_population[v$b,v$a]<-np
  }
  d$panels[[v$i,v$j]]<-p
  if(v$i<v$j){swap<-function(s)list(s0=t(s$s0),s1=t(s$s2),s2=t(s$s1))
   d$panels[[v$j,v$i]]<-list(phi1=p$phi2,phi2=p$phi1,first=swap(p$first),responsive=swap(p$responsive),scaler1=t(p$scaler2),scaler2=t(p$scaler1),residual=t(p$residual),status=t(p$status),endpoint_exact=t(p$endpoint_exact),refinement_population=t(p$refinement_population))}
  records[[row]]<-c(as.list(v),list(beta=bb,phi=pp,n_pop=np,residual=z$residual[1]))
  if(row%%10==0||row==nrow(cells)){message(row,'/',nrow(cells),' boundary cells resolved');flush.console()}
 }
 d$signature$boundary_refinement<-list(n_pop=n_pop,nx=nx,tol=tol,cells=cells)
 saveRDS(d,f)
 checks_file<-paste0('07_boundary_precision_checks_eta_',optout_eta_tag(eta),'.rds')
 saveRDS(records,checks_file)
 if(eta==.5)saveRDS(records,'07_boundary_precision_checks.rds')
 invisible(d)
}
