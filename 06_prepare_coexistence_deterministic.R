# Shared deterministic preparation for Script 06's beta and phi planes.
# Fixed E: uniform random allocations, E1=E2=B/2, separate calculated M1,M2.
# Responsive E: continuous-population equilibrium, including responsive M.
# HDR existence is tested geometrically, not by sampled strategy frequencies.
source('05_random_allocation_reference.R');source('05_energy_allocation_quadrature.R')
coex_names<-c('S2max','S2div','S1div','S1max')
coex_levels<-c(0,1e-4,1e-4,0)
coex_interp<-function(z,x,y,axis){
 i<-pmax(1L,pmin(length(axis)-1L,findInterval(x,axis)))
 j<-pmax(1L,pmin(length(axis)-1L,findInterval(y,axis)))
 tx<-(x-axis[i])/(axis[i+1]-axis[i]);ty<-(y-axis[j])/(axis[j+1]-axis[j])
 (1-tx)*(1-ty)*z[cbind(i,j)]+tx*(1-ty)*z[cbind(i+1,j)]+(1-tx)*ty*z[cbind(i,j+1)]+tx*ty*z[cbind(i+1,j+1)]
}
coex_crossing_nodes<-function(panel,n){
 hit<-matrix(FALSE,n-1,n-1)
 for(mode in c('fixed','responsive'))for(k in seq_along(coex_names)){
  f<-panel[[mode]][[coex_names[k]]]-coex_levels[k]
  zz<-list(f[-n,-n],f[-1,-n],f[-n,-1],f[-1,-1]);lo<-Reduce(pmin,zz);hi<-Reduce(pmax,zz)
  hit<-hit|(lo<=0 & hi>=0 & hi-lo>1e-12)
 }
 mask<-matrix(FALSE,n,n);cells<-which(hit,arr.ind=TRUE)
 for(k in seq_len(nrow(cells))){r<-cells[k,1];c<-cells[k,2];mask[r:(r+1),c:(c+1)]<-TRUE};mask
}
coex_swap<-function(p){
 out<-p
 for(mode in c('fixed','responsive')){
  out[[mode]]<-list(S2max=t(p[[mode]]$S1max),S2div=t(p[[mode]]$S1div),S1div=t(p[[mode]]$S2div),S1max=t(p[[mode]]$S2max))
 }
 out$rho<-1/t(p$rho);out$E1share<-1-t(p$E1share);out$direct<-t(p$direct)
 out
}
prepare_coexistence<-function(plane=c('beta','phi'),force=FALSE,threads=max(1L,min(6L,parallel::detectCores()-2L))){
 plane<-match.arg(plane);file<-paste0('06_coexistence_',plane,'_inputs.rds');migration<-NULL;previous_input_hash<-NULL
 if(file.exists(file)&&!force){
  saved<-readRDS(file)
  deps<-c(engine='06_coexistence_hdr_scores.cpp',population_engine='05_energy_allocation_population_phispace.cpp',
   population=if(plane=='beta')'05_energy_allocation_refined_inputs.rds'else '05_energy_allocation_phispace_refined_inputs.rds')
  valid<-all(vapply(names(deps),function(k)!file.exists(deps[k])||identical(saved$signature[[k]],unname(tools::md5sum(deps[k]))),logical(1)))
  if(valid){
   if(plane=='phi'&&is.null(saved$orientation_policy)){
    migration<-saved;previous_input_hash<-unname(tools::md5sum(file))
   }else return(saved)
  }
  message('Numerical inputs changed: rebuilding ',plane,' coexistence fields.')
 }
 Rcpp::sourceCpp('06_coexistence_hdr_scores.cpp',cacheDir='coexistence_cpp_cache')
 Rcpp::sourceCpp('05_energy_allocation_population_phispace.cpp',cacheDir='energy_phispace_cpp_cache')
 if(plane=='beta'){
  pop_file<-'05_energy_allocation_refined_inputs.rds';pop<-readRDS(pop_file)
  pop_axis<-pop$betas;panel_vals<-pop$phis;axes<-list(pop_axis[seq(1,length(pop_axis),by=2)],pop_axis)
 }else{
  pop_file<-'05_energy_allocation_phispace_refined_inputs.rds';pop<-readRDS(pop_file)
  pop_axis<-pop$phis;panel_vals<-pop$panel_betas
  coarse<-pop_axis[seq(1,length(pop_axis),by=2)]
  axes<-list(coarse,pop_axis,sort(c(pop_axis,(head(pop_axis,-1)+tail(pop_axis,-1))/2)))
 }
 B<-pop$B;sd<-pop$talent$sdlog;hc<- -sd^2;hr<-sd*sqrt(qchisq(.99,2,ncp=2*sd^2))
 q1<-energy_gh(48,sd);q2<-energy_gl(32)
 u<-seq(0,1,length.out=101);xt<-sort(unique(c(seq(0,B,length.out=181),B*u^2/2,B-B*u^2/2)))
 xl<-seq(0,B/2,length.out=121);xr<-seq(B/2,B,length.out=121);theta<-seq(-pi/2,pi/2,length.out=61)
 root<-file.path('coexistence_deterministic_cache',plane);dir.create(root,recursive=TRUE,showWarnings=FALSE)
 signature<-list(version=if(plane=='phi')2L else 1L,population=unname(tools::md5sum(pop_file)),
  engine=unname(tools::md5sum('06_coexistence_hdr_scores.cpp')),population_engine=unname(tools::md5sum('05_energy_allocation_population_phispace.cpp')))
 panels<-NULL;old_axis<-NULL
 for(stage in seq_along(axes)){
  axis<-axes[[stage]];n<-length(axis);gr<-expand.grid(x=axis,y=axis)
  moments<-sapply(panel_vals,function(v)vapply(axis,function(a)
   if(plane=='beta')random_allocation_reference(a,a,v,v,B)$M1 else random_allocation_reference(v,v,a,a,B)$M1,1))
  previous<-panels;panels<-vector('list',9);dim(panels)<-c(3,3)
  for(i in 1:3)for(j in i:3){
   ep<-pop$panels[[i,j]]
   if(!is.null(migration)&&!(i==2&&j==3)){
    ids<-match(round(axis,12),round(migration$axis,12));stopifnot(!anyNA(ids))
    p<-migration$panels[[i,j]]
    for(mode in c('fixed','responsive'))p[[mode]]<-lapply(p[[mode]],function(z)z[ids,ids,drop=FALSE])
    p$rho<-p$rho[ids,ids,drop=FALSE];p$E1share<-p$E1share[ids,ids,drop=FALSE];p$direct<-p$direct[ids,ids,drop=FALSE]
    panels[[i,j]]<-p;if(i!=j)panels[[j,i]]<-coex_swap(p)
    next
   }
   if(stage==1){
    p<-list(fixed=setNames(lapply(1:4,function(k)matrix(NA_real_,n,n)),coex_names),
     responsive=setNames(lapply(1:4,function(k)matrix(NA_real_,n,n)),coex_names),
     rho=matrix(NA_real_,n,n),E1share=matrix(NA_real_,n,n),direct=matrix(FALSE,n,n))
    need<-matrix(TRUE,n,n)
   }else{
    old<-previous[[i,j]]
    p<-old
    for(mode in c('fixed','responsive'))p[[mode]]<-lapply(old[[mode]],function(z)matrix(coex_interp(z,gr$x,gr$y,old_axis),n))
    p$rho<-exp(matrix(coex_interp(log(old$rho),gr$x,gr$y,old_axis),n))
    p$E1share<-matrix(coex_interp(old$E1share,gr$x,gr$y,old_axis),n)
    p$direct<-matrix(FALSE,n,n);ids<-seq(1,n,by=2);p$direct[ids,ids]<-old$direct
    need<-coex_crossing_nodes(p,n)&!p$direct
   }
   calc_rows<-function(mask){
    if(i==j)mask[lower.tri(mask)]<-FALSE
    for(r in which(rowSums(mask)>0)){
     wanted<-which(mask[r,]);fn<-file.path(root,sprintf('n%03d_panel%d_%d_row%03d.rds',n,i,j,r))
     rec<-if(file.exists(fn)&&!force)tryCatch(readRDS(fn),error=function(e)NULL)else NULL
     sig<-c(signature,list(stage=stage,axis=axis))
     pop_cache<-NULL
     if(!is.null(rec)&&!identical(rec$signature,sig)){
      old_sig<-rec$signature;new_sig<-sig;old_sig$engine<-NULL;new_sig$engine<-NULL
      if(plane=='phi'){old_sig$version<-NULL;new_sig$version<-NULL}
      if(identical(old_sig,new_sig))pop_cache<-rec
     }
     if(is.null(rec)||!identical(rec$signature,sig))rec<-list(signature=sig,cols=integer(),fixed=NULL,responsive=NULL,rho=numeric(),E1share=numeric(),direct=logical())
     missing<-setdiff(wanted,rec$cols)
     if(length(missing)){
      x<-rep(axis[r],length(missing));y<-axis[missing]
      b1<-if(plane=='beta')x else rep(panel_vals[i],length(x));b2<-if(plane=='beta')y else rep(panel_vals[j],length(x))
      f1<-if(plane=='phi')x else rep(panel_vals[i],length(x));f2<-if(plane=='phi')y else rep(panel_vals[j],length(x))
      fixed<-moments[missing,j]/moments[r,i]
      rho<-exp(coex_interp(log(ep$ratio_equilibrium),x,y,pop_axis))
      share<-coex_interp(ep$equilibrium/B,x,y,pop_axis)
      ix<-match(round(x,12),round(pop_axis,12));iy<-match(round(y,12),round(pop_axis,12))
      direct<-rep(FALSE,length(x));ok<-which(!is.na(ix)&!is.na(iy));direct[ok]<-ep$direct[cbind(ix[ok],iy[ok])]
      reverse_orientation<-plane=='phi' & i==2 & j==3 & f1>f2
      direct[reverse_orientation]<-FALSE
      if(!is.null(pop_cache)){
       cached<-match(missing,pop_cache$cols);reuse<-which(!is.na(cached)&pop_cache$direct[cached]&!(pop_cache$signature$version==1L&reverse_orientation))
       rho[reuse]<-pop_cache$rho[cached[reuse]];share[reuse]<-pop_cache$E1share[cached[reuse]];direct[reuse]<-TRUE
      }
      if(i==j){
       symmetric<-which(abs(x-y)<1e-12)
       rho[symmetric]<-1;share[symmetric]<-.5;direct[symmetric]<-TRUE
      }
      # Stage one supplies a coarse preview. All final boundary nodes must
      # subsequently use directly validated population moments and scores.
      if(stage>1||(plane=='phi'&&i==2&&j==3)){
       todo<-which(!direct)
       if(length(todo)){
        rev<-reverse_orientation[todo]
        cb1<-ifelse(rev,b2[todo],b1[todo]);cb2<-ifelse(rev,b1[todo],b2[todo])
        cf1<-ifelse(rev,f2[todo],f1[todo]);cf2<-ifelse(rev,f1[todo],f2[todo])
        cr<-ifelse(rev,1/fixed[todo],fixed[todo]);initial<-ifelse(rev,1/rho[todo],rho[todo])
        z<-energy_population_cpp(cb1,cb2,cf1,cf2,cr,initial,q1$a,q1$w,q2$x,q2$w,B,sd,nx=101,threads=threads)
        if(any(z$status!=1))stop('Population equilibrium check failed at ',plane,' stage ',stage,' panel ',i,',',j,' row ',r)
        rho[todo]<-ifelse(rev,1/z$ratio_equilibrium,z$ratio_equilibrium)
        share[todo]<-ifelse(rev,1-z$E1_equilibrium/B,z$E1_equilibrium/B);direct[todo]<-TRUE
       }
      }
      sfixed<-coexistence_scores_cpp(b1,b2,f1,f2,fixed,xt,xl,xr,theta,hc,hr,B,threads)
      sresp<-coexistence_scores_cpp(b1,b2,f1,f2,rho,xt,xl,xr,theta,hc,hr,B,threads)
      rec$cols<-c(rec$cols,missing);rec$rho<-c(rec$rho,rho);rec$E1share<-c(rec$E1share,share);rec$direct<-c(rec$direct,direct)
      if(is.null(rec$fixed)){rec$fixed<-sfixed;rec$responsive<-sresp}else for(k in coex_names){rec$fixed[[k]]<-c(rec$fixed[[k]],sfixed[[k]]);rec$responsive[[k]]<-c(rec$responsive[[k]],sresp[[k]])}
      saveRDS(rec,fn)
     }
     cols<-rec$cols
     for(mode in c('fixed','responsive'))for(k in coex_names)p[[mode]][[k]][r,cols]<<-rec[[mode]][[k]]
     p$rho[r,cols]<<-rec$rho;p$E1share[r,cols]<<-rec$E1share;p$direct[r,cols]<<-rec$direct
     if(i==j){
      for(mode in c('fixed','responsive'))for(k in seq_along(coex_names))p[[mode]][[coex_names[k]]][cols,r]<<-rec[[mode]][[coex_names[5-k]]]
      p$rho[cols,r]<<-1/rec$rho;p$E1share[cols,r]<<-1-rec$E1share;p$direct[cols,r]<<-rec$direct
     }
     if(r%%10==0||r==1)message('Coexistence ',plane,' n=',n,' panel ',i,',',j,' row ',r)
    }
   }
   for(pass in 1:8){
    calc_rows(need)
    if(stage==1)break
    need<-coex_crossing_nodes(p,n)&!p$direct
    if(i==j)need[lower.tri(need)]<-FALSE
    if(!any(need))break
   }
   if(stage>1)stopifnot(!any(need))
   panels[[i,j]]<-p;if(i!=j)panels[[j,i]]<-coex_swap(p)
   saveRDS(list(plane=plane,axis=axis,panel_vals=panel_vals,panels=panels,signature=signature),paste0('06_coexistence_',plane,'_partial.rds'))
  }
  old_axis<-axis
 }
 out<-list(plane=plane,axis=axis,panel_vals=panel_vals,panels=panels,B=B,talent=pop$talent,hdr=.99,levels=coex_levels,signature=signature,
  method='Separate random-allocation M_t at fixed E; deterministic population equilibrium; global HDR geometric presence scores; direct boundary refinement')
 if(plane=='phi'){
  out$orientation_policy<-'Lower-phi outer integration in the mixed beta=.5,.9 panel; exact exchange elsewhere'
  if(!is.null(previous_input_hash))out$compatible_geometry_input<-previous_input_hash
 }
 saveRDS(out,file);out
}
