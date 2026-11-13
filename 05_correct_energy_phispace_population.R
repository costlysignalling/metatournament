# Keep the companion phi-energy figure consistent with Script 06's converged
# integration direction in the mixed beta=.50,.90 panel. Other panels are reused.
source('06_prepare_coexistence_deterministic.R')
correct_energy_phispace_population<-function(threads=6L){
 coex<-prepare_coexistence('phi');base<-readRDS('05_energy_allocation_phispace_refined_inputs.rds')
 file<-'05_energy_allocation_phispace_population_patch.rds'
 sig<-list(base=unname(tools::md5sum('05_energy_allocation_phispace_refined_inputs.rds')),
  coex=unname(tools::md5sum('06_coexistence_phi_inputs.rds')),
  engine=unname(tools::md5sum('05_energy_allocation_population_phispace.cpp')),version=1L)
 if(file.exists(file)){z<-readRDS(file);if(identical(z$signature,sig))return(z)}
 Rcpp::sourceCpp('05_energy_allocation_population_phispace.cpp',cacheDir='energy_phispace_cpp_cache')
 B<-base$B;sd<-base$talent$sdlog;gh<-energy_gh(48,sd);gl<-energy_gl(32)
 coarse<-base$phis[seq(1,length(base$phis),by=2)];fine<-base$phis
 cp<-coex$panels[[2,3]];root<-'energy_phispace_population_patch_checkpoints';dir.create(root,showWarnings=FALSE)
 first<-eq<-rho<-NULL;oldaxis<-NULL
 for(stage in 1:2){
  axis<-if(stage==1)coarse else fine;n<-length(axis);gr<-expand.grid(x=axis,y=axis)
  if(stage==1){first<-eq<-rho<-matrix(NA_real_,n,n);direct<-matrix(FALSE,n,n);need<-matrix(TRUE,n,n)}else{
   first<-matrix(coex_interp(first,gr$x,gr$y,oldaxis),n);eq<-matrix(coex_interp(eq,gr$x,gr$y,oldaxis),n)
   rho<-exp(matrix(coex_interp(log(rho),gr$x,gr$y,oldaxis),n));direct<-matrix(FALSE,n,n);direct[seq(1,n,by=2),seq(1,n,by=2)]<-TRUE
  }
  moments<-sapply(c(.5,.9),function(beta)vapply(axis,function(phi)random_allocation_reference(beta,beta,phi,phi,B)$M1,1))
  boundary_nodes<-function(){
   d0<-first-B/2;d1<-eq-B/2;hit<-matrix(FALSE,n-1,n-1)
   for(f in list(d0,d1,d1-d0,d1+d0)){
    z<-list(f[-n,-n],f[-1,-n],f[-n,-1],f[-1,-1]);lo<-Reduce(pmin,z);hi<-Reduce(pmax,z)
    hit<-hit|(lo<=1e-8&hi>= -1e-8&hi-lo>1e-8)
   }
   mask<-matrix(FALSE,n,n);cells<-which(hit,arr.ind=TRUE)
   for(k in seq_len(nrow(cells))){r<-cells[k,1];c<-cells[k,2];mask[r:(r+1),c:(c+1)]<-TRUE};mask
  }
  if(stage==2)need<-boundary_nodes()&!direct
  for(pass in 1:8){
   for(r in which(rowSums(need)>0)){
    wanted<-which(need[r,]);fn<-file.path(root,sprintf('n%03d_row%03d.rds',n,r))
    rec<-if(file.exists(fn))tryCatch(readRDS(fn),error=function(e)NULL)else NULL
    rsig<-c(sig,list(axis=axis))
    if(is.null(rec)||!identical(rec$signature,rsig))rec<-list(signature=rsig,cols=integer(),result=NULL)
    missing<-setdiff(wanted,rec$cols)
    if(length(missing)){
     x<-rep(axis[r],length(missing));y<-axis[missing];rev<-x>y
     f1<-ifelse(rev,y,x);f2<-ifelse(rev,x,y);b1<-ifelse(rev,.9,.5);b2<-ifelse(rev,.5,.9)
     fixed<-moments[missing,2]/moments[r,1];fixed<-ifelse(rev,1/fixed,fixed)
     initial<-exp(coex_interp(log(cp$rho),x,y,coex$axis));initial<-ifelse(rev,1/initial,initial)
     z<-energy_population_cpp(b1,b2,f1,f2,fixed,initial,gh$a,gh$w,gl$x,gl$w,B,sd,nx=101,threads=threads)
     stopifnot(all(z$status==1))
     z$E1_first<-ifelse(rev,B-z$E1_first,z$E1_first);z$E1_equilibrium<-ifelse(rev,B-z$E1_equilibrium,z$E1_equilibrium)
     z$ratio_equilibrium<-ifelse(rev,1/z$ratio_equilibrium,z$ratio_equilibrium)
     rec$cols<-c(rec$cols,missing)
     if(is.null(rec$result))rec$result<-z else for(k in names(z))rec$result[[k]]<-c(rec$result[[k]],z[[k]])
     saveRDS(rec,fn)
    }
    c<-rec$cols;first[r,c]<-rec$result$E1_first;eq[r,c]<-rec$result$E1_equilibrium;rho[r,c]<-rec$result$ratio_equilibrium;direct[r,c]<-TRUE
    if(r%%10==0||r==1)message('Companion phi-energy correction n=',n,' row ',r)
   }
   if(stage==1)break
   need<-boundary_nodes()&!direct;if(!any(need))break
  }
  if(stage==2)stopifnot(!any(need))
  oldaxis<-axis
 }
 p<-list(beta1=.5,beta2=.9,first=first,equilibrium=eq,ratio_equilibrium=rho,direct=direct)
 out<-list(signature=sig,panel=p,axis=axis,method='Lower-phi outer integration; same continuous-population reference as Script 06')
 saveRDS(out,file);out
}
