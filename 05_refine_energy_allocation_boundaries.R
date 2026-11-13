# Higher precision around the A-D region boundaries.
# Refines the beta mesh locally, increases quadrature/optimization accuracy,
# and splits at specialist/diversifier transitions as well as at B/2.
# Untouched background nodes are interpolated from the base field and are
# explicitly flagged; every displayed category boundary is checked to have
# directly recalculated nodes on both sides. Rebuild is resumable by beta row.
source('05_energy_allocation_quadrature.R');source('05_random_allocation_reference.R')
Rcpp::sourceCpp('05_energy_allocation_population_refined.cpp',cacheDir='energy_refined_cpp_cache')
base_file<-'05_energy_allocation_inputs.rds'
base<-readRDS(base_file);B<-base$B;sd<-base$talent$sdlog;phis<-base$phis
coarse<-base$betas;nc<-length(coarse)
betas<-sort(c(coarse,(head(coarse,-1)+tail(coarse,-1))/2));n<-length(betas)
q1<-energy_gh(48,sd);q2<-energy_gl(32)
threads<-max(1L,min(6L,parallel::detectCores()-2L))
signature<-list(base_md5=unname(tools::md5sum(base_file)),
 engine=unname(tools::md5sum('05_energy_allocation_population_refined.cpp')),
 n_hermite=48L,n_legendre_per_interval=32L,nx=101L,n_beta=n)
checkpoint_dir<-'energy_allocation_refined_checkpoints'
dir.create(checkpoint_dir,showWarnings=FALSE)
interp<-function(z){
 gr<-expand.grid(x=betas,y=betas)
 i<-pmax(1L,pmin(nc-1L,findInterval(gr$x,coarse)));j<-pmax(1L,pmin(nc-1L,findInterval(gr$y,coarse)))
 tx<-(gr$x-coarse[i])/(coarse[i+1]-coarse[i]);ty<-(gr$y-coarse[j])/(coarse[j+1]-coarse[j])
 matrix((1-tx)*(1-ty)*z[cbind(i,j)]+tx*(1-ty)*z[cbind(i+1,j)]+
 (1-tx)*ty*z[cbind(i,j+1)]+tx*ty*z[cbind(i+1,j+1)],n)
}
crossing_cells<-function(first,eq){
 m<-nrow(first);d0<-first-B/2;d1<-eq-B/2;hit<-matrix(FALSE,m-1,m-1)
 for(f in list(d0,d1,d1-d0,d1+d0)){
  z<-list(f[-m,-m],f[-1,-m],f[-m,-1],f[-1,-1]);lo<-Reduce(pmin,z);hi<-Reduce(pmax,z)
  hit<-hit|(lo<=1e-8 & hi>= -1e-8 & hi-lo>1e-8)
 };hit
}
initial_selection<-function(p){
 hit<-crossing_cells(p$first,p$equilibrium);mask<-matrix(FALSE,n,n);cells<-which(hit,arr.ind=TRUE)
 for(k in seq_len(nrow(cells))){r<-2*cells[k,1]-1;c<-2*cells[k,2]-1;mask[r:(r+2),c:(c+2)]<-TRUE}
 mask
}
boundary_nodes<-function(first,eq){
 hit<-crossing_cells(first,eq);mask<-matrix(FALSE,n,n);cells<-which(hit,arr.ind=TRUE)
 for(k in seq_len(nrow(cells))){r<-cells[k,1];c<-cells[k,2];mask[r:(r+1),c:(c+1)]<-TRUE}
 mask
}
moments<-sapply(phis,function(phi)vapply(betas,function(beta)
 random_allocation_reference(beta,beta,phi,phi,B)$M1,1))
panels<-vector('list',9);dim(panels)<-c(3,3)
for(i in 1:3)for(j in i:3){
 p<-base$panels[[i,j]];first<-interp(p$first);eq<-interp(p$equilibrium)
 rho<-exp(interp(log(p$ratio_equilibrium)));direct<-matrix(FALSE,n,n)
 if(i==j){diag(first)<-B/2;diag(eq)<-B/2;diag(rho)<-1;diag(direct)<-TRUE}
 need<-initial_selection(p)
 run_rows<-function(mask){
  if(i==j)mask[lower.tri(mask,diag=TRUE)]<-FALSE
  rows<-which(rowSums(mask)>0)
  for(r in rows){
   wanted<-which(mask[r,]);file<-file.path(checkpoint_dir,sprintf('phi%d_%d_row%03d.rds',i,j,r))
   record<-if(file.exists(file))readRDS(file)else NULL
   if(is.null(record)||!identical(record$signature,signature))record<-list(cols=integer(),result=NULL)
   missing<-setdiff(wanted,record$cols)
   if(length(missing)){
    fixed<-moments[missing,j]/moments[r,i]
    z<-energy_population_cpp(rep(betas[r],length(missing)),betas[missing],rep(phis[i],length(missing)),
     rep(phis[j],length(missing)),fixed,rho[r,missing],q1$a,q1$w,q2$x,q2$w,B,sd,nx=101,threads=threads)
    if(any(z$status!=1)){
     bad<-which(z$status!=1);a<-energy_gh(96,sd);b<-energy_gl(64)
     retry<-energy_population_cpp(rep(betas[r],length(bad)),betas[missing[bad]],rep(phis[i],length(bad)),
      rep(phis[j],length(bad)),fixed[bad],rho[r,missing[bad]],a$a,a$w,b$x,b$w,B,sd,nx=201,threads=threads)
     for(name in names(z))z[[name]][bad]<-retry[[name]]
    }
    if(any(z$status!=1)){
     saveRDS(list(i=i,j=j,r=r,cols=missing,result=z),'energy_refined_failed_row.rds')
     stop('Refined equilibrium check failed; see energy_refined_failed_row.rds')
    }
    record$cols<-c(record$cols,missing)
    if(is.null(record$result))record$result<-z else for(name in names(z))record$result[[name]]<-c(record$result[[name]],z[[name]])
    record$signature<-signature;saveRDS(record,file)
   }
   cols<-record$cols;z<-record$result
   first[r,cols]<<-z$E1_first;eq[r,cols]<<-z$E1_equilibrium;rho[r,cols]<<-z$ratio_equilibrium;direct[r,cols]<<-TRUE
   if(i==j){first[cols,r]<<-B-z$E1_first;eq[cols,r]<<-B-z$E1_equilibrium;rho[cols,r]<<-1/z$ratio_equilibrium;direct[cols,r]<<-TRUE}
   if(r%%10==0||r==rows[1])message('Refined energy boundaries phi1=',phis[i],', phi2=',phis[j],': beta row ',r,'/',n)
  }
 }
 for(pass in 1:5){
  run_rows(need)
  need<-boundary_nodes(first,eq)&!direct
  if(i==j)need[lower.tri(need,diag=TRUE)]<-FALSE
  if(!any(need))break
  message('Additional boundary nodes needed: ',sum(need))
 }
 stopifnot(!any(need))
 panels[[i,j]]<-list(phi1=phis[i],phi2=phis[j],first=first,equilibrium=eq,ratio_equilibrium=rho,direct=direct)
 if(i!=j)panels[[j,i]]<-list(phi1=phis[j],phi2=phis[i],first=B-t(first),equilibrium=B-t(eq),ratio_equilibrium=1/t(rho),direct=t(direct))
 saveRDS(list(B=B,betas=betas,phis=phis,panels=panels,signature=signature),'05_energy_allocation_refined_inputs_partial.rds')
}
output<-base;output$betas<-betas;output$panels<-panels;output$signature<-signature
output$method<-'Adaptive beta mesh; marginal-cutoff and endpoint/strategy-switch-split deterministic integrals'
output$note<-'Every category boundary uses directly refined nodes; background interpolation is flagged by direct.'
saveRDS(output,'05_energy_allocation_refined_inputs.rds')
gr<-expand.grid(beta1=betas,beta2=betas)
tab<-do.call(rbind,lapply(panels,function(p){d<-gr;d$phi1<-p$phi1;d$phi2<-p$phi2;
 d$E1_first<-as.vector(p$first);d$E2_first<-B-d$E1_first;
 d$E1_equilibrium<-as.vector(p$equilibrium);d$E2_equilibrium<-B-d$E1_equilibrium;
 d$direct_high_accuracy<-as.vector(p$direct)
 k<-energy_feedback_class(d$E1_first,d$E1_equilibrium,B);d$feedback<-c('Neutral',LETTERS[1:4])[k+1];d}))
con<-gzfile('energy_allocation_refined_parameter_space.csv.gz','wt');write.csv(tab,con,row.names=FALSE);close(con)
message('Refined energy boundaries and direct-node flags saved.')
