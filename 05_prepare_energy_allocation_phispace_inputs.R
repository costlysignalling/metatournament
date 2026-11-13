# Rebuild the phi-space data. Rows are checkpointed, so interrupted runs resume.
# Random allocations define the fixed reference; colours show its first response.
# All optimizations compare global allocations, including specialist endpoints.
source('05_random_allocation_reference.R');source('05_energy_allocation_quadrature.R')
Rcpp::sourceCpp('05_energy_allocation_population_phispace.cpp',cacheDir='energy_phispace_cpp_cache')
B<-50;sdlog<-sqrt(log((1+sqrt(2))/2))
phis<-seq(.10,3.10,length.out=31L);panel_betas<-c(.10,.50,.90);n<-length(phis)
n_hermite<-48L;n_legendre<-32L;nx<-101L
q1<-energy_gh(n_hermite,sdlog);q2<-energy_gl(n_legendre)
threads<-max(1L,min(6L,parallel::detectCores()-2L))
signature<-list(phis=phis,panel_betas=panel_betas,B=B,sdlog=sdlog,
 engine=unname(tools::md5sum('05_energy_allocation_population_phispace.cpp')),
 n_hermite=n_hermite,n_legendre=n_legendre,nx=nx)
checkpoint_dir<-'energy_allocation_phispace_checkpoints';dir.create(checkpoint_dir,showWarnings=FALSE)
moments<-sapply(panel_betas,function(beta)vapply(phis,function(phi)
 random_allocation_reference(beta,beta,phi,phi,B)$M1,1))
panels<-vector('list',9);dim(panels)<-c(3,3)
for(i in 1:3)for(j in i:3){
 first<-eq<-rho<-matrix(NA_real_,n,n);residual<-matrix(0,n,n)
 if(i==j){diag(first)<-B/2;diag(eq)<-B/2;diag(rho)<-1}
 for(r in seq_len(n)){
  cols<-if(i==j)which(seq_len(n)>r)else seq_len(n)
  if(!length(cols))next
  file<-file.path(checkpoint_dir,sprintf('beta%d_%d_row%03d.rds',i,j,r))
  record<-if(file.exists(file))readRDS(file)else NULL
  if(is.null(record)||!identical(record$signature,signature)){
   fixed<-moments[cols,j]/moments[r,i]
   initial<-if(r>1)rho[r-1,cols]else fixed
   initial[!is.finite(initial)|initial<=0]<-fixed[!is.finite(initial)|initial<=0]
   z<-energy_population_cpp(rep(panel_betas[i],length(cols)),rep(panel_betas[j],length(cols)),
    rep(phis[r],length(cols)),phis[cols],fixed,initial,q1$a,q1$w,q2$x,q2$w,B,sdlog,nx=nx,threads=threads)
   if(any(z$status!=1)){
    bad<-which(z$status!=1);a<-energy_gh(96,sdlog);b<-energy_gl(64)
    retry<-energy_population_cpp(rep(panel_betas[i],length(bad)),rep(panel_betas[j],length(bad)),
     rep(phis[r],length(bad)),phis[cols[bad]],fixed[bad],fixed[bad],a$a,a$w,b$x,b$w,B,sdlog,nx=201,threads=threads)
    for(name in names(z))z[[name]][bad]<-retry[[name]]
   }
   if(any(z$status!=1)){saveRDS(list(i=i,j=j,r=r,cols=cols,result=z),'energy_phispace_failed_row.rds');stop('Equilibrium check failed')}
   saveRDS(list(signature=signature,result=z),file)
  }else z<-record$result
  first[r,cols]<-z$E1_first;eq[r,cols]<-z$E1_equilibrium;rho[r,cols]<-z$ratio_equilibrium;residual[r,cols]<-z$residual
  if(i==j){first[cols,r]<-B-z$E1_first;eq[cols,r]<-B-z$E1_equilibrium;rho[cols,r]<-1/z$ratio_equilibrium;residual[cols,r]<- -z$residual}
  message('Phi space beta1=',panel_betas[i],', beta2=',panel_betas[j],': row ',r,'/',n)
 }
 panels[[i,j]]<-list(beta1=panel_betas[i],beta2=panel_betas[j],first=first,equilibrium=eq,ratio_equilibrium=rho,residual=residual)
 if(i!=j)panels[[j,i]]<-list(beta1=panel_betas[j],beta2=panel_betas[i],first=B-t(first),equilibrium=B-t(eq),ratio_equilibrium=1/t(rho),residual=-t(residual))
 saveRDS(list(B=B,phis=phis,panel_betas=panel_betas,panels=panels,signature=signature),'05_energy_allocation_phispace_inputs_partial.rds')
}
stopifnot(all(vapply(panels,function(p)all(is.finite(p$first))&&all(is.finite(p$equilibrium))&&max(abs(p$residual))<2e-6,TRUE)))
output<-list(B=B,phis=phis,panel_betas=panel_betas,panels=panels,signature=signature,
 talent=list(meanlog=0,sdlog=sdlog),method='Deterministic strategy-switch and endpoint-split population integrals')
saveRDS(output,'05_energy_allocation_phispace_inputs.rds')
gr<-expand.grid(phi1=phis,phi2=phis)
tab<-do.call(rbind,lapply(panels,function(p){d<-gr;d$beta1<-p$beta1;d$beta2<-p$beta2;d$E1_first<-as.vector(p$first);d$E2_first<-B-d$E1_first;d$E1_equilibrium<-as.vector(p$equilibrium);d$E2_equilibrium<-B-d$E1_equilibrium;d}))
con<-gzfile('energy_allocation_phispace_parameter_space.csv.gz','wt');write.csv(tab,con,row.names=FALSE);close(con)
