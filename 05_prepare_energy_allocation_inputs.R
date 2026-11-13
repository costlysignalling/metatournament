# Rebuild population energies for the identical-layout parameter-space figure.
# The fixed block is the FIRST talent-informed response to random allocations
# (E1=E2=B/2, with M1 and M2 calculated separately). The responsive block solves
# the population consistency equation. Both blocks show actual response means.
# Supplied inputs make redraws independent of Rcpp/compiler availability.
source('05_random_allocation_reference.R')
source('05_energy_allocation_quadrature.R')
if(!requireNamespace('Rcpp',quietly=TRUE))stop('Rcpp is required to rebuild energies.')
Rcpp::sourceCpp('05_energy_allocation_population.cpp',cacheDir='energy_cpp_cache')
B<-50;sdlog<-sqrt(log((1+sqrt(2))/2))
# Change beta_stride to 1 for the original complete 101 x 101 beta grid.
beta_stride<-2L
n_hermite<-32L;n_legendre<-48L;nx<-61L
threads<-max(1L,min(6L,parallel::detectCores()-2L))
input_file<-'05_equal_talent_geometry_inputs.rds'
inputs<-if(file.exists(input_file))readRDS(input_file)else NULL
# A clean run starts from the fixed reference and solves the same equilibrium
# equation. The refined energy outputs subsequently seed coexistence.
betas<-if(is.null(inputs))seq(0.01,0.99,length.out=51L) else inputs$betas[seq(1,length(inputs$betas),by=beta_stride)]
phis<-c(.66,1,1.5)
ids<-if(!is.null(inputs))seq(1,length(inputs$betas),by=beta_stride)
n<-length(betas)
q1<-energy_gh(n_hermite,sdlog);q2<-energy_gl(n_legendre)
signature<-list(input=if(is.null(inputs))'fixed-reference' else unname(tools::md5sum(input_file)),
 engine=unname(tools::md5sum('05_energy_allocation_population.cpp')),
 beta_stride=beta_stride,n_hermite=n_hermite,n_legendre=n_legendre,nx=nx)
checkpoint_dir<-'energy_allocation_checkpoints'
dir.create(checkpoint_dir,showWarnings=FALSE)
moments<-sapply(phis,function(phi)vapply(betas,function(beta)
 random_allocation_reference(beta,beta,phi,phi,B)$M1,1))
panels<-vector('list',9);dim(panels)<-c(3,3)
for(i in 1:3)for(j in i:3){
 first<-eq<-rho<-matrix(NA_real_,n,n)
 if(i==j){diag(first)<-B/2;diag(eq)<-B/2;diag(rho)<-1}
 for(r in seq_len(n)){
  cols<-if(i==j)which(seq_len(n)>r)else seq_len(n)
  if(!length(cols))next
  file<-file.path(checkpoint_dir,sprintf('phi%d_%d_row%03d.rds',i,j,r))
  record<-if(file.exists(file))readRDS(file)else NULL
  if(is.null(record)||!identical(record$signature,signature)){
   fixed_ratio<-moments[cols,j]/moments[r,i]
   initial<-if(is.null(inputs))fixed_ratio else inputs$ratios[[i,j]][ids[r],ids[cols]]
   z<-energy_population_cpp(rep(betas[r],length(cols)),betas[cols],rep(phis[i],length(cols)),
    rep(phis[j],length(cols)),fixed_ratio,initial,q1$a,q1$w,q2$x,q2$w,B,sdlog,nx=nx,threads=threads)
   if(any(z$status!=1)){
    bad<-which(z$status!=1);a<-energy_gh(64,sdlog);b<-energy_gl(96)
    retry<-energy_population_cpp(rep(betas[r],length(bad)),betas[cols[bad]],rep(phis[i],length(bad)),
     rep(phis[j],length(bad)),fixed_ratio[bad],initial[bad],a$a,a$w,b$x,b$w,B,sdlog,nx=121,threads=threads)
    for(name in names(z))z[[name]][bad]<-retry[[name]]
   }
   if(any(z$status!=1)){
    saveRDS(list(i=i,j=j,r=r,cols=cols,result=z),'energy_allocation_failed_row.rds')
    stop('Equilibrium check failed; see energy_allocation_failed_row.rds')
   }
   saveRDS(list(signature=signature,result=z),file)
  }else z<-record$result
  first[r,cols]<-z$E1_first;eq[r,cols]<-z$E1_equilibrium;rho[r,cols]<-z$ratio_equilibrium
  if(i==j){first[cols,r]<-B-z$E1_first;eq[cols,r]<-B-z$E1_equilibrium;rho[cols,r]<-1/z$ratio_equilibrium}
  if(r%%5==0||r==1)message('Energy fields phi1=',phis[i],', phi2=',phis[j],': row ',r,'/',n)
 }
 panels[[i,j]]<-list(phi1=phis[i],phi2=phis[j],first=first,equilibrium=eq,ratio_equilibrium=rho)
 if(i!=j)panels[[j,i]]<-list(phi1=phis[j],phi2=phis[i],first=B-t(first),equilibrium=B-t(eq),ratio_equilibrium=1/t(rho))
 saveRDS(list(B=B,betas=betas,phis=phis,panels=panels,signature=signature), '05_energy_allocation_inputs_partial.rds')
}
stopifnot(all(vapply(panels,function(p)all(is.finite(p$first))&&all(is.finite(p$equilibrium))&&
 min(p$first,p$equilibrium)>=0&&max(p$first,p$equilibrium)<=B,TRUE)))
output<-list(B=B,betas=betas,phis=phis,panels=panels,signature=signature,
 method='Deterministic marginal-cutoff / strategy-switch-split population integrals',
 talent=list(meanlog=0,sdlog=sdlog),
 definitions=c(A='Same leader, stronger tilt',B='Same leader, closer to equality',
 C='New leader, closer to equality',D='New leader, stronger tilt'))
saveRDS(output,'05_energy_allocation_inputs.rds')
grid<-expand.grid(beta1=betas,beta2=betas)
tab<-do.call(rbind,lapply(panels,function(p){d<-grid;d$phi1<-p$phi1;d$phi2<-p$phi2;
 d$E1_first<-as.vector(p$first);d$E2_first<-B-d$E1_first;
 d$E1_equilibrium<-as.vector(p$equilibrium);d$E2_equilibrium<-B-d$E1_equilibrium;
 k<-energy_feedback_class(d$E1_first,d$E1_equilibrium,B)
 d$feedback<-c('Neutral',LETTERS[1:4])[k+1];d}))
con<-gzfile('energy_allocation_parameter_space.csv.gz','wt');write.csv(tab,con,row.names=FALSE);close(con)
message('All population-energy fields saved.')
