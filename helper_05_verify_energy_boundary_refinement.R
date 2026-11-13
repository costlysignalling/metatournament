# Check transition-split integrals and adaptive boundary coverage.
Rcpp::sourceCpp('05_energy_allocation_population_refined.cpp',cacheDir='energy_refined_cpp_cache')
source('05_energy_allocation_quadrature.R');source('05_random_allocation_reference.R');source('05_linear_population_means.R')
x<-readRDS('05_energy_allocation_refined_inputs.rds');base<-readRDS('05_equal_talent_geometry_inputs.rds')
B<-x$B;sd<-x$talent$sdlog
p<-readRDS('05_energy_allocation_resolution_check.rds')$parameters
p<-rbind(p,c(.6,.65,1.5,1.5),c(.6078,.8528,.66,1.5),c(.6568,.7646,1.5,1.5))
fixed<-apply(p,1,function(v)random_allocation_reference(v[1],v[2],v[3],v[4],B)$s)
initial<-apply(p,1,function(v)base$ratios[[match(v[3],base$phis),match(v[4],base$phis)]][
 which.min(abs(base$betas-v[1])),which.min(abs(base$betas-v[2]))])
run<-function(n1,n2,nx){a<-energy_gh(n1,sd);b<-energy_gl(n2)
 energy_population_cpp(p[,1],p[,2],p[,3],p[,4],fixed,initial,a$a,a$w,b$x,b$w,B,sd,nx=nx,threads=4)}
a<-run(48,32,101);b<-run(96,64,201)
error_first<-abs(a$E1_first-b$E1_first);error_eq<-abs(a$E1_equilibrium-b$E1_equilibrium)
print(data.frame(p,error_first,error_eq))
stopifnot(all(a$status==1),all(b$status==1),max(error_first)<.001,max(error_eq)<.001)
stopifnot(abs(b$E1_first[2]-linear_population_means(.4,.6,fixed[2],B)$E1)<1e-5,
 abs(b$E1_equilibrium[2]-linear_population_equilibrium(.4,.6,B)['E1'])<1e-5)
coverage<-list()
for(i in 1:3)for(j in 1:3){
 panel<-x$panels[[i,j]];n<-nrow(panel$first);d0<-panel$first-B/2;d1<-panel$equilibrium-B/2
 hit<-matrix(FALSE,n-1,n-1)
 for(f in list(d0,d1,d1-d0,d1+d0)){
  corners<-list(f[-n,-n],f[-1,-n],f[-n,-1],f[-1,-1]);lo<-Reduce(pmin,corners);hi<-Reduce(pmax,corners)
  hit<-hit|(lo<=1e-8 & hi>= -1e-8 & hi-lo>1e-8)
 }
 good<-panel$direct[-n,-n]&panel$direct[-1,-n]&panel$direct[-n,-1]&panel$direct[-1,-1]
 stopifnot(all(good[hit]),max(abs(panel$first+t(x$panels[[j,i]]$first)-B))<1e-12,
  max(abs(panel$equilibrium+t(x$panels[[j,i]]$equilibrium)-B))<1e-12)
 coverage[[paste(i,j,sep='_')]]<-c(boundary_cells=sum(hit),direct_nodes=sum(panel$direct))
}
saveRDS(list(parameters=p,coarse=a,fine=b,error_first=error_first,error_eq=error_eq,coverage=coverage),
 '05_energy_boundary_refinement_check.rds')
cat('Higher-order integrals, independent linear checks, exchange symmetry and direct boundary coverage passed.\n')
