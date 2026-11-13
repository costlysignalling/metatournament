# Resolution checks for the phi-space population calculations.
source('05_energy_allocation_quadrature.R');source('05_random_allocation_reference.R')
Rcpp::sourceCpp('05_energy_allocation_population_phispace.cpp',cacheDir='energy_phispace_cpp_cache')
B<-50;sd<-sqrt(log((1+sqrt(2))/2))
# Concave, mixed-curvature, and strongly convex examples across the new domain.
g<-data.frame(beta1=c(.1,.1,.1,.5,.5,.9,.9),beta2=c(.1,.5,.9,.5,.9,.9,.1),
 phi1=c(.4,1,3.1,1.5,.7,2.8,2.5),phi2=c(2.5,1,1.5,2.7,3.1,1.8,.4))
fixed<-vapply(seq_len(nrow(g)),function(k)with(g[k,],random_allocation_reference(beta1,beta2,phi1,phi2,B)$s),1)
a<-energy_gh(48,sd);b<-energy_gl(32)
low<-with(g,energy_population_cpp(beta1,beta2,phi1,phi2,fixed,fixed,a$a,a$w,b$x,b$w,B,sd,nx=101,threads=3))
a<-energy_gh(96,sd);b<-energy_gl(64)
high<-with(g,energy_population_cpp(beta1,beta2,phi1,phi2,fixed,low$ratio_equilibrium,a$a,a$w,b$x,b$w,B,sd,nx=201,threads=3))
stopifnot(all(low$status==1),all(high$status==1))
g$first_difference<-abs(low$E1_first-high$E1_first)
g$equilibrium_difference<-abs(low$E1_equilibrium-high$E1_equilibrium)
print(g)
# Tournament exchange must reverse both mean investments.
rev<-with(g,energy_population_cpp(beta2,beta1,phi2,phi1,1/fixed,1/high$ratio_equilibrium,a$a,a$w,b$x,b$w,B,sd,nx=201,threads=3))
g$exchange_error<-abs(rev$E1_equilibrium+high$E1_equilibrium-B)
print(g$exchange_error)
stopifnot(all(rev$status==1),max(g$first_difference,g$equilibrium_difference)<.002,max(g$exchange_error)<.002)
saveRDS(list(points=g,standard=low,high=high,reversed=rev),'05_energy_allocation_phispace_accuracy_check.rds')
