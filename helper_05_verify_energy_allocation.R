# Verification at higher integration orders and a finer global allocation grid.
Rcpp::sourceCpp('05_energy_allocation_population.cpp',cacheDir='energy_cpp_cache')
source('05_energy_allocation_quadrature.R');source('05_random_allocation_reference.R')
source('05_linear_population_means.R')
x<-readRDS('05_energy_allocation_inputs.rds');base<-readRDS('05_equal_talent_geometry_inputs.rds')
B<-x$B;sd<-x$talent$sdlog
p<-rbind(c(.4,.6,.66,.66),c(.4,.6,1,1),c(.4,.6,1.5,1.5),
 c(.2,.8,.66,1.5),c(.8,.2,1.5,1),c(.6,.6,.66,1.5),c(.8,.9,1.5,1.5),
 c(.01,.99,1.5,.66),c(.99,.01,1,1.5))
# Pick one well-separated representative of each A-D regime from the data.
representatives<-list();best_margin<-rep(-Inf,4)
for(panel in x$panels){
 d0<-as.vector(panel$first)-B/2;d1<-as.vector(panel$equilibrium)-B/2
 category<-energy_feedback_class(as.vector(panel$first),as.vector(panel$equilibrium),B)
 gr<-expand.grid(beta1=x$betas,beta2=x$betas)
 margin<-pmin(abs(d0),abs(d1),abs(abs(d0)-abs(d1)))
 margin[gr$beta1<.1|gr$beta1>.9|gr$beta2<.1|gr$beta2>.9]<- -Inf
 for(k in 1:4){ids<-which(category==k);if(!length(ids))next
  z<-ids[which.max(margin[ids])]
  if(margin[z]>best_margin[k]){best_margin[k]<-margin[z];representatives[[k]]<-c(gr[z,1],gr[z,2],panel$phi1,panel$phi2)}
 }
}
stopifnot(all(is.finite(best_margin)))
p<-rbind(p,do.call(rbind,representatives))
fixed<-apply(p,1,function(v)random_allocation_reference(v[1],v[2],v[3],v[4],B)$s)
initial<-apply(p,1,function(v){i<-which.min(abs(base$betas-v[1]));j<-which.min(abs(base$betas-v[2]));
 base$ratios[[match(v[3],base$phis),match(v[4],base$phis)]][i,j]})
run<-function(ngh,ngl,nx){a<-energy_gh(ngh,sd);b<-energy_gl(ngl)
 energy_population_cpp(p[,1],p[,2],p[,3],p[,4],fixed,initial,a$a,a$w,b$x,b$w,B,sd,nx=nx,threads=4)}
a<-run(32,48,61);b<-run(64,96,121)
error_first<-abs(a$E1_first-b$E1_first);error_eq<-abs(a$E1_equilibrium-b$E1_equilibrium)
print(data.frame(p,error_first,error_eq))
stopifnot(all(a$status==1),all(b$status==1),max(error_first)<.03,max(error_eq)<.03)
last<-10:13
stopifnot(identical(energy_feedback_class(b$E1_first[last],b$E1_equilibrium[last],B),1:4))
# Independent adaptive integrals in the linear case.
linear_first<-linear_population_means(.4,.6,fixed[2],B)$E1
linear_eq<-linear_population_equilibrium(.4,.6,B)['E1']
stopifnot(abs(linear_first-b$E1_first[2])<1e-5,abs(linear_eq-b$E1_equilibrium[2])<1e-5)
for(i in 1:3)for(j in 1:3){
 stopifnot(max(abs(x$panels[[i,j]]$first+t(x$panels[[j,i]]$first)-B))<1e-12,
 max(abs(x$panels[[i,j]]$equilibrium+t(x$panels[[j,i]]$equilibrium)-B))<1e-12)
}
# Explicit algebraic unit examples, including all feedback categories.
stopifnot(identical(energy_feedback_class(c(30,30,30,30),c(35,28,23,18),B),1:4))
saveRDS(list(parameters=p,coarse=a,fine=b,error_first=error_first,error_eq=error_eq,
 classification_margins=best_margin),'05_energy_allocation_resolution_check.rds')
cat('Energy fields passed higher-resolution, exchange-symmetry, A-D and independent linear-integral checks.\n')
