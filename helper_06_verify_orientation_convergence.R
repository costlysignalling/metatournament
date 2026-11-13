source('06_prepare_coexistence_deterministic.R');Rcpp::sourceCpp('05_energy_allocation_population_phispace.cpp',cacheDir='energy_phispace_cpp_cache')
g<-data.frame(b1=c(.5,.5,.5,.5,.5,.9),b2=c(.9,.9,.9,.9,.9,.9),p1=c(3.1,2.675,2.625,2.5,.2,.2),p2=c(.2,.2,.325,.5,3.1,3.1))
rev<-g$p1>g$p2;n<-nrow(g);a1<-ifelse(rev,g$b2,g$b1);a2<-ifelse(rev,g$b1,g$b2);f1<-pmin(g$p1,g$p2);f2<-pmax(g$p1,g$p2)
fixed<-vapply(1:n,function(k)random_allocation_reference(a1[k],a2[k],f1[k],f2[k])$s,1);s<-sqrt(log((1+sqrt(2))/2));initial<-fixed;ans<-list()
for(order in c(48L,96L,128L,192L)){
 a<-energy_gh(order,s);b<-energy_gl(if(order==48)32 else 64)
 z<-energy_population_cpp(a1,a2,f1,f2,fixed,initial,a$a,a$w,b$x,b$w,50,s,nx=if(order==48)101 else 201,threads=6);stopifnot(all(z$status==1));ans[[as.character(order)]]<-ifelse(rev,1/z$ratio_equilibrium,z$ratio_equilibrium);initial<-z$ratio_equilibrium
}
g$rel96<-ans[['96']]/ans[['48']]-1;g$rel128<-ans[['128']]/ans[['96']]-1;g$rel192<-ans[['192']]/ans[['128']]-1
stopifnot(max(abs(g$rel128))<5e-5,max(abs(g$rel192))<5e-5)
print(g);saveRDS(list(points=g,results=ans),'06_coexistence_orientation_convergence.rds')
