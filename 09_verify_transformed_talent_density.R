source('09_transformed_talent_density.R')
for(phi in c(.66,1,1.25,2)){
 mass<-integrate(function(z)d_harvesting_talent(z,phi),0,Inf,rel.tol=1e-9)$value
 upper<-q_harvesting_talent(.995,phi)
 shown<-integrate(function(z)d_harvesting_talent(z,phi),0,upper,rel.tol=1e-9)$value
 stopifnot(abs(mass-1)<1e-8,abs(shown-.995)<1e-8)
 cat('phi=',phi,' total mass=',mass,' mass to 99.5% quantile=',shown,'\n')
}
x<-seq(.01,5,length.out=1001)
stopifnot(max(abs(d_harvesting_talent(x,1)-dlnorm(x,0,harvesting_slog())))<1e-12)
cat('At phi=1 the density coincides with the original lognormal.\n')
