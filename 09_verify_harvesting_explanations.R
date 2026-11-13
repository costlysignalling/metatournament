source('09_harvesting_explanations.R')
s<-harvesting_slog()
stopifnot(abs(sqrt(exp(s*s)*expm1(s*s))-.5)<1e-14)
for(phi in c(.1,.66,1,1.25,2,3.1)){
 p<-c(1e-8,.01,.1,.5,.9,.99,1-1e-8)
 z<-q_harvesting_talent(p,phi)
 stopifnot(max(abs(p_harvesting_talent(z,phi)-p))<1e-12,
           abs(q_harvesting_talent(.5,phi)-expm1(phi*log(2)))<1e-12)
}
moment<-integrate(function(z)(exp(2*s*z)+2*exp(s*z))*dnorm(z),-10,10,rel.tol=1e-11)$value
stopifnot(abs(moment-(exp(2*s*s)+2*exp(s*s/2)))<1e-10)
error<-0
for(beta in c(.1,.4,.8,1))for(alpha in c(.05,1,3))for(x in c(.1,1,25,50)){
 val<-integrate(function(u)alpha*expm1(beta*log1p(u)),0,x,rel.tol=1e-10)$value
 exact<-alpha*(expm1((beta+1)*log1p(x))/(beta+1)-x)
 error<-max(error,abs(val-exact)/(1+abs(exact)))
 stopifnot(error<1e-9)
}
for(beta in c(.1,.4,.8,1))for(phi in c(.1,.66,1.25,2))for(alpha in c(.1,1,3))for(x in c(.1,1,25)){
 rate<-function(u)expm1(phi*log1p(alpha*expm1(beta*log1p(u))))
 h<-1e-4*(1+x);numeric<-(rate(x+h)-rate(x-h))/(2*h)
 exact<-phi*alpha*beta*(1+x)^(beta-1)*(1+alpha*expm1(beta*log1p(x)))^(phi-1)
 stopifnot(exact>0,abs(numeric-exact)/(1+abs(exact))<1e-7)
}
cat('Transformed CDF/quantiles, phi=2 moment, integrated-claim formula and positive curvature verified.\n')
cat('Largest relative area-formula error:',error,'\n')
