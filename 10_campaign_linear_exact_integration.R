# Exact one-dimensional reduction when phi1=phi2=1 in the population limit.
# Talents are independent lognormal variables with median 1 and SD .5.
source('10_convincing_model.R')
campaign_linear_compile<-local({core<-NULL;function(){if(is.null(core)){
 core<<-new.env();Rcpp::sourceCpp('10_campaign_linear_core.cpp',env=core,
 cacheDir=file.path(tempdir(),'10_campaign_linear_cpp'))};core}})
campaign_linear_recipient_moments<-function(beta,B=50,delta=0,log_ratio=0,tolerance=1e-8){
 s<-boredom_slog();sd<-sqrt(2)*s
 lower<-log(beta[2]/beta[1])+(beta[2]-1)*log1p(B)-log_ratio-delta
 upper<-log(beta[2]/beta[1])-(beta[1]-1)*log1p(B)-log_ratio-delta
 demand<-function(z){q<-z+log_ratio+delta
  if(beta[1]==beta[2])return(pmax(0,pmin(B,(B+2)*plogis(q/(1-beta[1]))-1)))
  lo<-rep(0,length(z));hi<-rep(B,length(z))
  for(k in 1:40){mid<-(lo+hi)/2
   positive<-q+log(beta[1]/beta[2])+(beta[1]-1)*log1p(mid)-(beta[2]-1)*log1p(B-mid)>0
   lo[positive]<-mid[positive];hi[!positive]<-mid[!positive]}
  (lo+hi)/2
 }
 integ<-function(f)integrate(function(z)f(z,demand(z))*dnorm(z,sd=sd),lower,upper,
  rel.tol=tolerance,abs.tol=1e-12,subdivisions=400L)$value
 E1<-B*pnorm(upper,sd=sd,lower.tail=FALSE)+integ(function(z,x)x)
 E2<-B*pnorm(lower,sd=sd)+integ(function(z,x)B-x)
 M1<-((1+B)^beta[1]-1)*exp(s*s/2)*pnorm(upper,mean=s*s,sd=sd,lower.tail=FALSE)+
  integ(function(z,x)exp(z/2+s*s/4)*expm1(beta[1]*log1p(x)))
 M2<-((1+B)^beta[2]-1)*exp(s*s/2)*pnorm(lower,mean=-s*s,sd=sd)+
  integ(function(z,x)exp(-z/2+s*s/4)*expm1(beta[2]*log1p(B-x)))
 c(E1=E1,M1=M1,M2=M2,E2=E2)
}
campaign_linear_global_mean<-function(beta,B=50,efficiency=.1,npressure=161L,nx=48L,
 tolerance=2e-5,pilot=201L){
 stopifnot(length(beta)==2,all(beta>0&beta<1),npressure>=5,npressure%%2==1,
  B>0,efficiency>0,tolerance>0,pilot>=21,nx>=8)
 reference<-function(z){m<-campaign_linear_recipient_moments(beta,B,log_ratio=z)
  z-log(m[1]/m[2]*m[3]/m[4])}
 ratio<-if(beta[1]==beta[2])0 else uniroot(reference,c(-5,5),tol=1e-9)$root
 effort<-seq(-B,B,length.out=npressure)
 moments<-t(vapply(efficiency*effort,function(d)campaign_linear_recipient_moments(beta,B,d,ratio),numeric(4)))
 core<-campaign_linear_compile();sd<-sqrt(2)*boredom_slog();limit<-8*sd
 evaluate<-function(z)core$campaign_linear_focal_cpp(cbind(exp(z/2),exp(-z/2)),moments,effort,beta,B)
 zz<-seq(-limit,limit,length.out=pilot);o<-evaluate(zz)
 region<-ifelse(abs(o[,3])<1e-6*B,0,sign(o[,3]));changes<-which(diff(region)!=0)
 breaks<-c(-limit,sd*c(-6,-4,-2,0,2,4,6),limit)
 for(k in changes){lo<-zz[k];hi<-zz[k+1];left<-region[k]
  for(it in 1:25){mid<-(lo+hi)/2;cc<-evaluate(mid)[1,3];side<-if(abs(cc)<1e-6*B)0 else sign(cc)
   if(side==left)lo<-mid else hi<-mid}
  breaks<-c(breaks,(lo+hi)/2)
 }
 breaks<-sort(unique(breaks));cache<-new.env(hash=TRUE)
 values<-function(z){key<-sprintf('%.16g',z)
  missing<-!vapply(key,exists,FALSE,envir=cache,inherits=FALSE)
  if(any(missing)){ans<-evaluate(z[missing]);ii<-which(missing)
   for(j in seq_along(ii))assign(key[ii[j]],ans[j,],envir=cache)}
  do.call(rbind,lapply(key,get,envir=cache,inherits=FALSE))}
 fun<-list(function(o)abs(o[,3])/B,function(o)pmax(o[,3],0)/B,
  function(o)pmax(-o[,3],0)/B,function(o)o[,1]/B,function(o)o[,2]/B)
 result<-error<-numeric(5);messages<-character()
 for(j in 1:5)for(k in 1:(length(breaks)-1)){
  out<-integrate(function(z)fun[[j]](values(z))*dnorm(z,sd=sd),breaks[k],breaks[k+1],
   rel.tol=tolerance,abs.tol=tolerance/(length(breaks)-1),subdivisions=600L,stop.on.error=FALSE)
  if(out$message!="OK")messages<-c(messages,out$message)
  result[j]<-result[j]+out$value;error[j]<-error[j]+out$abs.error
 }
 if(any(!is.finite(result))||max(error)>4*tolerance)stop('Outer integration error estimate exceeds tolerance: ',max(error))
 names(result)<-c('campaign','campaign1','campaign2','training1','training2')
 direct_result<-result
 # Integrate each component as a cross-check, then use the exact expectation
 # identities instead of propagating independent quadrature roundoff.
 result['campaign']<-result['campaign1']+result['campaign2']
 result['training2']<-1-result['campaign']-result['training1']
 list(result=c(result,baseline_residual=unname(abs(reference(ratio)))),outer_error_estimate=error,direct_result=direct_result,
  omitted_tail_bound=2*pnorm(-8),log_ratio=ratio,beta=beta,phi=c(1,1),
  npressure=npressure,nx=nx,tolerance=tolerance,breaks=breaks,messages=unique(messages),method='one-dimensional adaptive talent-ratio integration')
}
