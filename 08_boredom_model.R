# Local locked-time model. Run from the scripts folder.
# Density reference: N-1 expected competitors, plus the focal actor exactly.
# This is not an exact expectation of the random finite-classroom payoff.
# Numerical model calculations are in C++ through Rcpp; redraws use saved RDS.
boredom_w<-function(x,a,beta,phi)expm1(phi*log1p(a*expm1(beta*log1p(x))))
boredom_slog<-function()sqrt(log((1+sqrt(2))/2)) # median=1, SD=.5
boredom_gauss<-function(n,normal=TRUE){
 J<-matrix(0,n,n);k<-seq_len(n-1);v<-if(normal)sqrt(k/2)else k/sqrt(4*k*k-1)
 J[cbind(k,k+1)]<-v;J[cbind(k+1,k)]<-v;e<-eigen(J,symmetric=TRUE);o<-order(e$values)
 list(x=e$values[o],w=e$vectors[1,o]^2)
}
boredom_density<-function(n=32L){
 gh<-boredom_gauss(n);a<-exp(sqrt(2)*boredom_slog()*gh$x)
 list(talent=as.matrix(expand.grid(a,a)),weights=as.vector(outer(gh$w,gh$w)),order=n)
}
boredom_qmc_density<-function(n=8192L){
 # Deterministic, exchange-symmetric Halton integration in lognormal quantiles.
 stopifnot(n%%2==0,n>=32)
 radical<-function(i,base){x<-numeric(length(i));f<-1/base;while(any(i>0)){x<-x+f*(i%%base);i<-i%/%base;f<-f/base};x}
 i<-seq_len(n/2)+113L;u<-cbind(radical(i,2),radical(i,3));u<-rbind(u,u[,2:1])
 list(talent=exp(boredom_slog()*qnorm(u)),weights=rep(1/n,n),order=n)
}
boredom_reference<-function(beta,phi,B=50,compulsory=.5,order=61L){
 gh<-boredom_gauss(order);gl<-boredom_gauss(order,FALSE)
 a<-exp(sqrt(2)*boredom_slog()*gh$x);free<-(1-compulsory)*B
 x1<-compulsory*B+free*(gl$x+1)/2;x2<-free*(gl$x+1)/2
 M<-vapply(1:2,function(t)sum(gh$w*(outer(a,if(t==1)x1 else x2,
  function(a,x)boredom_w(x,a,beta[t],phi[t]))%*%gl$w)),0.)
 c(E1=compulsory*B+free/2,M1=M[1],M2=M[2])
}
boredom_compile<-function(){
 signature<-paste(unname(tools::md5sum(c('08_boredom_core.cpp','08_boredom_R_fallback.R'))),getOption('metatournament.boredom.backend','auto'),collapse=':')
 if(!exists('.boredom_cpp',envir=.GlobalEnv)||!identical(attr(get('.boredom_cpp',envir=.GlobalEnv),'signature'),signature)){
  core<-NULL
  if(getOption('metatournament.boredom.backend','auto')!='R')core<-tryCatch({
   env<-new.env(parent=globalenv());Rcpp::sourceCpp('08_boredom_core.cpp',env=env,cacheDir=file.path(tempdir(),'08_boredom_cpp'));env
  },error=function(e){warning('Compiled model unavailable; using the equivalent R implementation. ',conditionMessage(e),call.=FALSE);NULL})
  if(is.null(core)){source('08_boredom_R_fallback.R');core<-boredom_R_core()}
  attr(core,'signature')<-signature;assign('.boredom_cpp',core,envir=.GlobalEnv)
 }
 get('.boredom_cpp',envir=.GlobalEnv)
}
boredom_population<-function(beta=c(.8,.8),phi=c(1,1),B=50,compulsory=.5,N=30,
 density=boredom_density(),grid_half=80L,tolerance=1e-7,start=NULL,maxit=120L){
 stopifnot(length(beta)==2,length(phi)==2,all(beta>0 & beta<1),all(phi>0),B>0,N>1,compulsory>0,compulsory<1)
 core<-boredom_compile();ref<-boredom_reference(beta,phi,B,compulsory)
 first<-core$boredom_density_cpp(density$talent,density$weights,ref,beta[1],beta[2],phi[1],phi[2],B,compulsory,N,grid_half,1L,tolerance,FALSE)
 eq<-core$boredom_density_cpp(density$talent,density$weights,if(is.null(start))first$mapped[1:3]else start,
  beta[1],beta[2],phi[1],phi[2],B,compulsory,N,grid_half,maxit,tolerance,TRUE)
 list(beta=beta,phi=phi,B=B,compulsory=compulsory,N=N,reference=ref,first=first,responsive=eq,order=density$order)
}
boredom_map<-function(pop,mode=c('first','responsive'),alpha=seq(0,2.5,length.out=251),grid_half=160L){
 mode<-match.arg(mode);moment<-if(mode=='first')pop$reference else pop$responsive$moment
 bg<-(pop$N-1)*c(moment[1],pop$B-moment[1],moment[2:3]);actors<-as.matrix(expand.grid(alpha,alpha))
 q<-boredom_compile()$boredom_map_cpp(actors,bg,pop$beta[1],pop$beta[2],pop$phi[1],pop$phi[2],pop$B,pop$compulsory,grid_half)
 matrix(q,length(alpha),length(alpha))
}
boredom_boundaries<-function(pop,mode=c('first','responsive'),alpha=seq(0,2.5,length.out=251),grid_half=160L){
 mode<-match.arg(mode);moment<-if(mode=='first')pop$reference else pop$responsive$moment
 bg<-(pop$N-1)*c(moment[1],pop$B-moment[1],moment[2:3]);core<-boredom_compile()
 q_at<-function(x,y)core$boredom_map_cpp(cbind(x,y),bg,pop$beta[1],pop$beta[2],pop$phi[1],pop$phi[2],pop$B,pop$compulsory,grid_half)
 loq<-q_at(rep(min(alpha),length(alpha)),alpha);hiq<-q_at(rep(max(alpha),length(alpha)),alpha)
 equal_q<-if(pop$compulsory>.5)(1-pop$compulsory)*pop$B else pop$B/2
 lapply(c(1e-5,equal_q-1e-5,equal_q+1e-5,pop$B-1e-5),function(target){
  valid<-which(loq<target & hiq>=target);x<-rep(NA_real_,length(alpha))
  if(length(valid)){
   lo<-rep(min(alpha),length(valid));hi<-rep(max(alpha),length(valid));y<-alpha[valid]
   for(k in 1:24){mid<-(lo+hi)/2;above<-q_at(mid,y)>=target;hi[above]<-mid[above];lo[!above]<-mid[!above]}
   x[valid]<-(lo+hi)/2
  }
  list(x=x,y=alpha)
 })
}
