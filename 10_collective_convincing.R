# Exploratory positive-mass coalition game, NOT an individual Nash equilibrium.
# Two mirrored talent groups coordinate their members' common convincing time.
# Common beta, phi=1; exact conditional training allocation, deterministic
# integration over the existing lognormal distribution. Uniform mean exposure.
source('10_convincing_model.R')
collective_prepare<-function(beta=.4,B=50,efficiency=.1){
 stopifnot(beta>0,beta<1,B>0,efficiency>=0)
 list(beta=beta,B=B,efficiency=efficiency,s=boredom_slog())
}
collective_outcome<-function(p,v1=0,v2=0){
 stopifnot(all(is.finite(c(v1,v2))),max(abs(c(v1,v2)))<p$B)
 # r=log(alpha1/alpha2) is N(0,2s^2). Given r,
 # E[alpha1|r]=exp(r/2+s^2/4), E[alpha2|r]=exp(-r/2+s^2/4).
 # This reduces all lognormal moments to one-dimensional integration.
 sd<-sqrt(2)*p$s;delta<-p$efficiency*(v1+v2)/2
 probability<-function(lo,hi,mu=0){
  if(lo>=hi)return(0)
  if(lo>mu)pnorm(lo,mu,sd,lower.tail=FALSE)-pnorm(hi,mu,sd,lower.tail=FALSE)
  else pnorm(hi,mu,sd)-pnorm(lo,mu,sd)
 }
 moments<-function(z,v,lo,hi){
  R<-p$B-abs(v);gR<-expm1(p$beta*log1p(R))
  lower<--z-delta-(1-p$beta)*log1p(R)
  upper<--z-delta+(1-p$beta)*log1p(R)
  mass<-probability(lo,hi)
  above<-probability(max(lo,upper),hi);below<-probability(lo,min(hi,lower))
  x<-c(R*above,R*below)
  M<-gR*exp(p$s^2/2)*c(probability(max(lo,upper),hi,p$s^2),
   probability(lo,min(hi,lower),-p$s^2))
  l<-max(lo,lower);u<-min(hi,upper)
  if(l<u){
   allocation<-function(r,t){h<-(z+delta+r)/(1-p$beta)
    pmax(0,pmin(R,(R+2)*plogis(if(t==1)h else -h)-1))}
   for(t in 1:2){
    x[t]<-x[t]+integrate(function(r)allocation(r,t)*dnorm(r,0,sd),l,u,
     rel.tol=2e-9,abs.tol=1e-13)$value
    M[t]<-M[t]+exp(p$s^2/2)*integrate(function(r)
     expm1(p$beta*log1p(allocation(r,t)))*dnorm(r,if(t==1)p$s^2 else -p$s^2,sd),
     l,u,rel.tol=2e-9,abs.tol=1e-13)$value
   }
  }
  list(E=x+mass*c(max(v,0),max(-v,0)),M=M)
 }
 at_ratio<-function(z,full=FALSE){
  g1<-moments(z,v1,0,Inf);g2<-moments(z,v2,-Inf,0)
  E<-g1$E+g2$E;M<-g1$M+g2$M
  if(any(c(E,M)<=0))stop('Nonpositive moment: no population solution established.')
  residual<-z-(log(E[1])-log(E[2])+log(M[2])-log(M[1]))
  if(!full)return(residual)
  payoff<-2*c(sum(E/M*g1$M),sum(E/M*g2$M))
  list(v=c(v1,v2),delta=delta,log_ratio=z,E=E,M=M,payoff=payoff,
   residual=residual,mean_utility=mean(payoff))
 }
 z<-uniroot(at_ratio,c(-12,12),tol=1e-11)$root
 ans<-at_ratio(z,TRUE)
 stopifnot(abs(ans$residual)<1e-7,abs(sum(ans$E)-p$B)<1e-8,
  abs(ans$mean_utility-p$B)<1e-7)
 ans
}
collective_best_response<-function(p,other=0,grid=41L){
 # Both promotion directions are searched, not just the group's initial leader.
 vv<-seq(-.98*p$B,.98*p$B,length.out=grid)
 value<-function(v)collective_outcome(p,v,other)$payoff[1]
 uu<-vapply(vv,value,0.);peaks<-which(uu>=c(-Inf,head(uu,-1))&uu>=c(tail(uu,-1),-Inf))
 candidate<-cbind(v=vv,payoff=uu)
 for(j in peaks){
  fit<-optimize(value,vv[c(max(1,j-1),min(grid,j+1))],maximum=TRUE,tol=1e-7)
  candidate<-rbind(candidate,c(fit$maximum,fit$objective))
 }
 candidate[which.max(candidate[,2]),]
}
collective_symmetric_equilibrium<-function(p,interval=c(0,.6*p$B)){
 gap<-function(y)unname(collective_best_response(p,other=-y)['v'])-y
 y<-uniroot(gap,interval,tol=2e-5)$root
 # Broaden the independent global grid near the proposed fixed point.
 check<-collective_best_response(p,other=-y,grid=161L)
 ans<-collective_outcome(p,y,-y)
 ans$best_response<-check;ans$effort_residual<-unname(check['v'])-y
 ans$utility_regret<-unname(check['payoff'])-ans$payoff[1]
 stopifnot(abs(ans$effort_residual)<2e-3,abs(ans$utility_regret)<1e-6)
 ans
}
collective_example<-function(beta=.4,B=50,efficiency=.1){
 p<-collective_prepare(beta,B,efficiency)
 br0<-collective_best_response(p,other=0)
 half<-br0['v']/2
 nextbr<-collective_best_response(p,other=-half)
 report<-data.frame(case=c('No convincing','Group 1 responds alone','Half-response, mirrored'),
  v1=c(0,br0['v'],half),v2=c(0,0,-half))
 results<-lapply(seq_len(nrow(report)),function(j)collective_outcome(p,report$v1[j],report$v2[j]))
 report$payoff1<-vapply(results,function(x)x$payoff[1],0.)
 report$payoff2<-vapply(results,function(x)x$payoff[2],0.)
 report$delta<-vapply(results,function(x)x$delta,0.)
 equilibrium<-collective_symmetric_equilibrium(p)
 report<-rbind(report,data.frame(case='Symmetric coalition equilibrium',
  v1=equilibrium$v[1],v2=equilibrium$v[2],payoff1=equilibrium$payoff[1],
  payoff2=equilibrium$payoff[2],delta=equilibrium$delta))
 print(report)
 cat('Best response to the mirrored half-response:',nextbr['v'],
  '; proposed own effort:',half,'\n')
 invisible(list(parameters=p,report=report,br_to_zero=br0,
  half_response=half,br_to_half=nextbr,equilibrium=equilibrium))
}
