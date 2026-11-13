# Thermodynamic limit of ONE population-wide popularizer, first response.
# A fixed positive fraction of the population hears this actor. This is
# macroscopic reach, not a finite personal audience or a coalition game.
# This exact reduction assumes phi1=phi2=1, equal training beta, independent
# lognormal talents (median=1, SD=.5). No finite-population sample is used.
# source('10_popularizing_population.R')
# popularizing_predict_phi1(actor=c(2,.5), beta=.4, B=50, efficiency=.1)
popularizing_slog<-function()sqrt(log((1+sqrt(2))/2))
popularizing_logadd<-function(a,b){m<-pmax(a,b);ifelse(is.infinite(m)&m<0,-Inf,m+log(exp(a-m)+exp(b-m)))}
popularizing_raw_moments_phi1<-function(delta,beta=.4,B=50,s=popularizing_slog(),rel.tol=1e-9){
 stopifnot(length(delta)==1,is.finite(delta),length(beta)==1,beta>0,beta<=1,B>0,s>0)
 sdD<-sqrt(2)*s
 if(beta==1){
  le<-log(B)+c(pnorm(delta/sdD,log.p=TRUE),pnorm(-delta/sdD,log.p=TRUE))
  lm<-log(B)+s*s/2+c(pnorm((delta+s*s)/sdD,log.p=TRUE),pnorm((s*s-delta)/sdD,log.p=TRUE))
 }else{
  h<-1-beta;L<-log1p(B);lo<--delta-h*L;hi<--delta+h*L
  g<-function(x)expm1(beta*log1p(x))
  x<-function(z,reverse=FALSE)pmax(0,pmin(B,(B+2)*plogis(if(reverse)-z else z)-1))
  # The Gaussian density is rescaled before integration; normal tails and
  # additions stay in log space, retaining small expected competitor claims.
  middle<-function(mu,claim=FALSE,reverse=FALSE){
   peak<-dnorm(pmax(lo,pmin(hi,mu)),mu,sdD,log=TRUE)
   integrand<-function(z){xx<-x(z,reverse);if(claim)xx<-g(xx)
    xx*exp(dnorm(h*z-delta,mu,sdD,log=TRUE)-peak)}
   ans<-integrate(integrand,-L,L,rel.tol=rel.tol,abs.tol=0,subdivisions=300L)$value
   if(ans<=0)return(-Inf)
   log(h)+peak+log(ans)
  }
  le<-c(popularizing_logadd(log(B)+pnorm(hi,0,sdD,lower.tail=FALSE,log.p=TRUE),middle(0)),
   popularizing_logadd(log(B)+pnorm(lo,0,sdD,log.p=TRUE),middle(0,reverse=TRUE)))
  lm<-s*s/2+c(
   popularizing_logadd(log(g(B))+pnorm(hi,s*s,sdD,lower.tail=FALSE,log.p=TRUE),middle(s*s,claim=TRUE)),
   popularizing_logadd(log(g(B))+pnorm(lo,-s*s,sdD,log.p=TRUE),middle(-s*s,claim=TRUE,reverse=TRUE)))
 }
 list(E=exp(le),M=exp(lm),scaler=exp(le-lm),logE=le,logM=lm)
}
popularizing_moments_phi1<-function(delta,beta=.4,B=50,audience_fraction=1,s=popularizing_slog(),rel.tol=1e-9){
 stopifnot(audience_fraction>=0,audience_fraction<=1)
 ans<-popularizing_raw_moments_phi1(delta,beta,B,s,rel.tol)
 if(audience_fraction<1){
  base<-popularizing_raw_moments_phi1(0,beta,B,s,rel.tol)
  ans$logE<-popularizing_logadd(log(audience_fraction)+ans$logE,log1p(-audience_fraction)+base$logE)
  ans$logM<-popularizing_logadd(log(audience_fraction)+ans$logM,log1p(-audience_fraction)+base$logM)
  ans$E<-exp(ans$logE);ans$M<-exp(ans$logM);ans$scaler<-exp(ans$logE-ans$logM)
 }
 ans
}
popularizing_recipient_phi1<-function(log_talent_ratio,delta=0,beta=.4,B=50){
 if(beta==1)return(ifelse(log_talent_ratio+delta>0,B,ifelse(log_talent_ratio+delta<0,0,B/2)))
 pmax(0,pmin(B,(B+2)*plogis((log_talent_ratio+delta)/(1-beta))-1))
}
popularizing_training_phi1<-function(actor,scaler,remaining,beta){
 if(remaining<=0)return(c(x1=0,x2=0,utility=0))
 ratio<-log(actor[1])+log(scaler[1])-log(actor[2])-log(scaler[2])
 if(beta==1)x1<-if(ratio>0)remaining else if(ratio<0)0 else remaining/2
 else x1<-pmax(0,pmin(remaining,(remaining+2)*plogis(ratio/(1-beta))-1))
 x2<-remaining-x1
 c(x1=x1,x2=x2,utility=sum(actor*scaler*expm1(beta*log1p(c(x1,x2)))))
}
popularizing_predict_phi1<-function(actor=c(2,.5),beta=.4,B=50,efficiency=.1,
 audience_fraction=1,s=popularizing_slog(),nsearch=161L,rel.tol=1e-9){
 stopifnot(length(actor)==2,all(is.finite(actor)),all(actor>0),efficiency>=0,nsearch>=9)
 # Homogeneity in talent gives exactly the same choices after rescaling actor.
 scale<-sum(actor);a<-actor/scale;cache<-new.env(parent=emptyenv())
 evaluate<-function(v){
  key<-sprintf('%.16g',v)
  if(!exists(key,envir=cache,inherits=FALSE)){
   m<-popularizing_moments_phi1(efficiency*v,beta,B,audience_fraction,s,rel.tol)
   assign(key,popularizing_training_phi1(a,m$scaler,B-abs(v),beta),envir=cache)
  }
  get(key,envir=cache,inherits=FALSE)
 }
 grid<-sort(unique(c(seq(-B,B,length.out=nsearch),0)))
 val<-vapply(grid,function(v)evaluate(v)['utility'],0.)
 candidates<-grid
 k<-which(val[2:(length(val)-1)]>=val[1:(length(val)-2)]&val[2:(length(val)-1)]>=val[3:length(val)])+1
 for(j in k)candidates<-c(candidates,optimize(function(v)evaluate(v)['utility'],grid[c(j-1,j+1)],maximum=TRUE,tol=1e-8)$maximum)
 vlu<-vapply(candidates,function(v)evaluate(v)['utility'],0.);best<-which.max(vlu)
 ties<-which(vlu>=max(vlu)-1e-10*(1+max(vlu)))
 if(length(ties)>1){side<-if(actor[1]>=actor[2])1 else -1;best<-ties[which.max(side*candidates[ties])]}
 v<-candidates[best];base<-evaluate(0)
 if(evaluate(v)['utility']<=base['utility']+1e-10*(1+base['utility']))v<-0
 ans<-evaluate(v)
 list(actor=actor,beta=beta,B=B,efficiency=efficiency,audience_fraction=audience_fraction,
  x1=unname(ans['x1']),x2=unname(ans['x2']),signed_convincing=v,convincing_time=abs(v),
  utility=unname(scale*ans['utility']),baseline_utility=unname(scale*base['utility']),
  relative_utility=unname(ans['utility']/base['utility']),
  population=popularizing_moments_phi1(efficiency*v,beta,B,audience_fraction,s,rel.tol))
}
popularizing_linear_threshold<-function(B=50,s=popularizing_slog(),audience_fraction=1){
 stopifnot(B>0,s>0,audience_fraction>=0,audience_fraction<=1)
 if(audience_fraction==0)return(Inf)
 response<-(dnorm(0)/pnorm(0)-dnorm(s/sqrt(2))/pnorm(s/sqrt(2)))/(sqrt(2)*s)
 1/(B*audience_fraction*response)
}
