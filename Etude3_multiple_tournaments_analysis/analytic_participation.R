#Exact strategy probabilities for phi=1 and 2 or 3 identical tournaments.
#Talents are independent lognormal draws. All probabilities are deterministic
#one-dimensional quadratures of analytically derived region boundaries.

talent_log_sd<-function(median=1,sd=0.5){
  sqrt(log((1+sqrt(1+4*(sd/median)^2))/2))
}

participation_phi1<-function(beta,B,T=3L,slog=talent_log_sd()){
  stopifnot(length(beta)==1L,beta>0,beta<1,B>0,T%in%c(2L,3L),slog>0)
  L<-(1-beta)*log1p(B)
  if(T==2L){
    p1<-2*pnorm(-L/(sqrt(2)*slog))
    return(c(one=p1,two=1-p1))
  }
  p1<-3*integrate(function(z){
    dnorm(z)*pnorm(z-L/slog)^2
  },-Inf,Inf,rel.tol=1e-11,abs.tol=1e-13)$value

  #Condition on tournament 3 having the smallest talent. For
  #D1=log(alpha1/alpha3), D2=log(alpha2/alpha3):
  #var(Di)=2*slog^2, corr(D1,D2)=1/2.
  #All three are active iff exp(D1/(1-beta))+exp(D2/(1-beta))<B+2.
  p3<-3*integrate(function(t){
    d<-L*t
    upper<-(1-beta)*log(pmax(1,B+2-exp(log1p(B)*t)))
    sdcond<-sqrt(3/2)*slog
    L*dnorm(d,sd=sqrt(2)*slog)*
      (pnorm((upper-d/2)/sdcond)-pnorm(-d/(2*sdcond)))
  },0,1,rel.tol=1e-11,abs.tol=1e-13)$value
  ans<-c(one=p1,two=1-p1-p3,three=p3)
  stopifnot(all(ans>=-1e-11),abs(sum(ans)-1)<1e-12)
  pmax(ans,0)
}

#Independent derivation for checking the three-active probability: integrate
#over the minimum log-talent and the second labelled log-talent, rather than
#using the bivariate distribution of log-ratios.
check_three_via_talents<-function(beta,B,slog=talent_log_sd()){
  L<-(1-beta)*log1p(B)
  inner<-function(z){
    integrate(function(y){
      upper<-z+(1-beta)/slog*
        log(pmax(1,B+2-exp(slog*(y-z)/(1-beta))))
      dnorm(y)*(pnorm(upper)-pnorm(z))
    },z,z+L/slog,rel.tol=1e-8,abs.tol=1e-12)$value
  }
  3*integrate(function(z){dnorm(z)*vapply(z,inner,1)},
              -8,8,rel.tol=1e-8,abs.tol=1e-10)$value
}

run_checks<-function(){
  for(beta in c(0.4,0.8)){
    a<-participation_phi1(beta,50)
    independent<-check_three_via_talents(beta,50)
    stopifnot(abs(a['three']-independent)<1e-8)
    cat('Independent integral check passed: beta=',beta,'\n')
  }
  for(beta in c(0.2,0.4,0.8,0.95)){
    p2<-participation_phi1(beta,50,T=2)
    p3<-participation_phi1(beta,50,T=3)
    stopifnot(p3['three']<=p2['two']+1e-12)
    stopifnot(participation_phi1(beta,75)['three']>=p3['three']-1e-12)
  }
  stopifnot(participation_phi1(0.4,1e-7)['one']>0.99999)
  stopifnot(participation_phi1(0.999999,50)['one']>0.99999)
}

if(sys.nframe()==0L){
  run_checks()
  beta_grid<-c(0.1,0.2,0.4,0.6,0.8,0.9,0.95,0.99)
  result<-do.call(rbind,lapply(beta_grid,function(beta){
    p2<-participation_phi1(beta,50,T=2)
    p3fixed<-participation_phi1(beta,50,T=3)
    p3scaled<-participation_phi1(beta,75,T=3)
    data.frame(beta=beta,phi=1,
      two_tournaments_B50_both=p2['two'],
      three_B50_one=p3fixed['one'],three_B50_two=p3fixed['two'],
      three_B50_three=p3fixed['three'],
      three_B75_one=p3scaled['one'],three_B75_two=p3scaled['two'],
      three_B75_three=p3scaled['three'],row.names=NULL)
  }))
  write.csv(result,'participation_phi1.csv',row.names=FALSE)
  print(result,digits=7,row.names=FALSE)
  crossover<-uniroot(function(beta){
    p<-participation_phi1(beta,50)
    p['two']-p['three']
  },c(0.8,0.95),tol=1e-10)$root
  cat('Two-active/three-active frequency crossover at B=50, phi=1: beta=',
      format(crossover,digits=10),'\n')
}
