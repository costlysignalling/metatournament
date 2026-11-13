source('analytic_participation.R')

#Deterministic integration of the exact strategy regions for 0<phi<=1.
#The three-active probability conditions on the largest two talents a>=b.
#The third talent enters iff its initial marginal return exceeds the common
#marginal return of the optimal two-tournament allocation.
participation_concave<-function(beta,phi,B=50,slog=talent_log_sd(),
                                m=0,force_nested=FALSE){
  stopifnot(beta>0,beta<1,phi>0,phi<=1,B>0,slog>0)
  if(phi==1 && !force_nested) return(participation_phi1(beta,B,3,slog))
  gb<-g(B,beta)
  log_boundary<-function(loga){
    loga+(beta-1)*log1p(B)+(phi-1)*log1p(exp(loga)*gb)
  }
  p1<-3*integrate(function(z){
    dnorm(z)*pnorm((log_boundary(m+slog*z)-m)/slog)^2
  },-9,9,rel.tol=1e-9,abs.tol=1e-12)$value

  inner<-function(zb){
    logb<-m+slog*zb
    b<-exp(logb)
    logupper<-logb+2
    while(log_boundary(logupper)<logb) logupper<-logupper+2
    logamax<-uniroot(function(t) log_boundary(t)-logb,
                    c(logb,logupper),tol=1e-11)$root
    za_max<-min(9,(logamax-m)/slog)
    if(za_max<=zb) return(0)
    integrand<-function(za){
      vapply(za,function(zz){
        a<-exp(m+slog*zz)
        log_derivative<-function(x,aa){
          log(aa)+(beta-1)*log1p(x)+(phi-1)*log1p(aa*g(x,beta))
        }
        delta<-function(x) log_derivative(x,a)-log_derivative(B-x,b)
        if(delta(B)>=0){
          c<-b
        } else {
          x<-uniroot(delta,c(B/2,B),tol=1e-10)$root
          c<-exp(log_derivative(x,a))
        }
        dnorm(zz)*max(0,pnorm(zb)-pnorm((log(c)-m)/slog))
      },1)
    }
    integrate(integrand,zb,za_max,rel.tol=1e-7,abs.tol=1e-10)$value
  }
  p3<-6*integrate(function(zb){dnorm(zb)*vapply(zb,inner,1)},
                  -9,9,rel.tol=1e-7,abs.tol=1e-9)$value
  ans<-c(one=p1,two=1-p1-p3,three=p3)
  stopifnot(all(ans>=-1e-7),abs(sum(ans)-1)<1e-12)
  pmax(ans,0)
}

g<-function(x,beta) expm1(beta*log1p(x))

if(sys.nframe()==0L){
  for(beta in c(0.4,0.8)){
    nested<-participation_concave(beta,1,force_nested=TRUE)
    direct<-participation_phi1(beta,50)
    stopifnot(max(abs(nested-direct))<2e-7)
    cat('General region integral agrees with phi=1 closed reduction:',beta,'\n')
  }
  out<-list()
  for(beta in c(0.4,0.8,0.9)) for(phi in c(0.5,0.66,1)){
    for(B in c(50,75)){
      p<-participation_concave(beta,phi,B)
      out[[length(out)+1]]<-data.frame(beta,phi,B,
        one=unname(p[1]),two=unname(p[2]),three=unname(p[3]))
      cat('beta=',beta,' phi=',phi,' B=',B,': ',paste(round(p,8),collapse=' '),'\n')
    }
  }
  write.csv(do.call(rbind,out),'participation_concave.csv',row.names=FALSE)
}
