# Shared numerical engine for 05_explore_asymmetric_tournaments.R.
# Extracted from Script 05 so exploratory plots use the same global-boundary
# and population calculations, without running its four-scenario figure.
source("05_linear_population_means.R")
new_asymmetric_model <- function(nquad=31L) {
  model <- new.env(parent=globalenv())
  evalq({
g<-function(x,beta){
  (1+x)^beta-1
}

w<-function(x,alpha,beta,phi){
  (1+alpha*g(x,beta))^phi-1
}

wprime<-function(x,alpha,beta,phi){
  phi*alpha*beta*(1+x)^(beta-1)*
    (1+alpha*g(x,beta))^(phi-1)
}

#------------------------------------------------------------
#Exact global maxxer boundaries
#------------------------------------------------------------

#For fixed alpha2 and candidate allocation x>0, this gives the alpha1
#at which x ties complete 2maxxing:
#
# s*w1(x,alpha1)+w2(B-x,alpha2)=w2(B,alpha2),
#
#where s=kappa1/kappa2.
#
#At x=0 the limiting value is the endpoint-slope boundary.
tie_alpha1<-function(x,alpha2,B,beta1,beta2,phi1,phi2,s){
  if(x<=0){
    return(
      phi2*alpha2*beta2*(1+B)^(beta2-1)*
        (1+alpha2*g(B,beta2))^(phi2-1)/
        (s*phi1*beta1)
    )
  }

  d<-w(B,alpha2,beta2,phi2)-w(B-x,alpha2,beta2,phi2)
  ((1+d/s)^(1/phi1)-1)/g(x,beta1)
}

#Complete 2maxxing is globally optimal below this alpha1 threshold.
#The threshold is the lowest talent at which ANY alternative allocation
#can tie the 2maxxer.
critical_alpha1<-function(alpha2,B,beta1,beta2,phi1,phi2,s,nroot=501){
  if(alpha2==0) return(0)

  #If both tournament-return functions are concave, the transition is continuous.
  if(phi1<=1 && phi2<=1){
    return(tie_alpha1(0,alpha2,B,beta1,beta2,phi1,phi2,s))
  }

  xx<-seq(0,B,length.out=nroot)
  yy<-vapply(xx,tie_alpha1,1,
             alpha2=alpha2,B=B,beta1=beta1,beta2=beta2,
             phi1=phi1,phi2=phi2,s=s)

  #Always include both endpoints as candidates.
  cand<-c(yy[1],yy[length(yy)])

  #Refine every interior local minimum of the tie threshold.
  k<-which(yy[2:(length(yy)-1)]<=yy[1:(length(yy)-2)] &
             yy[2:(length(yy)-1)]<=yy[3:length(yy)])+1

  if(length(k)){
    cand<-c(cand,vapply(k,function(i){
      optimize(function(x){
        tie_alpha1(x,alpha2,B,beta1,beta2,phi1,phi2,s)
      },c(xx[i-1],xx[i+1]),tol=1e-10)$objective
    },1))
  }

  min(cand,na.rm=TRUE)
}

#b2(alpha2): alpha1 threshold below which 2maxxing is globally optimal.
#b1(alpha1): alpha2 threshold below which 1maxxing is globally optimal.
boundary_curves<-function(alpha,B,beta1,beta2,phi1,phi2,s,nroot=501){
  b2<-vapply(alpha,critical_alpha1,1,
             B=B,beta1=beta1,beta2=beta2,
             phi1=phi1,phi2=phi2,s=s,nroot=nroot)

  #The 1maxxer boundary follows by swapping tournament labels.
  b1<-vapply(alpha,critical_alpha1,1,
             B=B,beta1=beta2,beta2=beta1,
             phi1=phi2,phi2=phi1,s=1/s,nroot=nroot)

  list(b2=b2,b1=b1)
}

#------------------------------------------------------------
#Optimal allocation in talent space
#------------------------------------------------------------

strategy_map<-function(alpha,B,beta1,beta2,phi1,phi2,s,nx=801,nroot=501){
  n<-length(alpha)
  bc<-boundary_curves(alpha,B,beta1,beta2,phi1,phi2,s,nroot)

  z<-matrix(NA_real_,n,n)
  max2<-matrix(FALSE,n,n)
  max1<-matrix(FALSE,n,n)

  xx<-seq(0,B,length.out=nx)
  w1<-outer(alpha,xx,function(a,x){w(x,a,beta1,phi1)})
  w2<-outer(alpha,xx,function(a,x){w(x,a,beta2,phi2)})
  w2r<-w2[,nx:1,drop=FALSE]

  tol<-1e-10

  for(j in seq_len(n)){
    is2<-alpha<=bc$b2[j]*(1+tol)
    is1<-alpha[j]<=bc$b1*(1+tol)
    both<-is1 & is2

    if(any(is2 & !both)) max2[is2 & !both,j]<-TRUE
    if(any(is1 & !both)) max1[is1 & !both,j]<-TRUE

    #If numerical tolerances make both endpoint criteria true, choose the
    #endpoint with the larger actual utility. Exact ties form only a boundary.
    if(any(both)){
      ii<-which(both)
      u2<-rep(w(B,alpha[j],beta2,phi2),length(ii))
      u1<-s*w(B,alpha[ii],beta1,phi1)
      choose1<-u1>u2
      max1[ii[choose1],j]<-TRUE
      max2[ii[!choose1],j]<-TRUE
    }

    z[max2[,j],j]<-0
    z[max1[,j],j]<-1

    #Only the region not analytically classified as a maxxer is searched
    #numerically. Endpoints are excluded because they are known not to be global.
    inside<-which(!max1[,j] & !max2[,j])

    if(length(inside)){
      uu<-s*w1[inside,2:(nx-1),drop=FALSE]
      uu<-sweep(uu,2,w2r[j,2:(nx-1)],"+")
      k<-max.col(uu,ties.method="first")+1
      z[inside,j]<-xx[k]/B
    }
  }

  list(z=z,b2=bc$b2,b1=bc$b1,max2=max2,max1=max1)
}

#------------------------------------------------------------
#Population expectations and responsive kappa
#------------------------------------------------------------

gauss_hermite<-function(n=31){
  J<-matrix(0,n,n)
  off<-sqrt((1:(n-1))/2)
  J[cbind(1:(n-1),2:n)]<-off
  J[cbind(2:n,1:(n-1))]<-off

  e<-eigen(J,symmetric=TRUE)
  o<-order(e$values)

  list(z=e$values[o],w=sqrt(pi)*e$vectors[1,o]^2)
}


population_cache<-function(beta1,beta2,phi1,phi2,B,nx=801){
  xx<-seq(0,B,length.out=nx)

  w1<-outer(aq,xx,function(a,x){w(x,a,beta1,phi1)})
  w2<-outer(aq,xx,function(a,x){w(x,a,beta2,phi2)})

  list(xx=xx,w1=w1,w2r=w2[,nx:1,drop=FALSE],
       beta1=beta1,beta2=beta2,phi1=phi1,phi2=phi2)
}

#For a trial s=kappa1/kappa2, calculate population-optimal allocations,
#E1, E2, M1, M2, and the resulting updated kappa ratio.
population_state<-function(s,cache,B){
  E1<-0
  M1<-0
  M2<-0

  for(i in seq_along(aq)){
    uu<-sweep(cache$w2r,2,s*cache$w1[i,],"+")
    k<-max.col(uu,ties.method="first")

    pairw<-wq[i]*wq
    x<-cache$xx[k]

    E1<-E1+sum(pairw*x)
    M1<-M1+sum(pairw*cache$w1[i,k])
    M2<-M2+sum(pairw*cache$w2r[cbind(seq_along(aq),k)])
  }

  E2<-B-E1

  if(E1<=0 || E2<=0 || M1<=0 || M2<=0){
    return(list(E1=E1,E2=E2,M1=M1,M2=M2,snew=NA_real_))
  }

  list(E1=E1,E2=E2,M1=M1,M2=M2,
       snew=(E1/M1)/(E2/M2))
}

#Solve log(G(s))-log(s)=0. Stop rather than silently selecting one
#if the chosen scenario produces multiple responsive-allocation equilibria.
equilibrium_s<-function(cache,B){
  if(cache$phi1==1 && cache$phi2==1)
    return(linear_population_equilibrium(cache$beta1,cache$beta2,B,m,slog))
  f<-function(logs){
    st<-population_state(exp(logs),cache,B)
    if(!is.finite(st$snew) || st$snew<=0) return(NA_real_)
    log(st$snew)-logs
  }

  zz<-seq(-8,8,length.out=121)
  ff<-vapply(zz,f,1)

  k<-which(is.finite(ff[-length(ff)]) &
             is.finite(ff[-1]) &
             ff[-length(ff)]*ff[-1]<=0)

  if(length(k)!=1){
    stop("Expected one responsive-allocation equilibrium; found ",length(k),".")
  }

  s<-exp(uniroot(f,c(zz[k],zz[k+1]),tol=1e-9)$root)
  st<-population_state(s,cache,B)

  c(s=s,E1=st$E1,E2=st$E2,M1=st$M1,M2=st$M2)
}

#------------------------------------------------------------
#Exact 50:50 stationary curve
#------------------------------------------------------------

#For each alpha2, solve the marginal-return equality at x=B/2.
#The curve is shown only where this stationary allocation is globally optimal.
half_curve<-function(alpha2,B,beta1,beta2,phi1,phi2,s,alim,nx=801){
  mid<-B/2
  xx<-seq(0,B,length.out=nx)

  f<-function(alpha1,alpha2){
    s*wprime(mid,alpha1,beta1,phi1)-
      wprime(mid,alpha2,beta2,phi2)
  }

  alpha1<-rep(NA_real_,length(alpha2))
  global<-logical(length(alpha2))

  for(j in seq_along(alpha2)){
    a2<-alpha2[j]

    if(a2==0){
      alpha1[j]<-0
      global[j]<-TRUE
      next
    }

    lo<-f(alim[1],a2)
    hi<-f(alim[2],a2)

    if(!is.finite(lo) || !is.finite(hi) || lo*hi>0) next

    alpha1[j]<-uniroot(function(a1){
      f(a1,a2)
    },alim,tol=1e-10)$root

    u<-s*w(xx,alpha1[j],beta1,phi1)+
      w(B-xx,a2,beta2,phi2)

    kmid<-which.min(abs(xx-mid))
    global[j]<-abs(u[kmid]-max(u))<
      1e-8*max(1,abs(max(u)))
  }

  list(x=alpha1,y=alpha2,global=global)
}


talent_hdr_curves<-function(median=1,sd=0.5,rho=0,probs=c(0.50,0.90,0.99)){
  m<-log(median)
  s<-sqrt(log((1+sqrt(1+4*(sd/median)^2))/2))

  mu<-c(m,m)
  Sigma<-s^2*matrix(c(1,rho,rho,1),2,2)

  one<-c(1,1)
  centre<-mu-Sigma%*%one
  ncp<-drop(t(one)%*%Sigma%*%one)

  qs<-qchisq(probs,df=2,ncp=ncp)

  th<-seq(0,2*pi,length.out=1001)
  L<-t(chol(Sigma))

  lapply(qs,function(q){
    y<-drop(centre)+L%*%rbind(sqrt(q)*cos(th),sqrt(q)*sin(th))
    rbind(exp(y[1,]),exp(y[2,]))
  })
}

  },model)
  gh <- model$gauss_hermite(nquad)
  model$m <- 0
  model$slog <- sqrt(log((1+sqrt(2))/2))
  model$aq <- exp(sqrt(2)*sqrt(log((1+sqrt(2))/2))*gh$z)
  model$wq <- gh$w/sqrt(pi)
  model
}
