source('10_collective_convincing.R')
comparison_price_compile<-local({env<-NULL;function(){
 if(is.null(env)){env<<-new.env();Rcpp::sourceCpp('10_comparison_core.cpp',env=env,
  cacheDir=file.path(tempdir(),'10_comparison_core'))};env
}})
comparison_gh<-function(n){
 J<-matrix(0,n,n);J[cbind(1:(n-1),2:n)]<-sqrt(1:(n-1));J<-J+t(J)
 e<-eigen(J,symmetric=TRUE);list(node=e$values,weight=e$vectors[1,]^2)
}
comparison_gl<-function(n){
 J<-matrix(0,n,n);i<-1:(n-1);J[cbind(i,i+1)]<-i/sqrt(4*i*i-1);J<-J+t(J)
 e<-eigen(J,symmetric=TRUE);list(node=e$values,weight=2*e$vectors[1,]^2,split=TRUE)
}
comparison_collective_prepare<-function(beta,phi,B=50,efficiency=.1,order=10L,tolerance=1e-4,split_outer=FALSE){
 stopifnot(length(beta)==2,length(phi)==2,all(is.finite(c(beta,phi,B,efficiency))),
  all(beta>0&beta<1),all(phi>0),B>0,efficiency>=0,order>=2)
 list(beta=beta,phi=phi,B=B,efficiency=efficiency,
  gh=if(split_outer)comparison_gl(order)else comparison_gh(order),
  tolerance=tolerance,cache=new.env(parent=emptyenv()),
  analytic=all(phi==1)&&identical(beta[1],beta[2]))
}
comparison_collective_outcome<-function(p,v1=0,v2=0){
 key<-sprintf('%.7f_%.7f',v1,v2)
 if(exists(key,p$cache,inherits=FALSE))return(get(key,p$cache))
 analytic_result<-if(p$analytic)tryCatch(
  collective_outcome(collective_prepare(p$beta[1],p$B,p$efficiency),v1,v2),
  error=function(e)NULL)else NULL
 if(!is.null(analytic_result)){
  a<-analytic_result
  a$H<-a$log_ratio+a$delta;a$collapsed<-FALSE
 }else{
  R<-p$B-abs(c(v1,v2));delta<-p$efficiency*(v1+v2)/2
  evaluate<-function(H){
   m<-comparison_price_compile()$comparison_cohort_moments_cpp(p$gh$node,p$gh$weight,
    R[1],R[2],H,p$beta[1],p$beta[2],p$phi[1],p$phi[2],boredom_slog(),p$tolerance,
    nx=if(is.null(p$nx))32L else p$nx,split_outer=isTRUE(p$gh$split))
   E<-colSums(m[,1:2,drop=FALSE])+c(max(v1,0)+max(v2,0),max(-v1,0)+max(-v2,0))/2
   M<-colSums(m[,3:4,drop=FALSE])
   z<-log(pmax(E[1],1e-300))-log(pmax(E[2],1e-300))+
    log(pmax(M[2],1e-300))-log(pmax(M[1],1e-300))
   list(H=H,E=E,M=M,cohortM=m[,3:4,drop=FALSE],residual=H-delta-z)
  }
  symmetric_profile<-identical(p$beta[1],p$beta[2])&&identical(p$phi[1],p$phi[2])&&abs(v1+v2)<1e-10
  if(symmetric_profile){
   a<-evaluate(0);a$E[]<-p$B/2;a$M[]<-mean(a$M);a$payoff<-c(p$B,p$B)
   a$collapsed<-FALSE;a$residual<-0
  }else{
  interval<-delta+c(-1,1);lo<-evaluate(interval[1]);hi<-evaluate(interval[2])
  for(width in c(4,8,16))if(lo$residual*hi$residual>0){
   interval<-delta+c(-width,width);lo<-evaluate(interval[1]);hi<-evaluate(interval[2])
  }
  collapsed<-FALSE
  if(lo$residual*hi$residual>0){
   # Extremely small unvisited pools have a negligible MATERIAL payoff.
   # Only use their limiting payoff if no campaign donates to that pool.
   candidates<-list(lo,hi)
   for(test in candidates){
    minor<-which.min(test$E)
    donated<-if(minor==1)max(v1,0)+max(v2,0)else max(-v1,0)+max(-v2,0)
    if(test$E[minor]<1e-10&&donated==0){
     a<-test;a$E[minor]<-0;a$E[3-minor]<-p$B
     a$payoff<-2*a$E[3-minor]*a$cohortM[,3-minor]/a$M[3-minor]
     a$collapsed<-TRUE;a$residual<-NA_real_;collapsed<-TRUE;break
    }
   }
   if(!collapsed)stop('No bracketed population response for campaigns ',v1,', ',v2)
  }else{
   root<-uniroot(function(H)evaluate(H)$residual,interval,
    f.lower=lo$residual,f.upper=hi$residual,tol=max(1e-7,p$tolerance/10))$root
   a<-evaluate(root)
   if(abs(a$residual)>max(2e-4,10*p$tolerance)&&min(a$E)>1e-8)
    stop('Collective population residual too large: ',a$residual,
     '; campaigns ',v1,', ',v2,'; E ',paste(signif(a$E,5),collapse=', '))
   a$payoff<-2*as.vector(a$cohortM%*%(a$E/a$M));a$collapsed<-FALSE
  }
  }
  a$v<-c(v1,v2);a$delta<-delta;a$log_ratio<-a$H-delta
 }
 stopifnot(abs(sum(a$E)-p$B)<1e-6,abs(sum(a$payoff)-2*p$B)<1e-6)
 assign(key,a,p$cache);a
}
comparison_collective_br<-function(p,other=0,probability=1,group=1,grid=13L){
 value<-function(v){
  sum(vapply(seq_along(other),function(k){
   a<-if(group==1)comparison_collective_outcome(p,v,other[k])else
    comparison_collective_outcome(p,other[k],v)
   a$payoff[group]*probability[k]
  },0.))
 }
 if((length(p$gh$node)>=16||isTRUE(p$gh$split))&&!p$analytic){
  # Search the whole interval at a broader integration resolution, then
  # refine every competitive peak using the fine population calculation.
  # This avoids spending fine quadrature on grossly dominated campaigns
  # that leave almost no training time and an astronomically small pool.
  coarse<-p;coarse$gh<-comparison_gh(8);coarse$tolerance<-2e-4
  coarse$cache<-new.env(parent=emptyenv())
  first<-comparison_collective_br(coarse,other,probability,group,grid)
  peaks<-attr(first,'peaks')
  peaks<-peaks[peaks[,2]>=max(peaks[,2])-.5,,drop=FALSE]
  candidates<-matrix(numeric(),0,2,dimnames=list(NULL,c('v','payoff')))
  for(j in seq_len(nrow(peaks))){
   center<-peaks[j,1];width<-2*p$B/(grid-1)
   interval<-pmax(-.995*p$B,pmin(.995*p$B,center+c(-width,width)))
   fit<-optimize(value,interval,maximum=TRUE,tol=max(1e-4,20*p$tolerance))
   candidates<-rbind(candidates,c(fit$maximum,fit$objective),c(center,value(center)))
  }
  ans<-candidates[which.max(candidates[,2]),]
  attr(ans,'peaks')<-candidates;attr(ans,'coarse_scan')<-attr(first,'peaks')
  return(ans)
 }
 actions<-seq(-.995*p$B,.995*p$B,length.out=grid)
 values<-vapply(actions,value,0.)
 best<-which(values>=c(-Inf,head(values,-1))&values>=c(tail(values,-1),-Inf))
 candidates<-cbind(v=actions,payoff=values)
 peaks<-cbind(v=actions[c(1,length(actions))],payoff=values[c(1,length(actions))])
 for(j in best){
  fit<-optimize(value,actions[c(max(1,j-1),min(grid,j+1))],maximum=TRUE,
   tol=max(1e-4,20*p$tolerance))
  candidates<-rbind(candidates,c(fit$maximum,fit$objective))
  peaks<-rbind(peaks,c(fit$maximum,fit$objective))
 }
 ans<-candidates[which.max(candidates[,2]),]
 attr(ans,'peaks')<-peaks;ans
}
comparison_collective_pure<-function(p,start=c(0,0),maxit=18L){
 v<-start;history<-matrix(v,1,2)
 symmetric<-identical(p$beta[1],p$beta[2])&&identical(p$phi[1],p$phi[2])
 if(isTRUE(p$gh$split)&&!p$analytic)return(comparison_split_equilibrium(p,start))
 if(!symmetric&&(length(p$gh$node)>=16||isTRUE(p$gh$split)))return(comparison_collective_broyden(p,start))
 if(symmetric){
  # A symmetric saddle has v2=-v1. Find stationary candidates cheaply,
  # then verify against GLOBAL deviations in both campaign directions.
  h<-if(p$tolerance>1e-5).15 else .04
  derivative<-function(y){
   (comparison_collective_outcome(p,y+h,-y)$payoff[1]-
    comparison_collective_outcome(p,y-h,-y)$payoff[1])/(2*h)
  }
  yy<-sort(unique(c(seq(-.9*p$B,.9*p$B,length.out=9),start[1])))
  gg<-vapply(yy,derivative,0.);roots<-numeric()
  for(j in 1:(length(yy)-1))if(gg[j]*gg[j+1]<0){
   r<-tryCatch(uniroot(derivative,yy[c(j,j+1)],tol=.002)$root,error=function(e)NA_real_)
   if(is.finite(r))roots<-c(roots,r)
  }
  candidates<-unique(c(start[1],0,roots[roots>=0],roots[roots<0]))
  for(y in candidates){
   a<-comparison_collective_outcome(p,y,-y)
   br<-comparison_collective_br(p,-y,group=1,grid=17)
   regret<-unname(br['payoff'])-p$B
   history<-rbind(history,c(y,-y),c(unname(br['v']),-unname(br['v'])))
   message('Symmetric candidate ',round(y,3),'; regret ',signif(regret,3))
   if(regret<max(.003,100*p$tolerance))return(list(v=c(y,-y),outcome=a,
    regret=c(regret,regret),history=history,pure=TRUE))
  }
 }
 for(it in seq_len(maxit)){
  br1<-comparison_collective_br(p,v[2],group=1)
  br2<-if(symmetric)c(v=-br1['v'],payoff=br1['payoff'])else
   comparison_collective_br(p,v[1],group=2)
  target<-c(unname(br1['v']),unname(br2[1]))
  gap<-max(abs(target-v))
  message('Collective iteration ',it,': ',paste(round(v,3),collapse=', '),' ; gap ',signif(gap,3))
  history<-rbind(history,target)
  if(gap<max(.01,100*p$tolerance))break
  v<-(v+target)/2
 }
 a<-comparison_collective_outcome(p,v[1],v[2])
 # Evaluate deviations against the FINAL opponent policy, including when
 # the iteration limit was reached immediately after a damped update.
 br1<-comparison_collective_br(p,v[2],group=1)
 br2<-if(symmetric)br1 else comparison_collective_br(p,v[1],group=2)
 regret<-c(br1['payoff'],if(symmetric)br1['payoff']else br2['payoff'])-a$payoff
 list(v=v,outcome=a,regret=as.numeric(regret),history=history,
  pure=max(regret)<max(.003,100*p$tolerance))
}
comparison_split_deviations<-function(p,v){
 state<-comparison_collective_outcome(p,v[1],v[2]);best<-v;gain<-numeric(2)
 coarse<-p;coarse$gh<-comparison_gh(8);coarse$tolerance<-2e-4
 coarse$cache<-new.env(parent=emptyenv())
 for(group in 1:2){
  br<-comparison_collective_br(coarse,v[3-group],group=group,grid=17)
  peaks<-attr(br,'peaks');peaks<-peaks[peaks[,2]>=max(peaks[,2])-.5,,drop=FALSE]
  candidates<-unique(c(v[group],0,peaks[,1]));payoff<-numeric(length(candidates))
  for(j in seq_along(candidates)){
   trial<-v;trial[group]<-candidates[j]
   payoff[j]<-comparison_collective_outcome(p,trial[1],trial[2])$payoff[group]
  }
  best[group]<-candidates[which.max(payoff)];gain[group]<-max(payoff)-state$payoff[group]
 }
 list(target=best,regret=gain,outcome=state)
}
comparison_split_equilibrium<-function(p,start){
 v<-start;history<-matrix(v,1,2)
 symmetric<-identical(p$beta[1],p$beta[2])&&identical(p$phi[1],p$phi[2])
 F<-function(x)comparison_collective_outcome(p,x[1],x[2])$payoff[1]
 for(phase in 1:2){
  if(symmetric){
   y<-v[1];previous_y<-previous_g<-NULL;h<-.25;slope<--.03
   for(it in 1:5){
    g<-(F(c(y+h,-y))-F(c(y-h,-y)))/(2*h)
    message('Split symmetric ',it,': effort ',round(y,4),'; gradient ',signif(g,3))
    if(abs(g)<.001)break
    if(!is.null(previous_y)&&abs(y-previous_y)>1e-6)slope<-(g-previous_g)/(y-previous_y)
    if(slope>=0||!is.finite(slope))slope<--.03
    step<-max(-1,min(1,-g/slope));previous_y<-y;previous_g<-g;y<-y+step
   }
   v<-c(y,-y)
  }else{
   h<-.15
   for(it in 1:4){
    if(any(abs(v)<2*h))break
    center<-F(v);xp<-F(v+c(h,0));xm<-F(v-c(h,0));yp<-F(v+c(0,h));ym<-F(v-c(0,h))
    g<-c((xp-xm)/(2*h),-(yp-ym)/(2*h))
    xx<-(xp-2*center+xm)/h^2;yy<-(yp-2*center+ym)/h^2
    xy<-(F(v+c(h,h))-F(v+c(h,-h))-F(v+c(-h,h))+F(v-c(h,h)))/(4*h*h)
    J<-matrix(c(xx,-xy,xy,-yy),2,2)
    step<-tryCatch(as.vector(solve(J,-g)),error=function(e)c(0,0))
    if(max(abs(step))>2)step<-step*2/max(abs(step))
    message('Split saddle ',it,': ',paste(round(v,4),collapse=', '),
     '; gradients ',paste(signif(g,3),collapse=', '))
    if(max(abs(g))<.001||max(abs(step))<.01)break
    v<-pmax(-.99*p$B,pmin(.99*p$B,v+step));history<-rbind(history,v)
   }
  }
  check<-comparison_split_deviations(p,v);history<-rbind(history,check$target)
  message('Split deviation gains: ',paste(signif(check$regret,4),collapse=', '))
  if(max(check$regret)<.003)return(list(v=v,outcome=check$outcome,regret=check$regret,
   history=history,pure=TRUE))
  v<-(v+check$target)/2
 }
 list(v=v,outcome=comparison_collective_outcome(p,v[1],v[2]),regret=check$regret,history=history,pure=FALSE)
}
comparison_collective_broyden<-function(p,start){
 v<-start;history<-matrix(v,1,2);J<-matrix(c(-1,-1,1,-1),2,2)
 previous_v<-previous_gap<-NULL
 for(it in 1:10){
  a<-comparison_collective_outcome(p,v[1],v[2])
  b1<-comparison_collective_br(p,v[2],group=1,grid=17)
  b2<-comparison_collective_br(p,v[1],group=2,grid=17)
  target<-c(unname(b1['v']),unname(b2['v']));gap<-target-v
  regret<-c(unname(b1['payoff']),unname(b2['payoff']))-a$payoff
  history<-rbind(history,target)
  message('Coupled response ',it,': ',paste(round(v,4),collapse=', '),
   '; target ',paste(round(target,4),collapse=', '),
   '; gains ',paste(signif(regret,4),collapse=', '))
  if(max(regret)<.003)return(list(v=v,outcome=a,regret=regret,history=history,pure=TRUE))
  if(!is.null(previous_v)){
   s<-v-previous_v
   if(sum(s*s)>1e-10)J<-J+as.vector(gap-previous_gap-J%*%s)%o%s/sum(s*s)
  }
  step<-tryCatch(as.vector(solve(J,-gap)),error=function(e)gap/3)
  if(max(abs(step))>2)step<-step*2/max(abs(step))
  previous_v<-v;previous_gap<-gap
  v<-pmax(-.99*p$B,pmin(.99*p$B,v+step));history<-rbind(history,v)
 }
 list(v=previous_v,outcome=a,regret=regret,history=history,pure=FALSE)
}
# Linear-programming barrier for a finite zero-sum payoff matrix. No packages.
comparison_zero_sum<-function(A){
 solve_min<-function(M){
  n<-ncol(M)
  if(n==1)return(list(probability=1,value=max(M)))
  ui<-rbind(cbind(-sweep(M[,-n,drop=FALSE],1,M[,n],'-'),1),
   cbind(diag(n-1),0),c(rep(-1,n-1),0))
  ci<-c(M[,n],rep(0,n-1),-1)
  theta<-c(rep(1/n,n-1),max(rowMeans(M))+1)
  fit<-constrOptim(theta,function(z)z[n],function(z)c(rep(0,n-1),1),ui,ci,
   mu=1e-7,control=list(reltol=1e-11,maxit=2000),outer.iterations=40)
  q<-c(fit$par[1:(n-1)],1-sum(fit$par[1:(n-1)]));q<-pmax(q,0);q<-q/sum(q)
  list(probability=q,value=max(as.vector(M%*%q)))
 }
 q<-solve_min(A);p<-solve_min(-t(A))
 list(p=p$probability,q=q$probability,gap=q$value+p$value)
}
comparison_collective_solve<-function(beta,phi,B=50,efficiency=.1,start=NULL,campaign_method=c('auto','mixed')){
 campaign_method<-match.arg(campaign_method)
 coarse<-comparison_collective_prepare(beta,phi,B,efficiency,order=8,tolerance=2e-4)
 first<-if(is.null(start))comparison_collective_pure(coarse,maxit=6)else list(v=start)
 fine<-comparison_collective_prepare(beta,phi,B,efficiency,order=8,tolerance=1e-5,split_outer=TRUE)
 message('Refining collective equilibrium')
 second<-if(campaign_method=='mixed')list(v=first$v,pure=FALSE,
  history=rbind(first$v-.3,first$v+.3))else comparison_collective_pure(fine,first$v,maxit=14)
 if(second$pure){
  return(list(beta=beta,phi=phi,B=B,efficiency=efficiency,pure=TRUE,
   v1=second$v[1],v2=second$v[2],p1=1,p2=1,outcomes=list(second$outcome),
   regret=second$regret,coarse_v=first$v,order=8,tolerance=1e-5,integration='split-legendre'))
 }
 message('Pure fixed point not verified; solving a mixed campaign game')
 action1<-sort(unique(c(second$v[1],tail(second$history[,1],2))))
 action2<-sort(unique(c(second$v[2],tail(second$history[,2],2))))
 for(it in 1:12){
  A<-outer(action1,action2,Vectorize(function(a,b)comparison_collective_outcome(fine,a,b)$payoff[1]))
  mix<-comparison_zero_sum(A-B)
  s1<-which(mix$p>1e-5);s2<-which(mix$q>1e-5)
  p1<-mix$p[s1]/sum(mix$p[s1]);p2<-mix$q[s2]/sum(mix$q[s2])
  br1<-comparison_collective_br(fine,action2[s2],p2,1,grid=25)
  br2<-comparison_collective_br(fine,action1[s1],p1,2,grid=25)
  value<-as.numeric(t(p1)%*%A[s1,s2,drop=FALSE]%*%p2)
  regret<-c(br1['payoff']-value,br2['payoff']-(2*B-value))
  message('Campaign game ',it,': gains ',paste(signif(regret,4),collapse=', '),
   '; supports ',length(s1),', ',length(s2))
  if(max(regret)<.01){
   # Explicit flat order: campaign 2 varies fastest.
   outcomes<-lapply(seq_along(action1[s1]),function(i)lapply(seq_along(action2[s2]),
    function(j)comparison_collective_outcome(fine,action1[s1][i],action2[s2][j])))
   outcomes<-do.call(c,outcomes)
   return(list(beta=beta,phi=phi,B=B,efficiency=efficiency,pure=length(s1)==1&&length(s2)==1,
    v1=action1[s1],v2=action2[s2],p1=p1,p2=p2,outcomes=outcomes,
    regret=as.numeric(regret),order=8,tolerance=1e-5,integration='split-legendre'))
  }
  action1<-sort(unique(c(action1,unname(br1['v']))));action2<-sort(unique(c(action2,unname(br2['v']))))
 }
 stop('Mixed collective campaign regret remains above tolerance.')
}

