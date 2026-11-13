# Vectorized R implementation of the same finite-mass best-response problem.
# Used when Rcpp cannot be loaded, including Windows Application Control blocks.
boredom_R_core<-function(){
 core<-new.env(parent=globalenv());cache<-new.env(parent=emptyenv())
 weights<-function(x,a,b,p)expm1(p*log1p(a*expm1(b*log1p(x))))
 prepare<-function(actor,beta1,beta2,phi1,phi2,B,compulsory,grid_half){
  key<-paste(beta1,beta2,phi1,phi2,B,compulsory,grid_half,nrow(actor),sum(actor),sum(actor^2),collapse=':')
  if(exists(key,cache,inherits=FALSE)){old<-get(key,cache);if(identical(old$actor,actor))return(old)}
  g<-c(seq(0,compulsory*B,length.out=grid_half+1),seq(compulsory*B,B,length.out=grid_half+1)[-1]);free<-(1-compulsory)*B
  obj<-list(g=g,w1=outer(actor[,1],g,function(a,x)weights(x,a,beta1,phi1)),w2=outer(actor[,2],pmin(free,B-g),function(a,x)weights(x,a,beta2,phi2)),actor=actor,b1=beta1,b2=beta2,p1=phi1,p2=phi2,B=B,c=compulsory)
  # Keep only the current table; this avoids accumulating parameter-grid data.
  rm(list=ls(cache),envir=cache);assign(key,obj,cache);obj
 }
 best<-function(obj,bg){
  n<-nrow(obj$actor);g<-obj$g;ng<-length(g);B<-obj$B;c<-obj$c;free<-(1-c)*B
  u<-sweep(obj$w1,2,bg[1]+g,'*')/(bg[3]+obj$w1)+sweep(obj$w2,2,bg[2]+B-g,'*')/(bg[4]+obj$w2)
  ix<-max.col(u,ties.method='first');q<-g[ix];bv<-u[cbind(seq_len(n),ix)]
  peaks<-which(u[,2:(ng-1),drop=FALSE]>=u[,1:(ng-2),drop=FALSE]&u[,2:(ng-1),drop=FALSE]>=u[,3:ng,drop=FALSE],arr.ind=TRUE)
  mid<-(ng+1L)%/%2L
  derivative<-function(j,right){
   q<-g[j];u<-obj$w1[,j];v<-obj$w2[,j];s1<-u/(bg[3]+u);s2<-v/(bg[4]+v)
   wp1<-obj$p1*obj$actor[,1]*obj$b1*(1+q)^(obj$b1-1)*(1+u)^((obj$p1-1)/obj$p1)
   wp2<-obj$p2*obj$actor[,2]*obj$b2*(1+min(free,B-q))^(obj$b2-1)*(1+v)^((obj$p2-1)/obj$p2)
   s1+(bg[1]+q)/(bg[3]+u)*(1-s1)*wp1-s2-if(right)(bg[2]+B-q)/(bg[4]+v)*(1-s2)*wp2 else 0
  }
  extra<-list(which(derivative(1,FALSE)>0 & derivative(2,FALSE)<0),which(derivative(ng,TRUE)<0 & derivative(ng-1,TRUE)>0),which(derivative(mid,FALSE)<0 & derivative(mid-1,FALSE)>0),which(derivative(mid,TRUE)>0 & derivative(mid+1,TRUE)<0))
  ids<-c(peaks[,1],unlist(extra));left<-c(g[peaks[,2]],rep(c(g[1],g[ng-1],g[mid-1],g[mid]),lengths(extra)));right<-c(g[peaks[,2]+2],rep(c(g[2],g[ng],g[mid],g[mid+1]),lengths(extra)))
  if(!length(ids))return(list(q=q,value=bv))
  a1<-obj$actor[ids,1];a2<-obj$actor[ids,2]
  value<-function(x){w1<-weights(x,a1,obj$b1,obj$p1);w2<-weights(pmin(free,B-x),a2,obj$b2,obj$p2);(bg[1]+x)*w1/(bg[3]+w1)+(bg[2]+B-x)*w2/(bg[4]+w2)}
  ratio<-(sqrt(5)-1)/2;x<-right-ratio*(right-left);y<-left+ratio*(right-left);vx<-value(x);vy<-value(y)
  for(k in 1:38){
   keep<-vx>vy;ox<-x;oy<-y;ovx<-vx;ovy<-vy
   right[keep]<-oy[keep];left[!keep]<-ox[!keep]
   xx<-ifelse(keep,right-ratio*(right-left),left+ratio*(right-left));nv<-value(xx)
   x<-ifelse(keep,xx,oy);y<-ifelse(keep,ox,xx);vx<-ifelse(keep,nv,ovy);vy<-ifelse(keep,ovx,nv)
  }
  xx<-ifelse(vx>vy,x,y);vv<-pmax(vx,vy)
  # Several candidate maxima can belong to the same actor.
  ord<-order(ids,-vv);win<-ord[!duplicated(ids[ord])];rows<-ids[win];improve<-vv[win]>bv[rows]+1e-12
  q[rows[improve]]<-xx[win[improve]];bv[rows[improve]]<-vv[win[improve]]
  list(q=q,value=bv)
 }
 core$boredom_map_cpp<-function(actor,background,beta1,beta2,phi1,phi2,B=50,compulsory=.5,grid_half=160L){
  out<-numeric(nrow(actor))
  for(rows in split(seq_len(nrow(actor)),ceiling(seq_len(nrow(actor))/1024))){
   obj<-prepare(actor[rows,,drop=FALSE],beta1,beta2,phi1,phi2,B,compulsory,grid_half);out[rows]<-best(obj,background)$q
  }
  out
 }
 core$boredom_density_cpp<-function(talent,weights,moment,beta1,beta2,phi1,phi2,B=50,compulsory=.5,N=30,grid_half=80L,maxit=120L,tolerance=1e-7,responsive=TRUE){
  obj<-prepare(talent,beta1,beta2,phi1,phi2,B,compulsory,grid_half)
  evaluate<-function(m){q<-best(obj,(N-1)*c(m[1],B-m[1],m[2:3]))$q;list(q=q,mapped=c(sum(weights*q),sum(weights*expm1(phi1*log1p(talent[,1]*expm1(beta1*log1p(q))))),sum(weights*expm1(phi2*log1p(talent[,2]*expm1(beta2*log1p(pmin((1-compulsory)*B,B-q)))))),sum(weights*pmax(compulsory*B-q,0))))}
  encode<-function(m)c(qlogis(pmin(pmax(m[1]/B,1e-12),1-1e-12)),log(pmax(m[2:3],1e-20)))
  decode<-function(z)c(B*plogis(z[1]),exp(z[2:3]))
  residual<-function(z){m<-decode(z);g<-evaluate(m)$mapped;c((g[1]-m[1])/B,log(pmax(g[2:3],1e-20)/m[2:3]))}
  state<-moment;best_state<-state;error<-Inf;converged<-!responsive;iterations<-0L
  newton<-function(seed,limit){
   z<-encode(seed)
   for(it in seq_len(limit)){
    r<-residual(z);er<-max(abs(r));if(er<error){error<<-er;best_state<<-decode(z)};if(er<tolerance)return(TRUE)
    J<-vapply(1:3,function(j){zz<-z;zz[j]<-zz[j]+1e-4;(residual(zz)-r)/1e-4},numeric(3))
    step<-tryCatch(as.vector(solve(J,-r)),error=function(e)NULL);if(is.null(step)||any(!is.finite(step)))break
    accepted<-FALSE
    for(rate in 2^-(0:14)){zz<-z+rate*pmin(pmax(step,-2),2);if(max(abs(residual(zz)))<er){z<-zz;accepted<-TRUE;break}}
    if(!accepted)break
   }
   FALSE
  }
  if(responsive){
   # Warm population states usually need only a few Newton corrections.
   converged<-newton(moment,8L)
   if(!converged)for(damping in c(.7,.35,.15)){
    state<-best_state
    for(it in seq_len(maxit)){
     g<-evaluate(state)$mapped;r<-c((g[1]-state[1])/B,log(pmax(g[2:3],1e-20)/pmax(state[2:3],1e-20)));er<-max(abs(r))
     if(er<error){best_state<-state;error<-er;iterations<-it}
     if(er<tolerance){converged<-TRUE;break}
     state<-c(state[1]+damping*(g[1]-state[1]),exp((1-damping)*log(pmax(state[2:3],1e-20))+damping*log(pmax(g[2:3],1e-20))))
    }
    if(converged)break
   }
   if(!converged)converged<-newton(best_state,35L)
  }
  final<-evaluate(best_state)
  c(list(moment=best_state),final,list(converged=converged,residual=if(responsive)error else 0,iterations=iterations))
 }
 core
}
