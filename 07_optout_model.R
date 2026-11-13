# Shared model for the Script 07 opt-out figures. Run from the scripts folder.
# Requires Rcpp and a C++ compiler (Rtools on Windows).
# Compile in the local R temporary folder, avoiding Dropbox DLL-policy blocks.
if(!exists('.optout_cpp',envir=.GlobalEnv,inherits=FALSE)) {
  if(!requireNamespace('Rcpp',quietly=TRUE)) stop('Install Rcpp first.')
  core<-new.env(parent=globalenv())
  Rcpp::sourceCpp('07_optout_core.cpp',env=core,
    cacheDir=file.path(tempdir(),'07_optout_cpp'),rebuild=TRUE)
  assign('.optout_cpp',core,envir=.GlobalEnv)
}

optout_eta_tag <- function(eta) {
  if(length(eta)!=1L || !is.finite(eta) || eta<0 || eta>1)
    stop('Require a single eta between 0 and 1.')
  # Preserve distinct names for values with more than two decimal places.
  z<-format(eta,scientific=FALSE,trim=TRUE,digits=12)
  if(!grepl('.',z,fixed=TRUE)) z<-paste0(z,'.0')
  z
}
optout_file <- function(stem,eta,ext='png') paste0(stem,'_eta_',optout_eta_tag(eta),'.',ext)
optout_save <- function(draw,stem,eta,width=3640/300,height=1740/300) {
  if(isTRUE(getOption('metatournament.png',TRUE))) {
    png(optout_file(stem,eta),width=width,height=height,units='in',res=600)
    tryCatch(draw(),finally=dev.off())
  }
  if(isTRUE(getOption('metatournament.pdf',TRUE))) {
    pdf(optout_file(stem,eta,'pdf'),width=width,height=height)
    tryCatch(draw(),finally=dev.off())
  }
  if(interactive() && isTRUE(getOption('metatournament.plot',TRUE))) draw()
}

optout_slog <- function(median=1,sd=.5) sqrt(log((1+sqrt(1+4*(sd/median)^2))/2))
optout_w <- function(x,a,beta,phi) expm1(phi*log1p(a*expm1(beta*log1p(x))))
optout_gaussian <- function(n,normal=TRUE) {
  J<-matrix(0,n,n); k<-seq_len(n-1L)
  off<-if(normal) sqrt(k/2) else k/sqrt(4*k*k-1)
  J[cbind(k,k+1L)]<-off; J[cbind(k+1L,k)]<-off
  e<-eigen(J,symmetric=TRUE);o<-order(e$values)
  list(x=e$values[o],w=e$vectors[1,o]^2)
}

# Reference BEFORE opting out: X1~Uniform(0,B), X2=B-X1, independent of
# talents. E1=E2=B/2; M1 and M2 are calculated separately, not equated.
optout_fixed_scalers <- function(beta,phi,B=50,median=1,sd=.5,nquad=61L) {
  if(length(beta)==1) beta<-rep(beta,2)
  if(length(phi)==1) phi<-rep(phi,2)
  gh<-optout_gaussian(nquad);gl<-optout_gaussian(nquad,FALSE)
  a<-exp(log(median)+sqrt(2)*optout_slog(median,sd)*gh$x)
  x<-B*(gl$x+1)/2
  M<-vapply(1:2,function(t) {
    if(phi[t]==1) return(exp(log(median)+optout_slog(median,sd)^2/2)*
      (expm1((beta[t]+1)*log1p(B))/((beta[t]+1)*B)-1))
    sum(gh$w * (outer(a,x,function(a,x) optout_w(x,a,beta[t],phi[t])) %*% gl$w))
  },1)
  list(E=c(B/2,B/2),M=M,scaler=(B/2)/M)
}

optout_population <- function(n=getOption('metatournament.optout.npop',8192L),median=1,sd=.5) {
  n<-as.integer(n);if(n<32 || n%%2L) stop('Population size must be even and at least 32.')
  inverse<-function(i,base) {
    out<-numeric(length(i));f<-1/base
    while(any(i>0L)) {out<-out+f*(i%%base);i<-i%/%base;f<-f/base};out
  }
  i<-64L+seq_len(n/2L)
  a<-exp(log(median)+optout_slog(median,sd)*qnorm(inverse(i,2L)))
  b<-exp(log(median)+optout_slog(median,sd)*qnorm(inverse(i,3L)))
  # Exact exchange symmetry without imposing symmetry on tournament parameters.
  rbind(cbind(a,b),cbind(b,a))
}

optout_allocate <- function(talent,beta,phi,scaler,eta=.5,B=50,
    nscan=getOption('metatournament.optout.nscan',129L)) {
  talent<-as.matrix(talent)
  if(ncol(talent)!=2L || any(!is.finite(talent)) || any(talent<0) ||
     length(beta)!=2L || length(phi)!=2L || length(scaler)!=2L ||
     any(!is.finite(c(beta,phi,scaler,B))) || any(beta<=0|beta>=1) ||
     any(phi<=0) || any(scaler<0) || B<=0 || nscan<17)
    stop('Require two nonnegative talents, 0<beta<1, phi>0, scalers>=0 and B>0.')
  optout_eta_tag(eta)
  get('.optout_cpp',envir=.GlobalEnv)$optout_alloc_cpp(talent,beta,phi,scaler,B,eta,as.integer(nscan))
}

optout_state <- function(scaler,beta,phi,eta,B=50,population=optout_population()) {
  z<-optout_allocate(population,beta,phi,scaler,eta,B)
  E<-colMeans(z[,2:3,drop=FALSE])
  M<-vapply(1:2,function(t) mean(optout_w(z[,t+1],population[,t],beta[t],phi[t])),1)
  list(E0=mean(z[,1]),E1=unname(E[1]),E2=unname(E[2]),M1=M[1],M2=M[2],
       Utotal=mean(z[,4]),scaler=ifelse(M>0,E/M,NA_real_))
}

optout_talent_map <- function(alpha,beta,phi,scaler,eta,B=50,with_bounds=FALSE) {
  if(is.null(.GlobalEnv$.optout_cpp$maps)) .GlobalEnv$.optout_cpp$maps<-new.env(parent=emptyenv())
  key<-paste(c(beta,phi,scaler,eta,B,length(alpha),range(alpha),
    getOption('metatournament.optout.nscan',129L),with_bounds),collapse='_')
  cache<-.GlobalEnv$.optout_cpp$maps
  if(exists(key,cache,inherits=FALSE)) {
    previous<-get(key,cache)
    if(identical(previous$alpha,alpha)) return(previous$map)
  }
  n<-length(alpha);aa<-as.matrix(expand.grid(alpha,alpha))
  z<-optout_allocate(aa,beta,phi,scaler,eta,B)
  x0<-matrix(z[,1]/B,n,n);x1<-matrix(z[,2]/B,n,n);x2<-matrix(z[,3]/B,n,n)
  result<-list(u=matrix(z[,4],n,n),x0=x0,x1=x1,x2=x2,
    selfmax=x0>1-1e-9,max1=x1>1-1e-9,max2=x2>1-1e-9)
  if(with_bounds) result$bounds<-optout_maxxer_bounds(seq(min(alpha),max(alpha),length.out=1201L),
    beta,phi,scaler,eta,B)
  assign(key,list(alpha=alpha,map=result),cache);result
}

# Smooth global-deviation boundaries, independent of the talent raster.
# A pure tournament-1 allocation must beat every split between tournament 2
# and the best possible remaining allocation to tournament 1 / outside time.
optout_maxxer_bounds <- function(alpha,beta,phi,scaler,eta,B=50,nx=4001L) {
  xx<-seq(0,B,length.out=nx);xp<-xx[-1]
  self<-vapply(1:2,function(t) {
    if(scaler[t]<=0) return(Inf)
    min(eta/(scaler[t]*phi[t]*beta[t]),
      expm1(log1p(eta*xp/scaler[t])/phi[t])/expm1(beta[t]*log1p(xp)))
  },1)
  branch<-function(t) {
    u<-3-t
    if(scaler[t]<=0) return(rep(NA_real_,length(alpha)))
    vapply(alpha,function(a) {
      v<-scaler[t]*optout_w(xx,a,beta[t],phi[t])-eta*xx
      if(v[nx]<max(v)-1e-10) return(NA_real_)
      if(scaler[u]<=0) return(Inf)
      gap<-v[nx]-cummax(v)[nx:1][-1]
      max(0,min(expm1(log1p(pmax(0,gap+eta*xp)/scaler[u])/phi[u])/
        expm1(beta[u]*log1p(xp))))
    },1)
  }
  list(alpha=alpha,self=self,one=branch(1),two=branch(2))
}

optout_draw_maxxers <- function(bounds,alim,col1='#E06',col2='#08F',col0='black',outline=FALSE) {
  draw<-function(y,swap,col) {
    x<-bounds$alpha;y<-pmin(y,alim[2]);good<-which(!is.na(y)&y>=alim[1])
    if(!length(good)) return(invisible(NULL))
    groups<-split(good,cumsum(c(1,diff(good)>1)))
    for(ii in groups) if(length(ii)>1) {
      px<-c(x[ii],rev(x[ii]));py<-c(rep(alim[1],length(ii)),rev(y[ii]))
      if(swap) {tmp<-px;px<-py;py<-tmp}
      polygon(px,py,col=if(outline) NA else col,
        border=if(outline) col else NA,lwd=1.1,lty=2)
    }
  }
  draw(bounds$one,FALSE,col1);draw(bounds$two,TRUE,col2)
  a<-pmin(bounds$self,alim[2])
  if(all(a>alim[1])) rect(alim[1],alim[1],a[1],a[2],
    col=if(outline) NA else col0,border=if(outline) col0 else NA,lty=2,lwd=1.1)
}

# Equal positive tournament investments with a positive outside allocation.
# A raster sign-contour is insufficient: it can trace a jump between unequal
# optima. Construct stationary equal allocations and verify global optimality.
optout_equal_partial <- function(alpha2,beta,phi,scaler,eta,B,alim) {
  out<-matrix(NA_real_,length(alpha2),2)
  if(any(scaler<=0) || eta<=0) return(out)
  slope<-function(x,a,t) scaler[t]*phi[t]*a*beta[t]*(1+x)^(beta[t]-1)*
    (1+a*expm1(beta[t]*log1p(x)))^(phi[t]-1)
  xx<-sort(unique(c(seq(0,B/2,length.out=129),B/2*.6^(1:30))))
  for(i in seq_along(alpha2)) {
    a2<-alpha2[i];ff<-slope(xx,a2,2)-eta
    ii<-which(ff[-length(ff)]*ff[-1]<0)
    if(!length(ii)) next
    roots<-vapply(ii,function(k) uniroot(function(x) slope(x,a2,2)-eta,
      xx[k+0:1],tol=1e-10)$root,1)
    for(j in seq_along(roots)) {
      x<-roots[j]
      if(x<=1e-8 || x>=B/2-1e-8) next
      f<-function(a) slope(x,a,1)-eta
      if(prod(vapply(alim,f,1))>0) next
      a1<-uniroot(f,alim,tol=1e-10)$root
      z<-optout_allocate(matrix(c(a1,a2),1),beta,phi,scaler,eta,B)
      if(z[1,'x0']>1e-7 && abs(z[1,'x1']-x)<1e-5 && abs(z[1,'x2']-x)<1e-5)
        out[i,min(j,2)]<-a1
    }
  }
  out
}

optout_solve_responsive <- function(beta,phi,eta=.5,B=50,population=optout_population(),tol=2e-5) {
  ref<-optout_fixed_scalers(beta,phi,B)$scaler
  if(eta==1) return(list(status='inactive',scaler=c(0,0),
    state=list(E0=B,E1=0,E2=0,M1=0,M2=0,Utotal=B),err=0))
  residual<-function(logc,return_state=FALSE) {
    if(any(!is.finite(logc)) || any(abs(logc)>45)) return(rep(100,2))
    st<-optout_state(exp(logc),beta,phi,eta,B,population)
    if(any(!is.finite(st$scaler)) || any(st$scaler<=0)) return(rep(100,2))
    if(return_state) return(st)
    log(st$scaler)-logc
  }
  symmetric<-all(abs(diff(beta))<1e-14) && all(abs(diff(phi))<1e-14)
  if(symmetric) {
    f<-function(lc) mean(residual(rep(lc,2)))
    grid<-log(ref[1])+seq(-5,3,length.out=41)
    yy<-vapply(grid,f,1);ii<-which(yy[-length(yy)]*yy[-1]<0 &
      abs(yy[-length(yy)])<50 & abs(yy[-1])<50)
    if(length(ii)>1) stop('Multiple symmetric active roots detected; inspect separately.')
    if(length(ii)) lc<-rep(uniroot(f,grid[ii[1]+0:1],tol=1e-9)$root,2) else lc<-log(ref)
  } else lc<-log(ref)
  # Damped best responses give a useful starting point. Newton refinement
  # then solves the two absolute scaler equations, with a residual line search.
  for(i in 1:45) {
    rr<-residual(lc);if(max(abs(rr))<tol || max(abs(rr))>50) break
    lc<-lc+.25*rr
  }
  for(i in 1:35) {
    rr<-residual(lc);if(max(abs(rr))<tol) break
    delta<-1e-3
    J<-sapply(1:2,function(j) {
      d<-numeric(2);d[j]<-delta
      (residual(lc+d)-residual(lc-d))/(2*delta)
    })
    step<-tryCatch(solve(J,-rr),error=function(e) .15*rr)
    step<-step/max(1,max(abs(step)))
    trial<-lapply(2^-(0:9),function(d) lc+d*step)
    loss<-vapply(trial,function(z) sum(residual(z)^2),1)
    if(min(loss)>=sum(rr^2)) break
    lc<-trial[[which.min(loss)]]
  }
  err<-max(abs(residual(lc)))
  if(err>tol) {
    fit<-optim(lc,function(z) sum(residual(z)^2),method='Nelder-Mead',
      control=list(maxit=350,reltol=1e-12))
    lc<-fit$par;err<-max(abs(residual(lc)))
  }
  if(!is.finite(err) || err>tol) stop('Active responsive solution failed: residual=',signif(err,4),
    '. Increase population or allocation search resolution; no equilibrium figure was saved.')
  list(status='active',scaler=exp(lc),state=residual(lc,TRUE),err=err)
}

# Reuse identical numerical populations between strategy and utility figures.
optout_responsive <- function(beta,phi,eta=.5,B=50,population=optout_population(),tol=2e-5) {
  if(is.null(.GlobalEnv$.optout_cpp$solutions)) .GlobalEnv$.optout_cpp$solutions<-new.env(parent=emptyenv())
  cache<-.GlobalEnv$.optout_cpp$solutions
  key<-paste(c(beta,phi,eta,B,tol,getOption('metatournament.optout.nscan',129L)),collapse='_')
  if(exists(key,cache,inherits=FALSE)) {
    previous<-get(key,cache)
    if(identical(previous$population,population)) return(previous$solution)
  }
  result<-optout_solve_responsive(beta,phi,eta,B,population,tol)
  assign(key,list(population=population,solution=result),cache);result
}
