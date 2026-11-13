# Geometry of the stationary 50:50 curve and alpha1=alpha2.
# Positive log(R) means that the curve lies above the diagonal in the plot.
# The supported phi pairs are equal phis, or phis on opposite sides of 1
# (including 1). For these pairs log(R) is monotone, so there is at most one
# positive crossing. Other pairs can have two crossings and require extension.
equal_talent_geometry <- function(beta1,beta2,phi1,phi2,ratio,B=50,tol=1e-11) {
  if(phi1!=phi2 && (phi1-1)*(phi2-1)>0)
    stop("This phi pair can have two crossings; use equal phis or phis straddling 1.")
  n <- max(length(beta1),length(beta2),length(ratio))
  b1 <- rep_len(beta1,n); b2 <- rep_len(beta2,n); r <- rep_len(ratio,n)
  if(any(!is.finite(c(b1,b2,r))) || min(b1,b2,r)<=0 || max(b1,b2)>=1)
    stop("Invalid beta or population ratio.")
  L <- 1+B/2; G1 <- expm1(b1*log(L)); G2 <- expm1(b2*log(L))
  low <- log(r)+log(phi1*b1/(phi2*b2))+(b1-b2)*log(L)
  high <- if(phi1==phi2) low+(phi1-1)*log(G1/G2) else rep(sign(phi1-phi2)*Inf,n)
  low[abs(low)<tol] <- 0
  high[is.finite(high)&abs(high)<tol] <- 0
  coincident <- (phi1==phi2)&(low==0)&(high==0)
  crossing <- !coincident & low*sign(high)<0
  sign_high <- sign(high); sign_high[sign_high==0] <- sign(low[sign_high==0])
  side <- sign(low); side[side==0] <- sign_high[side==0]
  # 0: coincident; 1: always above; 2: always below;
  # 3: crosses, above at large talents; 4: crosses, below at large talents.
  type <- ifelse(coincident,0L,ifelse(side>0,1L,2L))
  type[crossing] <- ifelse(sign_high[crossing]>0,3L,4L)
  alpha_cross <- rep(NA_real_,n)
  k <- which(crossing)
  if(length(k)) {
    lo <- rep(-64,length(k)); hi <- rep(64,length(k))
    f <- function(z) low[k]+(phi1-1)*log1p(exp(z)*G1[k])-
      (phi2-1)*log1p(exp(z)*G2[k])
    flo <- f(lo); fhi <- f(hi)
    if(any(flo*fhi>0)) stop("Crossing outside the supported floating-point bracket.")
    for(i in 1:75) {
      mid <- (lo+hi)/2; fm <- f(mid)
      same <- sign(fm)==sign(flo)
      lo[same] <- mid[same]; flo[same] <- fm[same]
      hi[!same] <- mid[!same]
    }
    alpha_cross[k] <- exp((lo+hi)/2)
  }
  data.frame(beta1=b1,beta2=b2,ratio=r,type=as.integer(type),
             logR_low=low,logR_high=high,alpha_cross=alpha_cross,
             sum_cross=2*alpha_cross)
}

# A stationary equal split need not be a global optimum when phi>1.
# Test it against the whole investment interval, including both endpoints.
check_equal_split_global <- function(field,phi1,phi2,B=50,nx=801L) {
  ans <- rep(NA,length(field$type))
  ids <- which(is.finite(field$alpha_cross))
  if(!length(ids)) return(ans)
  if(phi1<=1 && phi2<=1) { ans[ids] <- TRUE; return(ans) }
  xx <- seq(0,B,length.out=nx)
  for(ii in split(ids,ceiling(seq_along(ids)/200L))) {
    a <- field$alpha_cross[ii]
    G1 <- expm1(outer(field$beta1[ii],log1p(xx)))
    G2 <- expm1(outer(field$beta2[ii],log1p(B-xx)))
    U <- field$ratio[ii]*expm1(phi1*log1p(a*G1))+expm1(phi2*log1p(a*G2))
    best <- U[cbind(seq_along(ii),max.col(U,ties.method='first'))]
    midpoint <- field$ratio[ii]*expm1(phi1*log1p(a*expm1(field$beta1[ii]*log1p(B/2))))+
      expm1(phi2*log1p(a*expm1(field$beta2[ii]*log1p(B/2))))
    ans[ii] <- best<=midpoint+1e-8*pmax(1,abs(midpoint))
  }
  ans
}
