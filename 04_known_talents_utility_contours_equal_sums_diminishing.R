source("00_helper_functions.R")

#Optimal total utility in alpha1-alpha2 space when talents are known
#Six panels at the local-curvature and endpoint-tie thresholds, beta = 0.4
#
#The plotted quantity is normalized expected total utility:
# U*(alpha1,alpha2)=B*F*(alpha1,alpha2)/E[F*],
# where F*=max_x [w(x,alpha1)+w(B-x,alpha2)].
#With this normalization, mean total utility over the population talent
#distribution is exactly B in every panel.

B<-50
beta<-0.4
alim<-c(0,2.5)

nalpha<-301
nx<-801
nroot<-501
contour_levels<-seq(10,150,10)  #same levels in every panel; edit freely
utility_zlim<-NULL  #set, e.g., c(30,170) for a manual common colour scale
talent_sums<-c(1,2,3)

#Same lognormal talent distribution as in the previous scripts
median<-1
sd<-0.5
m<-log(median)
s<-sqrt(log((1+sqrt(1+4*(sd/median)^2))/2))
nquad<-31

g<-function(x,beta){
  (1+x)^beta-1
}

w<-function(x,alpha,beta,phi){
  (1+alpha*g(x,beta))^phi-1
}

wprime<-function(x,alpha,beta,phi){
  phi*alpha*beta*(1+x)^(beta-1)*(1+alpha*g(x,beta))^(phi-1)
}


gauss_hermite<-function(n=31){
  J<-matrix(0,n,n)
  off<-sqrt((1:(n-1))/2)
  J[cbind(1:(n-1),2:n)]<-off
  J[cbind(2:n,1:(n-1))]<-off
  e<-eigen(J,symmetric=TRUE)
  o<-order(e$values)
  list(z=e$values[o],w=sqrt(pi)*e$vectors[1,o]^2)
}

gh<-gauss_hermite(nquad)
aq<-exp(m+sqrt(2)*s*gh$z)
wq<-gh$w/sqrt(pi)

#Population mean of the globally optimal total spoils weight.
#By symmetry Mphi=E[F*]/2, so normalized total utility is B*F*/E[F*].
mean_optimal_weight<-function(B,beta,phi,nx=801){
  xx<-seq(0,B,length.out=nx)
  ww<-outer(aq,xx,function(a,x){
    (1+a*((1+x)^beta-1))^phi-1
  })
  ww2<-ww[,nx:1,drop=FALSE]

  mu<-0
  for(i in seq_along(aq)){
    uu<-sweep(ww2,2,ww[i,],"+")
    best<-apply(uu,1,max)
    mu<-mu+wq[i]*sum(wq*best)
  }
  mu
}

continuous_boundary<-function(alpha2,B,beta,phi){
  alpha2*(1+B)^(beta-1)*(1+alpha2*g(B,beta))^(phi-1)
}

tie_alpha1<-function(x,alpha2,B,beta,phi){
  if(x<=0) return(continuous_boundary(alpha2,B,beta,phi))
  if(x>=B) return(alpha2)

  rhs<-(1+alpha2*g(B,beta))^phi-
    (1+alpha2*g(B-x,beta))^phi+1

  (rhs^(1/phi)-1)/g(x,beta)
}

stationary_tie<-function(x,alpha2,B,beta,phi){
  alpha1<-tie_alpha1(x,alpha2,B,beta,phi)
  wprime(x,alpha1,beta,phi)-wprime(B-x,alpha2,beta,phi)
}

critical_alpha1<-function(alpha2,B,beta,phi,nroot=501){
  if(alpha2==0) return(c(alpha1=0,x=0))

  if(phi<=1){
    return(c(alpha1=continuous_boundary(alpha2,B,beta,phi),x=0))
  }

  ca<-c(continuous_boundary(alpha2,B,beta,phi),alpha2)
  cx<-c(0,B)

  xx<-seq(1e-7,B-1e-7,length.out=nroot)
  yy<-vapply(xx,stationary_tie,1,
             alpha2=alpha2,B=B,beta=beta,phi=phi)

  k<-which(is.finite(yy[-length(yy)]) &
             is.finite(yy[-1]) &
             yy[-length(yy)]*yy[-1]<0)

  if(length(k)){
    roots<-vapply(k,function(i){
      uniroot(stationary_tie,c(xx[i],xx[i+1]),
              alpha2=alpha2,B=B,beta=beta,phi=phi,
              tol=1e-10)$root
    },1)

    aa<-vapply(roots,tie_alpha1,1,
               alpha2=alpha2,B=B,beta=beta,phi=phi)

    ok<-is.finite(aa) & aa>=0 & aa<=alpha2*(1+1e-8)
    ca<-c(ca,aa[ok])
    cx<-c(cx,roots[ok])
  }

  i<-which.min(ca)
  c(alpha1=ca[i],x=cx[i])
}

boundary_curve<-function(alpha,B,beta,phi,nroot=501){
  if(phi<=1){
    return(continuous_boundary(alpha,B,beta,phi))
  }

  vapply(alpha,function(alpha2){
    critical_alpha1(alpha2,B,beta,phi,nroot)["alpha1"]
  },1)
}

#Optimal total spoils weight. Maxxer regions are evaluated analytically;
#only the diversifier band is searched numerically.
utility_map<-function(alpha,B,beta,phi,nx=801,nroot=501){
  n<-length(alpha)
  bound<-boundary_curve(alpha,B,beta,phi,nroot)
  u<-matrix(NA_real_,n,n)

  xx<-seq(0,B,length.out=nx)
  ww<-outer(alpha,xx,function(a,x){
    (1+a*((1+x)^beta-1))^phi-1
  })
  ww2<-ww[,nx:1,drop=FALSE]

  for(j in seq_len(n)){
    if(j>1){
      ii<-1:(j-1)

      #2maxxers: utility depends only on alpha2
      maxxer<-ii[alpha[ii]<=bound[j]]
      if(length(maxxer)){
        u[maxxer,j]<-w(B,alpha[j],beta,phi)
      }

      #Diversifiers: numerically maximize over allocation
      inside<-ii[alpha[ii]>bound[j]]
      if(length(inside)){
        uu<-sweep(ww[inside,,drop=FALSE],2,ww2[j,],"+")
        u[inside,j]<-apply(uu,1,max)
      }

      #Symmetry under exchange of tournament labels
      u[j,ii]<-u[ii,j]
    }

    #Equal talents
    u[j,j]<-max(ww[j,]+ww2[j,])
  }

  list(u=u,bound=bound)
}

#Tournament colours used in the corresponding investment figures
col2_max<-"#08F"
col1_max<-"#E06"
coldiag_ref<-"#D0D0D0"

#Neutral shading: darker areas have higher optimal utility
#ucols<-gray.colors(256,start=0.97,end=0.35)
vir <- hcl.colors(256, "Viridis")
ucols <- colorRampPalette(vir)(256)

i <- round((seq(0, 1, length.out = 256)^0.65) * 255) + 1
ucols <- vir[i]


#Compute thresholds at full precision; round only the panel titles.
training_curvature_factor<-(1-beta)/beta*(1-(1+B/2)^(-beta))
phi_curvature<-(1+training_curvature_factor-2/(talent_sums*g(B/2,beta)))/(2+training_curvature_factor)
phi_tie<-vapply(talent_sums,function(S){
  uniroot(function(phi){
    2*expm1(phi*log1p(S/2*g(B/2,beta)))-
      expm1(phi*log1p(S*g(B,beta)))
  },c(1e-6,1),tol=1e-12)$root
},1)
phis<-c(phi_curvature,phi_tie)
print(data.frame(S=talent_sums,phi_curvature,phi_tie),digits=10)

draw_utility_map<-function(){
  alpha<-seq(alim[1],alim[2],length.out=nalpha)
  panels<-vector("list",length(phis))
  zlim_global<-c(Inf,-Inf)
  for(k in seq_along(phis)){
    phi<-phis[k]
    um<-utility_map(alpha,B,beta,phi,nx,nroot)
    meanF<-mean_optimal_weight(B,beta,phi,nx)
    u<-B*um$u/meanF
    S<-talent_sums[(k-1)%%3+1]
    #For these phi<1 panels the equal-talent optimum is the 50:50 split.
    balanced_utility<-B*2*w(B/2,S/2,beta,phi)/meanF
    balanced_index<-which.min(abs(alpha-S/2))
    stopifnot(abs(u[balanced_index,balanced_index]-balanced_utility)<1e-10)
    panels[[k]]<-list(phi=phi,u=u,bound=um$bound,
                       balanced_utility=balanced_utility)
    zlim_global<-range(zlim_global,u,finite=TRUE)
  }
  if(!is.null(utility_zlim)) zlim_global<-utility_zlim

  #Rotated row descriptions use the wider left outer margin.
  layout(matrix(1:6,nrow=2,byrow=TRUE))
  par(oma=c(2,6,2.8,0.8))
  for(k in seq_along(phis)){
    i<-(k-1)%/%3+1
    j<-(k-1)%%3+1
    phi<-phis[k]
    u<-panels[[k]]$u
    bound<-panels[[k]]$bound
    par(mar=c(3.1,2.75,2.8,0.6),
        mgp=c(1.75,0.48,0),tcl=-0.23,pty="s",xpd=FALSE)
    image(alpha,alpha,u,col=ucols,zlim=zlim_global,
          xlim=alim,ylim=alim,xaxs="i",yaxs="i",
          axes=FALSE,xlab="",ylab="",useRaster=TRUE)
    for(ts in talent_sums) abline(ts,-1,col="white",lty=3,lwd=1)
    abline(0,1,col=coldiag_ref,lty=2,lwd=0.9)
    lev<-contour_levels[contour_levels>min(u,na.rm=TRUE) &
                        contour_levels<max(u,na.rm=TRUE)]
    contour2(alpha,alpha,u,levels=lev,add=TRUE,
             drawlabels=TRUE,labcex=1.1,lwd=0.8,
             method="flattest",label.shift=0.15,col=1)
    #Highlight the utility contour through alpha1=alpha2=S/2.
    contour(alpha,alpha,u,levels=panels[[k]]$balanced_utility,
            add=TRUE,drawlabels=FALSE,col=2,lwd=1.5)
    lines(bound,alpha,lty=2,lwd=1.1,col=adjustcolor(col2_max,0.9))
    lines(alpha,bound,lty=2,lwd=1.1,col=adjustcolor(col1_max,0.9))
    box()
    if(i==2){
      axis(1,at=seq(alim[1],alim[2],0.5),cex.axis=1.2)
      mtext(expression(alpha[1]~~"(talent in tournament 1)"),
            side=1,line=2.4,cex=0.95)
    }
    if(j==1){
      axis(2,at=seq(alim[1],alim[2],0.5),las=1,cex.axis=1.2)
      mtext(expression(alpha[2]~~"(talent in tournament 2)"),
            side=2,line=2.4,cex=0.95)
    }
    mtext(bquote(varphi==.(formatC(phi,format="f",digits=3))),
          side=3,line=0.95,cex=1.4)
  }
  for(j in seq_along(talent_sums)){
    mtext(bquote(alpha[1]+alpha[2]==.(talent_sums[j])),
          side=3,outer=TRUE,at=(j-0.5)/3,line=0.7,cex=1.2)
  }
  mtext("Local curvature thresholds",side=2,outer=TRUE,
        at=0.75,line=2.3,las=0,cex=1.05)
  mtext("Balanced and concentrated\ntalents have equal utility",
        side=2,outer=TRUE,at=0.25,line=2.3,las=0,cex=1.05)
  mtext(bquote(beta==.(beta)),side=2,outer=TRUE,at=0.5,line=2,cex=1.3)
}

png("Figure_known_talent_utility_contours_equal_sums_diminishing.png",
    width=8.4,height=7,units="in",res=300)
draw_utility_map()
dev.off()

# if(interactive()) draw_utility_map()
