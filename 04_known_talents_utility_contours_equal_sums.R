source("00_helper_functions.R")

#Optimal total utility in alpha1-alpha2 space when talents are known
#Eight panels: beta = 0.4, 0.8 and varphi = 0.79, 1, 1.21, 1.42
#
#The plotted quantity is normalized expected total utility:
# U*(alpha1,alpha2)=B*F*(alpha1,alpha2)/E[F*],
# where F*=max_x [w(x,alpha1)+w(B-x,alpha2)].
#With this normalization, mean total utility over the population talent
#distribution is exactly B in every panel.

B<-50
betas<-c(0.4,0.8)
phis<-c(0.66,1,1.25,1.5)
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

draw_utility_map<-function(){
  alpha<-seq(alim[1],alim[2],length.out=nalpha)

  #Precompute all panels so the colour scale can be shared across them
  panels<-vector("list",length(betas)*length(phis))
  k<-0
  zlim_global<-c(Inf,-Inf)

  for(i in seq_along(betas)){
    beta<-betas[i]

    for(j in seq_along(phis)){
      phi<-phis[j]
      um<-utility_map(alpha,B,beta,phi,nx,nroot)
      meanF<-mean_optimal_weight(B,beta,phi,nx)
      u<-B*um$u/meanF

      k<-k+1
      panels[[k]]<-list(beta=beta,phi=phi,u=u,bound=um$bound)

      zlim_global[1]<-min(zlim_global[1],min(u,na.rm=TRUE))
      zlim_global[2]<-max(zlim_global[2],max(u,na.rm=TRUE))
    }
  }

  if(!is.null(utility_zlim)) zlim_global<-utility_zlim

  layout(matrix(c(1,2,3,4,9,
                  5,6,7,8,9),nrow=2,byrow=TRUE),
         widths=c(1,1,1,1,0.46),heights=c(1,1))

  par(oma=c(2,3.5,3,0.0))

  k<-0
  for(i in seq_along(betas)){
    beta<-betas[i]

    for(j in seq_along(phis)){
      phi<-phis[j]
      k<-k+1
      u<-panels[[k]]$u
      bound<-panels[[k]]$bound

      par(mar=c(2.75,2.75,0.5,0.2),
          mgp=c(1.75,0.48,0),tcl=-0.23,pty="s",xpd=FALSE)

      image(alpha,alpha,u,col=ucols,zlim=zlim_global,
            xlim=alim,ylim=alim,xaxs="i",yaxs="i",
            axes=FALSE,xlab="",ylab="",useRaster=TRUE)

      #Lines connecting individuals with equal total talent.
      for(ts in talent_sums) abline(ts,-1,col="white",lty=3,lwd=1)

      #Equal-talent reference line, deliberately below contours and boundaries.
      abline(0,1,col=coldiag_ref,lty=2,lwd=0.9)

      #Same manually specified utility levels in every panel.
      #contour() simply ignores levels outside the range of the current panel.
      lev<-contour_levels[contour_levels>min(u,na.rm=TRUE) &
                          contour_levels<max(u,na.rm=TRUE)]

      contour2(alpha,alpha,u,levels=lev,add=TRUE,
               drawlabels=TRUE,labcex=1.1,lwd=0.8,
               method="flattest",label.shift=0.15,col=1)

      #Analytical specialization boundaries. Beyond these boundaries the
      #contours become exactly horizontal or vertical.
      lines(bound,alpha,lty=2,lwd=1.1,col=adjustcolor(col2_max,0.9))
      lines(alpha,bound,lty=2,lwd=1.1,col=adjustcolor(col1_max,0.9))

      box()

      if(i==length(betas)){
        axis(1,at=seq(alim[1],alim[2],0.5),cex.axis=1.2)
        mtext(expression(alpha[1]~~"(talent in tournament 1)"),
              side=1,line=2.4,cex=0.95)
      }

      if(j==1){
        axis(2,at=seq(alim[1],alim[2],0.5),las=1,cex.axis=1.2)
        mtext(expression(alpha[2]~~"(talent in tournament 2)"),
              side=2,line=2.4,cex=0.95)
      }

      if(i==1){
        f<-function(x) formatC(x,format="f",digits=2)
        mtext(bquote(varphi==.(f(phi))),side=3,line=0.95,cex=1.4)
      }
    }
  }

  mtext(expression(beta==0.4),side=2,outer=TRUE,at=0.76,line=1.25,cex=1.4)
  mtext(expression(beta==0.8),side=2,outer=TRUE,at=0.26,line=1.25,cex=1.4)

  #Keep the right panel empty so the layout matches the investment figures.
  par(mar=c(2.5,0.55,0.5,1.60),pty="m",xpd=NA)
  plot.new()
}


# pdf("Figure_known_talent_utility_contours.pdf",width=12.6,height=5.8)
# draw_utility_map()
# dev.off()

png("Figure_known_talent_utility_contours.png",width=3580,height=1740,res=300)
draw_utility_map()
dev.off()

# if(interactive()) draw_utility_map()
