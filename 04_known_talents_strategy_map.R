#Optimal allocation in alpha1-alpha2 space when talents are known
#Analytical maxxer boundaries + numerical allocation only inside the diversifier region

B<-50
betas<-c(0.4,0.8)
phis<-c(0.66,1,1.25,1.5)
alim<-c(0,2.5)

nalpha<-301
nx<-801
nroot<-501
noverlap<-2  #number of alpha-grid cells rasterized underneath maxxer polygons

g<-function(x,beta){
  (1+x)^beta-1
}

w<-function(x,alpha,beta,phi){
  (1+alpha*g(x,beta))^phi-1
}

wprime<-function(x,alpha,beta,phi){
  phi*alpha*beta*(1+x)^(beta-1)*(1+alpha*g(x,beta))^(phi-1)
}

#If the transition out of 2maxxing is continuous, its boundary follows from
#w'(0,alpha1)=w'(B,alpha2).
continuous_boundary<-function(alpha2,B,beta,phi){
  alpha2*(1+B)^(beta-1)*(1+alpha2*g(B,beta))^(phi-1)
}

#For any candidate interior allocation x, this is the talent alpha1 for which
#that allocation ties exactly with complete 2maxxing:
#w(x,alpha1)+w(B-x,alpha2)=w(B,alpha2).
tie_alpha1<-function(x,alpha2,B,beta,phi){
  if(x<=0) return(continuous_boundary(alpha2,B,beta,phi))
  if(x>=B) return(alpha2)

  rhs<-(1+alpha2*g(B,beta))^phi-
    (1+alpha2*g(B-x,beta))^phi+1

  (rhs^(1/phi)-1)/g(x,beta)
}

#At an interior minimum of the tie threshold, the tying allocation is also
#stationary in utility.
stationary_tie<-function(x,alpha2,B,beta,phi){
  alpha1<-tie_alpha1(x,alpha2,B,beta,phi)
  wprime(x,alpha1,beta,phi)-wprime(B-x,alpha2,beta,phi)
}

#Global 2maxxer boundary for a given alpha2.
#
#The critical alpha1 is the minimum talent needed for ANY alternative
#allocation to tie the 2maxxer. Candidate minima are:
#  x -> 0 : continuous transition,
#  0 < x < B : discontinuous transition,
#  x = B : direct tie between the two maxxer strategies.
critical_alpha1<-function(alpha2,B,beta,phi,nroot=501){
  if(alpha2==0) return(c(alpha1=0,x=0))

  #For phi<=1 the objective is concave, so only the continuous boundary matters.
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

#Compute the optimal relative allocation only inside the analytically delimited
#diversifier region. The two maxxer regions are never searched numerically.
strategy_map<-function(alpha,B,beta,phi,nx=801,nroot=501,noverlap=2){
  n<-length(alpha)
  bound<-boundary_curve(alpha,B,beta,phi,nroot)
  z<-matrix(NA_real_,n,n)
  overlap<-noverlap*(alpha[2]-alpha[1])

  xx<-seq(0,B,length.out=nx)
  ww<-outer(alpha,xx,function(a,x){
    (1+a*((1+x)^beta-1))^phi-1
  })
  ww2<-ww[,nx:1,drop=FALSE]

  #Lower triangle: alpha1 < alpha2.
  #Numerically evaluate a narrow strip just inside the analytical maxxer
  #region as well. The polygon is drawn on top, so this overlap removes
  #anti-aliased/raster seams along the exact boundary.
  for(j in 2:n){
    inside<-which(alpha>bound[j]-overlap & alpha<alpha[j]-1e-12)

    if(length(inside)){
      uu<-sweep(ww[inside,,drop=FALSE],2,ww2[j,],"+")
      k<-max.col(uu,ties.method="first")
      z[inside,j]<-xx[k]/B

      #Tournament-label symmetry gives the upper triangle for free.
      z[j,inside]<-1-z[inside,j]
    }
  }

  #Equal talents: select one global optimum by convention.
  #The separate split_global vector below records whether 50:50 itself is global.
  mid<-(nx+1)/2
  split_global<-logical(n)

  for(j in seq_len(n)){
    uu<-ww[j,]+ww2[j,]
    k<-which.max(uu)
    z[j,j]<-xx[k]/B
    split_global[j]<-abs(uu[mid]-max(uu))<
      1e-9*max(1,abs(max(uu)))
  }

  list(z=z,bound=bound,split_global=split_global)
}

#Diversifier gradient is deliberately more muted than the solid maxxer colours
# col2_div<-"#0A0"
# colmid<-"#FFF"
# col1_div<-"#C0C"
# col2_max<-"#0F0"
# col1_max<-"#F0F"
# coldiag<-"#FF0"
col2_div<-"#315A9B"
colmid<-"#FFF"
col1_div<-"#B54343"
col2_max<-"#08F"
col1_max<-"#E06"
coldiag<-"#FF0"
cols<-colorRampPalette(c(col2_div,colmid,col1_div))(801)

draw_strategy_map<-function(){
  alpha<-seq(alim[1],alim[2],length.out=nalpha)

  layout(matrix(c(1,2,3,4,9,
                  5,6,7,8,9),nrow=2,byrow=TRUE),
         widths=c(1,1,1,1,0.46),heights=c(1,1))

  par(oma=c(2,3.5,3,0.0))

  for(i in seq_along(betas)){
    beta<-betas[i]

    for(j in seq_along(phis)){
      phi<-phis[j]
      sm<-strategy_map(alpha,B,beta,phi,nx,nroot,noverlap)
      z<-sm$z
      bound<-sm$bound

      par(mar=c(2.75,2.75,0.5,0.2),
          mgp=c(1.75,0.48,0),tcl=-0.23,pty="s",xpd=FALSE)

      #Only diversifiers are rasterized.
      image(alpha,alpha,z,col=cols,zlim=c(0,1),xlim=alim,ylim=alim,
            xaxs="i",yaxs="i",axes=FALSE,xlab="",ylab="",useRaster=TRUE)
      
      #abline(0,1,col="#D0D0D0",lty=2,lwd=0.9)

      #Analytically delimited maxxer regions.
      polygon(c(bound,0,0),c(alpha,max(alpha),0),
              col=col2_max,border=NA)
      polygon(c(alpha,max(alpha),0),c(bound,0,0),
              col=col1_max,border=NA)

      #Show the exact 50:50 strategy only where it is globally optimal.
      splitidx<-which(sm$split_global)
      if(length(splitidx)>1){
        groups<-split(splitidx,cumsum(c(1,diff(splitidx)!=1)))
        for(r in groups){
          if(length(r)>1){
            lines(alpha[r],alpha[r],col=coldiag,lwd=3)
          }
        }
      }

      box()

      if(i==length(betas)){
        axis(1,at=seq(0,2.5,0.5),cex.axis=1.2)
        mtext(expression(alpha[1]~~"(talent in tournament 1)"),
              side=1,line=2.4,cex=0.95)
      }

      if(j==1){
        axis(2,at=seq(0,2.5,0.5),las=1,cex.axis=1.2)
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

  #Common vertical colour legend
  par(mar=c(2.5,0.55,0.5,1.60),pty="m",xpd=NA)
  plot.new()
  plot.window(xlim=c(-0.08,1.00),ylim=c(0,1),xaxs="i",yaxs="i")

  zz<-seq(0,1,length.out=801)
  for(k in 1:800){
    rect(0.34,zz[k],0.64,zz[k+1],col=cols[k],border=NA)
  }
  rect(0.34,0,0.64,1,border="black",lwd=0.8)

  wd<-0.006
  rect(0.34,0-wd,0.64,0+wd,col=col2_max,lwd=0.8,xpd=NA,border=1)
  rect(0.34,0.5-wd,0.64,0.5+wd,col=coldiag,lwd=0.8,xpd=NA,border=1)
  rect(0.34,1-wd,0.64,1+wd,col=col1_max,lwd=0.8,xpd=NA,border=1)

  ticks<-seq(0,1,0.25)
  segments(0.64,ticks,0.75,ticks,lwd=0.8)
  text(0.88,ticks,sprintf("%.2f",ticks),adj=c(0,0.5),cex=1.2,xpd=NA)

  text(0.10,0.5,"Relative time invested in tournament 1",
       srt=90,cex=1.7,xpd=NA)
  text(0.49,1.035,"1maxxers",cex=1.5,xpd=NA,col=col1_max)
  text(0.49,-0.035,"2maxxers",cex=1.5,xpd=NA,col=col2_max)
}

# pdf("Figure_known_talent_strategy_map.pdf",width=12.6,height=5.8)
# draw_strategy_map()
# dev.off()

png("Figure_known_talent_strategy_map.png",width=3580,height=1740,res=300)
draw_strategy_map()
dev.off()

# if(interactive()) draw_strategy_map()


# For arbitrary talents, the global optimum can be calculated directly in one dimension. A dense grid followed by local refinement is useful because for $\varphi>1$ the utility landscape can contain more than one local maximum.
# 
# ```{r optimal-allocation}
# g <- function(x,beta){
#   (1+x)^beta-1
# }
# 
# w <- function(x,alpha,beta,phi){
#   (1+alpha*g(x,beta))^phi-1
# }
# 
# best_allocation <- function(alpha1,alpha2,B,beta,phi){
#   f <- function(x){
#     w(x,alpha1,beta,phi)+w(B-x,alpha2,beta,phi)
#   }
#   
#   xx <- seq(0,B,length.out=5001)
#   yy <- f(xx)
#   i <- which.max(yy)
#   
#   if(i==1) return(0)
#   if(i==length(xx)) return(B)
#   
#   optimize(function(x){-f(x)},
#            c(xx[i-1],xx[i+1]),
#            tol=1e-10)$minimum
# }
# 
# best_allocation(1,2,B=50,beta=0.4,phi=1.5)
# best_allocation(5,10,B=50,beta=0.4,phi=1.5)
# ```