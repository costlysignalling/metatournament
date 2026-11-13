# Source this file from the scripts folder, then call plot_asymmetric_hdr().
# Both figures share asymmetric_examples from 05_asymmetric_examples.R.
# Draws in R's plotting window; export only with png=TRUE and/or pdf=TRUE.
# Population moments use the same switch-split solver as explore_asymmetric().
# Increase nquad to 96 or 128 to check population convergence (minimum 64).
source("05_asymmetric_examples.R")
source("05_random_allocation_reference.R")
source("05_linear_population_means.R")
source("05_explorer_population.R")

plot_asymmetric_hdr <- function(scenarios=asymmetric_examples,B=50,
    alim=c(0,2.5),nalpha=301L,nx=801L,nroot=501L,
    nquad=64L,nxpop=801L,population=NULL,
    png=FALSE,pdf=FALSE,file="Figure_asymmetric_tournaments.png") {
  validate_asymmetric_examples(scenarios,B,alim,nalpha,nx,nroot,nquad,nxpop,png,pdf)
  nx <- as.integer(nx); if(nx%%2==0) nx <- nx+1L
  nxpop <- as.integer(nxpop); if(nxpop%%2==0) nxpop <- nxpop+1L
  median <- 1; sd <- 0.5; m <- 0
  slog <- sqrt(log((1+sqrt(2))/2))
#------------------------------------------------------------
#Core functions
#------------------------------------------------------------

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

    # Refine local peaks before deciding whether the stationary split is global.
    peaks <- which(u[2:(nx-1)]>=u[1:(nx-2)] & u[2:(nx-1)]>=u[3:nx])+1L
    objective <- function(x) s*w(x,alpha1[j],beta1,phi1)+w(B-x,a2,beta2,phi2)
    best <- max(c(u, vapply(peaks,function(k)
      optimize(objective,c(xx[k-1],xx[k+1]),maximum=TRUE,tol=1e-10)$objective,1)))
    global[j] <- abs(objective(mid)-best)<1e-8*max(1,abs(best))
  }

  list(x=alpha1,y=alpha2,global=global)
}

#------------------------------------------------------------
#Prepare the four scenarios
#------------------------------------------------------------

prep <- prepare_asymmetric_examples(scenarios,B,nquad,nxpop,population)


cat("Population mean budget shares:\n")
for(j in seq_along(prep)){
  x<-prep[[j]]

  cat("\nScenario",j,"\n")
  cat("  first response:   E1/B =",
      formatC(x$response$E1/B,digits=3,format="f"),
      " E2/B =",formatC(x$response$E2/B,digits=3,format="f"),"\n")
  cat("  responsive E_t: E1/B =",
      formatC(x$eq["E1"]/B,digits=3,format="f"),
      " E2/B =",formatC(x$eq["E2"]/B,digits=3,format="f"),
      " (E1/M1)/(E2/M2) =",format(x$eq["s"],digits=6),"\n")
}

#------------------------------------------------------------
#Colours
#------------------------------------------------------------

col2_div<-"#315A9B"
colmid<-"#FFF"
col1_div<-"#B54343"
col2_max<-"#08F"
col1_max<-"#E06"
coldiag<-"#FF0"

cols<-colorRampPalette(c(col2_div,colmid,col1_div))(801)

#Stacked side bar showing the expected budget shares E1/B and E2/B.
#Tournament 1 is the bottom rectangle, tournament 2 is the top rectangle.
draw_budget_bar<-function(E1,E2,alim,
                          x0=NULL,x1=NULL,
                          col1=col1_max,col2=col2_max){
  dx<-diff(alim)
  dy<-diff(alim)

  if(is.null(x0)) x0<-alim[2]+0.015*dx
  if(is.null(x1)) x1<-alim[2]+0.075*dx

  y0<-alim[1]
  y1<-alim[1]+(E1/B)*dy
  y2<-alim[2]

  rect(x0,y0,x1,y1,col=col1,border=NA,lwd=0.7,xpd=NA)
  rect(x0,y1,x1,y2,col=col2,border=NA,lwd=0.7,xpd=NA)

  xm<-(x0+x1)/2
  text(xm,alim[1],expression(E[1]),pos=3,cex=0.75,xpd=NA,col="white",font=2)
  text(xm,alim[2],expression(E[2]),pos=1,cex=0.75,xpd=NA,col="white",font=2)
}


#Highest-density region contours of the bivariate lognormal talent distribution.
#These are the exact iso-density curves enclosing the requested probability mass.
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
#------------------------------------------------------------
#Figure
#------------------------------------------------------------

map_cache <- new.env(parent=emptyenv())
draw_asymmetric_maps<-function(){
  default_alim<-alim
  hdr_curves<-talent_hdr_curves(median,sd,probs=c(0.50,0.90,0.99))

  nc <- length(prep)
  layout(cbind(matrix(seq_len(2*nc),nrow=2,byrow=TRUE),2*nc+1L),
         widths=c(rep(1,nc),0.46),heights=c(1,1))

  par(oma=c(2,3.5,3,0.0))

  for(i in 1:2){
    for(j in seq_along(prep)){
      obj<-prep[[j]]
      sc<-obj$pars
      alim<-if(is.null(sc$alim))default_alim else sc$alim
      alpha<-seq(alim[1],alim[2],length.out=nalpha)
      tick_scale<-if(max(abs(alim))>=1000)1000 else 1

      if(i==1){
        s_use<-obj$fixed$s
        E1<-obj$response$E1
        E2<-obj$response$E2
      } else {
        s_use<-unname(obj$eq["s"])
        E1<-unname(obj$eq["E1"])
        E2<-unname(obj$eq["E2"])
      }

      key <- paste(i,j,sep=":")
      if(!exists(key,map_cache,inherits=FALSE))
        assign(key,strategy_map(alpha,B,sc$beta1,sc$beta2,sc$phi1,sc$phi2,
                                s_use,nx,nroot),map_cache)
      sm <- get(key,map_cache,inherits=FALSE)

      par(mar=c(2.75,2.75,0.5,0.2),
          mgp=c(1.75,0.48,0),tcl=-0.23,pty="s",xpd=FALSE)

      #Muted continuous gradient for diversifiers.
      image(alpha,alpha,sm$z,col=cols,zlim=c(0,1),
            xlim=alim,ylim=alim,xaxs="i",yaxs="i",
            axes=FALSE,xlab="",ylab="",useRaster=TRUE)

      #Bright solid colours only for analytically identified global maxxers.
      z2<-ifelse(sm$max2,1,NA_real_)
      z1<-ifelse(sm$max1,1,NA_real_)

      image(alpha,alpha,z2,col=col2_max,breaks=c(0.5,1.5),
            add=TRUE,axes=FALSE,useRaster=TRUE)
      image(alpha,alpha,z1,col=col1_max,breaks=c(0.5,1.5),
            add=TRUE,axes=FALSE,useRaster=TRUE)
      
      abline(0,1,col="#000",lty=2,lwd=0.9)

      #Reference HDR contours of the talent distribution.
      for(hc in hdr_curves){
        lines(hc[1,],hc[2,],col="#EEE",lwd=0.8)
      }

      #Exact global maxxer boundaries.
      lines(sm$b2,alpha,col=adjustcolor(col2_max,0.9),lwd=1.05)
      lines(alpha,sm$b1,col=adjustcolor(col1_max,0.9),lwd=1.05)

      #Exact 50:50 curve, but only where 50:50 is globally optimal.
      hc<-half_curve(alpha,B,
                     sc$beta1,sc$beta2,sc$phi1,sc$phi2,
                     s_use,alim,nx)

      ii<-which(is.finite(hc$x) & hc$global)

      if(length(ii)>1){
        groups<-split(ii,cumsum(c(1,diff(ii)!=1)))

        for(r in groups){
          if(length(r)>1){
            lines(hc$x[r],hc$y[r],col=coldiag,lwd=3)
          }
        }
      }

      box()
      draw_budget_bar(E1,E2,alim)

      if(i==2){
        axis(1,at=pretty(alim),labels=if(tick_scale==1) TRUE else pretty(alim)/tick_scale,cex.axis=1.2)
        mtext(expression(alpha[1]~~"(talent in tournament 1)"),
              side=1,line=2.4,cex=0.95)
      }

      if(j==1 || !identical(alim,default_alim)){
        axis(2,at=pretty(alim),labels=if(tick_scale==1) TRUE else pretty(alim)/tick_scale,las=1,cex.axis=1.2)
        if(j==1)mtext(expression(alpha[2]~~"(talent in tournament 2)"),
              side=2,line=2.4,cex=0.95)
      }

      if(tick_scale>1)mtext(expression("" %*% 10^3),side=3,line=.15,
        at=alim[1],adj=0,cex=.75,xpd=NA)

      if(i==1){
        # One bold label per example column; aligned with the parameter headers.
        mtext(LETTERS[(j-1L)%%26L+1L],side=3,line=1.55,
              at=alim[1],adj=0,cex=1.25,font=2,col="black",xpd=NA)
        tx1<-alim[1]+0.37*diff(alim)
        tx2<-alim[1]+0.77*diff(alim)

        mtext(bquote(beta[1]==.(formatC(sc$beta1,format="f",digits=2))),
              side=3,line=1.55,at=tx1,cex=1.12,col=col1_max)
        mtext(bquote(beta[2]==.(formatC(sc$beta2,format="f",digits=2))),
              side=3,line=1.55,at=tx2,cex=1.12,col=col2_max)

        mtext(bquote(varphi[1]==.(formatC(sc$phi1,format="f",digits=2))),
              side=3,line=0.15,at=tx1,cex=1.12,col=col1_max)
        mtext(bquote(varphi[2]==.(formatC(sc$phi2,format="f",digits=2))),
              side=3,line=0.15,at=tx2,cex=1.12,col=col2_max)
      }
    }
  }

  mtext(expression("Fixed "*E[1]==E[2]),
        side=2,outer=TRUE,at=0.76,line=1.25,cex=1.25)

  mtext(expression("Responsive "*E[t]),
        side=2,outer=TRUE,at=0.26,line=1.25,cex=1.25)

  #Common strategy legend
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
  text(0.88,ticks,sprintf("%.2f",ticks),
       adj=c(0,0.5),cex=1.2,xpd=NA)

  text(0.10,0.5,"Relative time invested in tournament 1",
       srt=90,cex=1.7,xpd=NA)

  text(0.49,1.035,"1maxxers",
       cex=1.5,xpd=NA,col=col1_max)

  text(0.49,-0.035,"2maxxers",
       cex=1.5,xpd=NA,col=col2_max)

}


  files <- asymmetric_figure_output(draw_asymmetric_maps,length(prep),file,png,pdf)
  invisible(list(population=prep,draw=draw_asymmetric_maps,files=files))
}

# Examples (both functions accept the identical list):
# result <- plot_asymmetric_hdr()
# result <- plot_asymmetric_hdr(scenarios=asymmetric_examples,png=TRUE,pdf=TRUE)

# Compatible entry point for previous interactive use.
draw_asymmetric_maps <- plot_asymmetric_hdr

