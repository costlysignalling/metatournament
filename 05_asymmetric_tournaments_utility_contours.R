# Source this file from the scripts folder, then call plot_asymmetric_utility().
# Both figures share asymmetric_examples from 05_asymmetric_examples.R.
# Draws in R's plotting window; export only with png=TRUE and/or pdf=TRUE.
# Population moments use the same switch-split solver as explore_asymmetric().
# Increase nquad to 96 or 128 to check population convergence (minimum 64).
source("05_asymmetric_examples.R")
source("05_random_allocation_reference.R")
source("05_linear_population_means.R")
source("05_explorer_population.R")
source("00_helper_functions.R")
source("05_asymmetric_utility_refinement.R")

plot_asymmetric_utility <- function(scenarios=asymmetric_examples,B=50,
    alim=c(0,2.5),nalpha=301L,nx=801L,nroot=501L,
    nquad=64L,nxpop=801L,population=NULL,
    png=FALSE,pdf=FALSE,file="Figure_asymmetric_tournaments_utility_contours.png",
    contour_levels=seq(10,150,10),utility_zlim=NULL) {
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

critical_alpha1<-function(alpha2,B,beta1,beta2,phi1,phi2,s,nroot=501){
  if(alpha2==0) return(0)

  if(phi1<=1 && phi2<=1){
    return(tie_alpha1(0,alpha2,B,beta1,beta2,phi1,phi2,s))
  }

  xx<-seq(0,B,length.out=nroot)
  yy<-vapply(xx,tie_alpha1,1,
             alpha2=alpha2,B=B,beta1=beta1,beta2=beta2,
             phi1=phi1,phi2=phi2,s=s)

  cand<-c(yy[1],yy[length(yy)])

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

boundary_curves<-function(alpha,B,beta1,beta2,phi1,phi2,s,nroot=501){
  b2<-vapply(alpha,critical_alpha1,1,
             B=B,beta1=beta1,beta2=beta2,
             phi1=phi1,phi2=phi2,s=s,nroot=nroot)

  b1<-vapply(alpha,critical_alpha1,1,
             B=B,beta1=beta2,beta2=beta1,
             phi1=phi2,phi2=phi1,s=1/s,nroot=nroot)

  list(b2=b2,b1=b1)
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
#Utility calculations
#------------------------------------------------------------

#Population mean of the globally optimal objective F*=max_x[s w1+w2].
#Optimal normalized utility map. Maxxer regions are evaluated analytically;
#only the diversifier band is searched numerically.
utility_map_asym<-function(alpha,B,beta1,beta2,phi1,phi2,s,nx=801,nroot=501){
  n<-length(alpha)
  bc<-boundary_curves(alpha,B,beta1,beta2,phi1,phi2,s,nroot)

  u<-matrix(NA_real_,n,n)
  max2<-matrix(FALSE,n,n)
  max1<-matrix(FALSE,n,n)

  tol<-1e-10

  for(j in seq_len(n)){
    is2<-alpha<=bc$b2[j]*(1+tol)
    is1<-alpha[j]<=bc$b1*(1+tol)
    both<-is1 & is2

    if(any(is2 & !both)){
      ii<-which(is2 & !both)
      max2[ii,j]<-TRUE
      u[ii,j]<-w(B,alpha[j],beta2,phi2)
    }

    if(any(is1 & !both)){
      ii<-which(is1 & !both)
      max1[ii,j]<-TRUE
      u[ii,j]<-s*w(B,alpha[ii],beta1,phi1)
    }

    if(any(both)){
      ii<-which(both)
      u2<-rep(w(B,alpha[j],beta2,phi2),length(ii))
      u1<-s*w(B,alpha[ii],beta1,phi1)
      choose1<-u1>u2

      if(any(choose1)){
        jj<-ii[choose1]
        max1[jj,j]<-TRUE
        u[jj,j]<-u1[choose1]
      }
      if(any(!choose1)){
        jj<-ii[!choose1]
        max2[jj,j]<-TRUE
        u[jj,j]<-u2[!choose1]
      }
    }

  }

  inside <- which(!max1 & !max2,arr.ind=TRUE)
  if(nrow(inside)) u[inside] <- asymmetric_utility_values(
    alpha[inside[,1]],alpha[inside[,2]],B,beta1,beta2,phi1,phi2,s,nx)

  list(u=u,b2=bc$b2,b1=bc$b1,max2=max2,max1=max1)
}

#------------------------------------------------------------
#Prepare the four scenarios
#------------------------------------------------------------

prep <- prepare_asymmetric_examples(scenarios,B,nquad,nxpop,population)


#------------------------------------------------------------
#Colours and utility shading
#------------------------------------------------------------

col2_div<-"#315A9B"
colmid<-"#FFF"
col1_div<-"#B54343"
col2_max<-"#08F"
col1_max<-"#E06"
coldiag<-"#FF0"

#Viridis shading, shifted slightly towards lighter tones
vir<-hcl.colors(256,"Viridis")
i<-round((seq(0,1,length.out=256)^0.65)*255)+1
ucols<-vir[i]

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

#------------------------------------------------------------
#Figure
#------------------------------------------------------------

panel_cache <- NULL
draw_asymmetric_utility<-function(){
  if(is.null(panel_cache)) {
  alpha<-seq(alim[1],alim[2],length.out=nalpha)

  #Precompute all panels so the contour levels and colour scale are shared.
  panels<-vector("list",2*length(prep))
  k<-0
  zlim_global<-c(Inf,-Inf)

  for(i in 1:2){
    for(j in seq_along(prep)){
      obj<-prep[[j]]
      sc<-obj$pars

      if(i==1){
        s_use<-obj$fixed$s
        E1<-obj$response$E1
        E2<-obj$response$E2
        M2<-obj$fixed$M2
        utility_scale<-obj$fixed$E2/M2
      } else {
        s_use<-unname(obj$eq["s"])
        E1<-unname(obj$eq["E1"])
        E2<-unname(obj$eq["E2"])
        M2<-unname(obj$eq["M2"])
        utility_scale<-unname(obj$eq["E2"])/M2
      }

      um<-utility_map_asym(alpha,B,
                           sc$beta1,sc$beta2,sc$phi1,sc$phi2,
                           s_use,nx,nroot)

      moments <- if(i==1) obj$response else as.list(obj$eq)
      meanF <- s_use*moments$M1+moments$M2
      u<-utility_scale*um$u # Actual utility under this population reference.

      hc<-half_curve(alpha,B,
                     sc$beta1,sc$beta2,sc$phi1,sc$phi2,
                     s_use,alim,nx)

      k<-k+1
      panels[[k]]<-list(beta1=sc$beta1,beta2=sc$beta2,
                        phi1=sc$phi1,phi2=sc$phi2,
                        s=s_use,E1=E1,E2=E2,M2=M2,mean_utility=utility_scale*meanF,
                        u=u,b2=um$b2,b1=um$b1,hc=hc)

      zlim_global[1]<-min(zlim_global[1],min(u,na.rm=TRUE))
      zlim_global[2]<-max(zlim_global[2],max(u,na.rm=TRUE))
    }
  }

  if(!is.null(utility_zlim)) zlim_global<-utility_zlim

  panel_cache <<- list(panels=panels,zlim=zlim_global,alpha=alpha)
  } else {
    panels <- panel_cache$panels; zlim_global <- panel_cache$zlim
    alpha <- panel_cache$alpha
  }

  nc <- length(prep)
  layout(cbind(matrix(seq_len(2*nc),nrow=2,byrow=TRUE),2*nc+1L),
         widths=c(rep(1,nc),0.46),heights=c(1,1))

  par(oma=c(2,3.5,3,0.0))

  k<-0
  for(i in 1:2){
    for(j in seq_along(prep)){
      k<-k+1
      pan<-panels[[k]]

      par(mar=c(2.75,2.75,0.5,0.2),
          mgp=c(1.75,0.48,0),tcl=-0.23,pty="s",xpd=FALSE)

      image(alpha,alpha,pan$u,col=ucols,zlim=zlim_global,
            xlim=alim,ylim=alim,xaxs="i",yaxs="i",
            axes=FALSE,xlab="",ylab="",useRaster=TRUE)
      
      abline(0,1,col="#D0D0D0",lty=2,lwd=0.9)

      lev<-contour_levels[contour_levels>min(pan$u,na.rm=TRUE) &
                          contour_levels<max(pan$u,na.rm=TRUE)]

      contour2(alpha,alpha,pan$u,levels=lev,add=TRUE,
               drawlabels=TRUE,labcex=1.05,lwd=0.8,
               method="flattest",label.shift=0.15,col=1)

      #Exact global specialization boundaries.
      lines(pan$b2,alpha,col=adjustcolor(col2_max,0.9),lwd=1.05,lty=2)
      lines(alpha,pan$b1,col=adjustcolor(col1_max,0.9),lwd=1.05,lty=2)

      # Stationary 50:50 curve, including sections hidden by specialization
      # in the strategy map. Dashed sections do not imply optimal equal splitting.
      dashed <- pan$hc$x
      dashed[pan$hc$global] <- NA_real_
      lines(dashed,pan$hc$y,col=coldiag,lwd=0.9,lty=2)
      # Solid yellow only where 50:50 is globally optimal.
      ii<-which(is.finite(pan$hc$x) & pan$hc$global)
      if(length(ii)>1){
        groups<-split(ii,cumsum(c(1,diff(ii)!=1)))
        for(r in groups){
          if(length(r)>1){
            lines(pan$hc$x[r],pan$hc$y[r],col=coldiag,lwd=2.5)
          }
        }
      }

      box()
      # Explicitly reset margin clipping after contour labels and prior devices.
      par(xpd=NA)
      draw_budget_bar(pan$E1,pan$E2,alim)

      if(i==2){
        axis(1,at=pretty(alim),cex.axis=1.2)
        mtext(xpd=NA,expression(alpha[1]~~"(talent in tournament 1)"),
              side=1,line=2.4,cex=0.95)
      }

      if(j==1){
        axis(2,at=pretty(alim),las=1,cex.axis=1.2)
        mtext(xpd=NA,expression(alpha[2]~~"(talent in tournament 2)"),
              side=2,line=2.4,cex=0.95)
      }

      if(i==1){
        # One bold label per example column; aligned with the parameter headers.
        mtext(LETTERS[(j-1L)%%26L+1L],side=3,line=1.55,
              at=alim[1],adj=0,cex=1.25,font=2,col="black",xpd=NA)
        tx1<-alim[1]+0.37*diff(alim)
        tx2<-alim[1]+0.77*diff(alim)

        mtext(xpd=NA,bquote(beta[1]==.(formatC(pan$beta1,format="f",digits=2))),
              side=3,line=1.55,at=tx1,cex=1.12,col=col1_max)
        mtext(xpd=NA,bquote(beta[2]==.(formatC(pan$beta2,format="f",digits=2))),
              side=3,line=1.55,at=tx2,cex=1.12,col=col2_max)

        mtext(xpd=NA,bquote(varphi[1]==.(formatC(pan$phi1,format="f",digits=2))),
              side=3,line=0.15,at=tx1,cex=1.12,col=col1_max)
        mtext(xpd=NA,bquote(varphi[2]==.(formatC(pan$phi2,format="f",digits=2))),
              side=3,line=0.15,at=tx2,cex=1.12,col=col2_max)
      }
    }
  }

  mtext(xpd=NA,expression("Fixed "*E[1]==E[2]),
        side=2,outer=TRUE,at=0.76,line=1.25,cex=1.25)

  mtext(xpd=NA,expression("Responsive "*E[t]),
        side=2,outer=TRUE,at=0.26,line=1.25,cex=1.25)

  #Keep the right panel empty so the layout matches the investment figure.
  par(mar=c(2.5,0.55,0.5,1.60),pty="m",xpd=NA)
  plot.new()
}


  files <- asymmetric_figure_output(draw_asymmetric_utility,length(prep),file,png,pdf)
  invisible(list(population=prep,draw=draw_asymmetric_utility,files=files))
}

# Examples (both functions accept the identical list):
# result <- plot_asymmetric_utility()
# result <- plot_asymmetric_utility(scenarios=asymmetric_examples,png=TRUE,pdf=TRUE)

# Compatible entry point for previous interactive use.
draw_asymmetric_utility <- plot_asymmetric_utility
