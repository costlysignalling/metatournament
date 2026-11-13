# Opting out with random-allocation fixed E or self-consistent responsive E.
# The absolute E_t/M_t matters because eta supplies an outside return.
# Set options(metatournament.eta=.5) or options(metatournament.eta=1).
# PNG and PDF include eta automatically; defaults: both TRUE, PNG 600 dpi.
# optout model uses continuous allocation candidates and deterministic QMC.
source("07_optout_model.R")
source("00_helper_functions.R")
B<-50
betas<-c(0.4,0.8)
phis<-c(0.66,1,1.25,1.5)
eta<-getOption("metatournament.eta",.5)
alim<-c(0,2.5)

nalpha<-getOption("metatournament.optout.nalpha",251L)
nx<-801
nxb<-4001
contour_levels<-seq(10,150,10)  #same levels in every panel; edit freely
utility_zlim<-NULL               #set manually if desired

#Baseline talent distribution
median<-1
sd<-0.5
m<-log(median)
s<-sqrt(log((1+sqrt(1+4*(sd/median)^2))/2))

nxpop<-801

#------------------------------------------------------------
#Core functions
#------------------------------------------------------------

g<-function(x,beta){
  (1+x)^beta-1
}

w<-function(x,alpha,beta,phi){
  (1+alpha*g(x,beta))^phi-1
}

#------------------------------------------------------------
#Fixed equal-E random-allocation reference from the canonical eta=0 population
#------------------------------------------------------------


# The fixed reference is random allocation, independent of talents.
benchmark_scaler<-function(beta,phi,B,nx=801){
  optout_fixed_scalers(beta,phi,B)$scaler[1]
}
scaler_grid<-outer(betas,phis,Vectorize(function(beta,phi){
  benchmark_scaler(beta,phi,B,nxpop)
}))

#------------------------------------------------------------
#Global optimizer in talent space
#------------------------------------------------------------

#Running maximum along each row. This allows the triangular constraint
#x1+x2<=B to be solved without evaluating the full allocation triangle
#separately for every talent pair.
prefix_max<-function(v){
  n<-ncol(v)
  val<-v
  idx<-matrix(1L,nrow(v),n)

  for(k in 2:n){
    take<-v[,k]>val[,k-1]
    val[,k]<-ifelse(take,v[,k],val[,k-1])
    idx[,k]<-ifelse(take,k,idx[,k-1])
  }

  list(val=val,idx=idx)
}

#Optimal focal utility U*(alpha1,alpha2) against the fixed equal-E random-allocation reference.
#is not changed when self-directed time is introduced. The constant eta*B is added after maximizing
#V1(x1)+V2(x2), where Vt(x)=scaler*w(x,alpha)-eta*x.
utility_map_optout<-function(alpha,B,beta,phi,scaler,eta,nx=801){
  optout_talent_map(alpha,rep(beta,2),rep(phi,2),rep(scaler,2),eta,B)$u
}
#------------------------------------------------------------
#Pure maxxer boundaries, matching the strategy-map figure
#------------------------------------------------------------

maxxer_boundaries<-function(alpha_boundary,B,beta,phi,scaler,eta,nxb=4001){
  xx<-seq(0,B,length.out=nxb)
  gx<-g(xx,beta)

  xpos<-xx[-1]
  self_alpha<-min(eta/(scaler*phi*beta),
    ((1+eta*xpos/scaler)^(1/phi)-1)/g(xpos,beta))

  other_max<-rep(NA_real_,length(alpha_boundary))
  cap<-nxb:1

  for(i in seq_along(alpha_boundary)){
    a1<-alpha_boundary[i]
    v1<-scaler*w(xx,a1,beta,phi)-eta*xx
    v1B<-v1[nxb]

    if(v1B < max(v1)-1e-11) next

    p1<-cummax(v1)
    x2<-xx[-1]
    best1<-p1[cap[-1]]
    gap<-v1B-best1

    rhs<-1+(gap+eta*x2)/scaler
    aa<-rep(Inf,length(x2))
    ok<-rhs>0 & gx[-1]>0
    aa[ok]<-(rhs[ok]^(1/phi)-1)/gx[-1][ok]
    other_max[i]<-max(0,min(aa,na.rm=TRUE))
  }

  list(alpha=alpha_boundary,self=self_alpha,max_other=other_max)
}

#Population-average budget allocation, copied in spirit from the strategy-map figure.
population_summary_optout<-function(beta,phi,scaler,eta,B,nx=801,npop=201){
  st<-optout_state(rep(scaler,2),rep(beta,2),rep(phi,2),eta,B)
  c(E0=st$E0,E1=st$E1,E2=st$E2,M1=st$M1,M2=st$M2,
    Utour=scaler*(st$M1+st$M2),Uself=eta*st$E0,Utotal=st$Utotal,
    tournament_sink=B-scaler*(st$M1+st$M2))
}
#------------------------------------------------------------
#Colours and side bars
#------------------------------------------------------------

col2_max<-"#08F"
col1_max<-"#E06"
self_max<-"#000000"
coldiag_ref<-"#D0D0D0"

vir<-hcl.colors(256,"Viridis")
i<-round((seq(0,1,length.out=256)^0.65)*255)+1
ucols<-vir[i]

draw_population_bar<-function(E0,E1,E2,alim,
                               x0=NULL,x1=NULL,
                               col0=self_max,col1=col1_max,col2=col2_max){
  dx<-diff(alim)
  dy<-diff(alim)

  if(is.null(x0)) x0<-alim[2]+0.015*dx
  if(is.null(x1)) x1<-alim[2]+0.075*dx

  y0<-alim[1]
  y1<-y0+(E0/B)*dy
  y2<-y1+(E1/B)*dy
  y3<-alim[2]

  rect(x0,y0,x1,y1,col=col0,border=NA,xpd=NA)
  rect(x0,y1,x1,y2,col=col1,border=NA,xpd=NA)
  rect(x0,y2,x1,y3,col=col2,border=NA,xpd=NA)
  rect(x0,y0,x1,y3,border="black",col=NA,lwd=0.7,xpd=NA)

  xm<-(x0+x1)/2
  if(y1-y0>0.055*dy)
    text(xm,(y0+y1)/2,expression(E[0]),cex=0.67,xpd=NA,col="white",font=2)
  if(y2-y1>0.055*dy)
    text(xm,(y1+y2)/2,expression(E[1]),cex=0.67,xpd=NA,col="white",font=2)
  if(y3-y2>0.055*dy)
    text(xm,(y2+y3)/2,expression(E[2]),cex=0.67,xpd=NA,col="white",font=2)
}

#Complete outlines of the three pure maxxer regions.
#
#The curved line max_other(alpha1) is only one part of a tournament-maxxer
#boundary. For sufficiently winner-take-more cases, complete maxxing can appear
#discontinuously: the finite maxxer branch then begins at alpha1>0 with
#max_other(alpha1)>0. This creates an additional vertical boundary segment for
#the 1maxxer region and, by symmetry, a horizontal segment for the 2maxxer
#region. These edges are visually obvious in the filled strategy map but must
#be drawn explicitly in a utility-contour plot.
draw_maxxer_boundaries<-function(bounds,alim,lty=2,lwd=1.1){
  a<-bounds$self
  x<-bounds$alpha
  y<-pmin(bounds$max_other,alim[2])

  valid<-is.finite(y) & y>=alim[1] &
    x>=alim[1] & x<=alim[2]

  idx<-which(valid)

  if(length(idx)){
    groups<-split(idx,cumsum(c(1,diff(idx)!=1)))
    eps<-1e-10*max(1,diff(alim))

    for(ii in groups){
      if(length(ii)<2) next

      xx<-x[ii]
      yy<-y[ii]

      #Curved global-tie boundaries.
      lines(xx,yy,lty=lty,lwd=lwd,
            col=adjustcolor(col1_max,0.9))
      lines(yy,xx,lty=lty,lwd=lwd,
            col=adjustcolor(col2_max,0.9))

      #If a finite maxxer branch starts inside the plotted talent range,
      #complete the outline by drawing its discontinuous onset edge.
      if(xx[1]>alim[1]+eps && yy[1]>alim[1]+eps){
        segments(xx[1],alim[1],xx[1],yy[1],
                 lty=lty,lwd=lwd,col=adjustcolor(col1_max,0.9))
        segments(alim[1],xx[1],yy[1],xx[1],
                 lty=lty,lwd=lwd,col=adjustcolor(col2_max,0.9))
      }

      #Likewise close a branch if it terminates before the upper plotting limit.
      if(xx[length(xx)]<alim[2]-eps && yy[length(yy)]>alim[1]+eps){
        segments(xx[length(xx)],alim[1],
                 xx[length(xx)],yy[length(yy)],
                 lty=lty,lwd=lwd,col=adjustcolor(col1_max,0.9))
        segments(alim[1],xx[length(xx)],
                 yy[length(yy)],xx[length(xx)],
                 lty=lty,lwd=lwd,col=adjustcolor(col2_max,0.9))
      }
    }
  }

  #The remaining two interior sides of the self-maxxer square.
  #Its lower and left sides coincide with the plot frame.
  if(is.finite(a) && a>alim[1] && a<alim[2]){
    segments(alim[1],a,a,a,lty=lty,lwd=lwd,col=self_max)
    segments(a,alim[1],a,a,lty=lty,lwd=lwd,col=self_max)
  }
}

#------------------------------------------------------------
#Figure
#------------------------------------------------------------

draw_utility_map_optout<-function(){
  alpha<-seq(alim[1],alim[2],length.out=nalpha)
  alpha_boundary<-seq(alim[1],alim[2],length.out=1201)

  panels<-vector("list",length(betas)*length(phis))
  k<-0
  zlim_global<-c(Inf,-Inf)

  for(i in seq_along(betas)){
    beta<-betas[i]

    for(j in seq_along(phis)){
      phi<-phis[j]
      scaler<-scaler_grid[i,j]

      u<-utility_map_optout(alpha,B,beta,phi,scaler,eta,nx)
      bounds<-maxxer_boundaries(alpha_boundary,B,beta,phi,scaler,eta,nxb)
      pa<-population_summary_optout(beta,phi,scaler,eta,B,nxpop)

      k<-k+1
      panels[[k]]<-list(beta=beta,phi=phi,u=u,bounds=bounds,pa=pa)
      zlim_global[1]<-min(zlim_global[1],min(u,na.rm=TRUE))
      zlim_global[2]<-max(zlim_global[2],max(u,na.rm=TRUE))
    }
  }

  if(!is.null(utility_zlim)) zlim_global<-utility_zlim

  cat("\nPopulation summary under the fixed equal-E random-allocation reference:\n")
  for(k in seq_along(panels)){
    pan<-panels[[k]]
    pa<-pan$pa
    cat("beta =",pan$beta," phi =",pan$phi,
        " E0/B =",formatC(pa["E0"]/B,digits=3,format="f"),
        " Utour/B =",formatC(pa["Utour"]/B,digits=3,format="f"),
        " Uself/B =",formatC(pa["Uself"]/B,digits=3,format="f"),
        " Utotal/B =",formatC(pa["Utotal"]/B,digits=3,format="f"),"\n")
  }

  layout(matrix(c(1,2,3,4,9,
                  5,6,7,8,9),nrow=2,byrow=TRUE),
         widths=c(1,1,1,1,0.52),heights=c(1,1))

  eta_at<-(4+0.52/2)/(4+0.52)
  par(oma=c(2,3.5,3,0.0))

  k<-0
  for(i in seq_along(betas)){
    for(j in seq_along(phis)){
      k<-k+1
      u<-panels[[k]]$u
      bounds<-panels[[k]]$bounds
      pa<-panels[[k]]$pa

      par(mar=c(2.75,2.75,0.5,0.2),
          mgp=c(1.75,0.48,0),tcl=-0.23,pty="s",xpd=FALSE)

      image(alpha,alpha,u,col=ucols,zlim=zlim_global,
            xlim=alim,ylim=alim,xaxs="i",yaxs="i",
            axes=FALSE,xlab="",ylab="",useRaster=TRUE)

      abline(0,1,col=coldiag_ref,lty=2,lwd=0.9)

      lev<-contour_levels[contour_levels>min(u,na.rm=TRUE) &
                            contour_levels<max(u,na.rm=TRUE)]
      if(length(lev)) contour2(alpha,alpha,u,levels=lev,add=TRUE,
               drawlabels=TRUE,labcex=1.1,lwd=0.8,
               method="flattest",label.shift=0.15,col=1)

      draw_maxxer_boundaries(bounds,alim)
      box()
      draw_population_bar(pa["E0"],pa["E1"],pa["E2"],alim)

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
        mtext(bquote(varphi==.(f(phis[j]))),side=3,line=0.95,cex=1.4)
      }
    }
  }

  mtext(expression(beta==0.4),side=2,outer=TRUE,at=0.76,line=1.25,cex=1.4)
  mtext(expression(beta==0.8),side=2,outer=TRUE,at=0.26,line=1.25,cex=1.4)

  par(mar=c(2.5,0.55,0.5,1.60),pty="m",xpd=NA)
  plot.new()
  mtext(bquote(eta==.(formatC(eta,format="f",digits=2))),
        side=3,outer=TRUE,at=eta_at,line=0.50,cex=1.5,font=2)
  mtext(expression("Fixed "*E[1]==E[2]),
        side=3,outer=TRUE,at=-0.03,line=0.5,cex=1.3,adj=0,col="#802")
}

optout_save(draw_utility_map_optout,"Figure_known_talent_optout_utility_contours_fixed_E",eta)
