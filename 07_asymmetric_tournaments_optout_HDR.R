# Opting out with random-allocation fixed E or self-consistent responsive E.
# The absolute E_t/M_t matters because eta supplies an outside return.
# Set options(metatournament.eta=.5) or options(metatournament.eta=1).
# PNG and PDF include eta automatically; defaults: both TRUE, PNG 600 dpi.
# optout model uses continuous allocation candidates and deterministic QMC.
source("07_optout_model.R")
source("00_helper_functions.R")
B<-50
eta<-getOption("metatournament.eta",.5)

# EXACTLY the four scenarios used in 05_asymmetric_tournaments_figure_HDR.R
source("05_asymmetric_examples.R")
scenarios<-asymmetric_examples

alim<-c(0,2.5)
nalpha<-getOption("metatournament.optout.nalpha",251L)
nx<-801

# Talent distribution
median<-1
sd<-0.5
m<-log(median)
slog<-sqrt(log((1+sqrt(1+4*(sd/median)^2))/2))

# Population integration

nxpop<-801

# Responsive two-scaler solver.
# Several damping values are tried because the global best-response map is
# piecewise smooth and can jump when strategies switch.
eq_tol<-1e-4
eq_maxit<-400
eq_damps<-c(0.35,0.20,0.10,0.05)
inactive_tol<-1e-12

#------------------------------------------------------------
# Core functions
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

#------------------------------------------------------------
# Population expectations
#------------------------------------------------------------


population_cache<-function(beta1,beta2,phi1,phi2,B,nx=801){
  list(beta=c(beta1,beta2),phi=c(phi1,phi2))
}
# Random-allocation reference, with a separate E_t/M_t for each tournament.
benchmark_scaler<-function(cache,B){
  optout_fixed_scalers(cache$beta,cache$phi,B)$scaler
}
# Population best responses to arbitrary scaler1,scaler2.
population_state<-function(scaler1,scaler2,cache,eta,B){
  optout_state(c(scaler1,scaler2),cache$beta,cache$phi,eta,B)
}
# One damped attempt at scaler_t = E_t/M_t.
# Responsive solver is shared in 07_optout_model.R.
# Multistart version for manual exploration.
responsive_state<-function(cache,eta,B,k0){
  eq<-optout_responsive(cache$beta,cache$phi,eta,B)
  list(status=eq$status,scaler=eq$scaler,state=eq$state,err=eq$err)
}
#------------------------------------------------------------
# Talent-space strategy map
#------------------------------------------------------------

strategy_map_optout<-function(alpha,B,beta1,beta2,phi1,phi2,scaler1,scaler2,eta,nx=801){
  optout_talent_map(alpha,c(beta1,beta2),c(phi1,phi2),c(scaler1,scaler2),eta,B,with_bounds=TRUE)
}
#------------------------------------------------------------
# Exact 50:50 stationary curve, retained from Figure 05 logic
#------------------------------------------------------------

# At x0=0 and x1=x2=B/2, the full-budget reallocation derivative is
# scaler1*w1' - scaler2*w2' = 0. We plot this curve only where that candidate
# is also globally optimal against BOTH tournaments and the outside option.
best_utility_pair<-function(alpha1,alpha2,B,beta1,beta2,phi1,phi2,scaler1,scaler2,eta,nx=801){
  optout_allocate(matrix(c(alpha1,alpha2),1),c(beta1,beta2),c(phi1,phi2),
    c(scaler1,scaler2),eta,B)[1,'utility']
}
half_curve<-function(alpha2,B,beta1,beta2,phi1,phi2,
                     scaler1,scaler2,eta,alim,nx=801){

  if(scaler1<=0 || scaler2<=0) return(list(x=rep(NA_real_,length(alpha2)),y=alpha2,global=rep(FALSE,length(alpha2))))
  mid<-B/2

  f<-function(alpha1,a2){
    scaler1*wprime(mid,alpha1,beta1,phi1)-
      scaler2*wprime(mid,a2,beta2,phi2)
  }

  alpha1<-rep(NA_real_,length(alpha2))
  global<-logical(length(alpha2))

  for(j in seq_along(alpha2)){
    a2<-alpha2[j]

    if(a2==0){
      alpha1[j]<-0
    } else {
      lo<-f(alim[1],a2)
      hi<-f(alim[2],a2)

      if(!is.finite(lo) || !is.finite(hi) || lo*hi>0) next

      alpha1[j]<-uniroot(function(a1){
        f(a1,a2)
      },alim,tol=1e-10)$root
    }

    a1<-alpha1[j]

    uhalf<-scaler1*w(mid,a1,beta1,phi1)+
      scaler2*w(mid,a2,beta2,phi2)

    umax<-best_utility_pair(
      a1,a2,B,beta1,beta2,phi1,phi2,
      scaler1,scaler2,eta,nx
    )

    global[j]<-abs(uhalf-umax)<1e-8*max(1,abs(umax))
  }

  list(x=alpha1,y=alpha2,global=global)
}

#------------------------------------------------------------
# HDR contours
#------------------------------------------------------------

talent_hdr_curves<-function(median=1,sd=0.5,rho=0,
                            probs=c(0.50,0.90,0.99)){
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
# Colours -- COPIED from Figure 07
#------------------------------------------------------------

# When x0=0 this reproduces the old blue-white-red diversifier gradient.
# As x0/B increases, the colour is blended toward self_col. Complete
# self-maxxers and the two tournament maxxer types are overlaid in solid
# colours, and the global 50:50 splitter line is highlighted in yellow.
col2_div<-"#315A9B"
colmid<-"#FFF"
col1_div<-"#B54343"
col2_max<-"#08F"
col1_max<-"#E06"
self_col<-"#555555"
self_max<-"#000000"
coldiag<-"#FF0"

base_cols<-colorRampPalette(c(col2_div,colmid,col1_div))(801)
base_rgb<-col2rgb(base_cols)/255
self_rgb<-drop(col2rgb(self_col)/255)

allocation_colours<-function(x0,x1,x2){
  qt<-pmax(1-x0,0)
  r1<-matrix(0.5,nrow(x0),ncol(x0))
  used<-qt>1e-12
  r1[used]<-x1[used]/qt[used]

  idx<-pmin(pmax(round(800*r1)+1,1),801)
  out<-matrix(NA_character_,nrow(x0),ncol(x0))

  for(i in seq_len(nrow(x0))){
    br<-base_rgb[,idx[i,],drop=FALSE]
    q<-qt[i,]
    rr<-sweep(br,2,q,"*")+self_rgb%o%(1-q)
    out[i,]<-rgb(rr[1,],rr[2,],rr[3,])
  }

  out
}

# Draw a raster whose matrix coordinates follow image(x=alpha,y=alpha,z=...).
draw_colour_raster<-function(alpha,colmat){
  ras<-as.raster(t(colmat)[ncol(colmat):1,,drop=FALSE])
  rasterImage(ras,min(alpha),min(alpha),max(alpha),max(alpha),interpolate=FALSE)
}

# Solid overlays for the three complete maxxer classes.
draw_maxxer_masks<-function(alpha,sm){
  optout_draw_maxxers(sm$bounds,range(alpha),col1_max,col2_max,self_max)
}
# Stacked side bar -- COPIED from Figure 07.
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

# Two-dimensional legend -- COPIED from Figure 07.
draw_allocation_legend<-function(){
  ns<-181
  nr<-181
  s0<-seq(0,1,length.out=ns)
  r1<-seq(0,1,length.out=nr)

  cm<-matrix(NA_character_,ns,nr)
  idx<-pmin(pmax(round(800*r1)+1,1),801)

  for(i in seq_along(s0)){
    q<-1-s0[i]
    br<-base_rgb[,idx,drop=FALSE]
    rr<-q*br+self_rgb%o%rep(s0[i],nr)
    cm[i,]<-rgb(rr[1,],rr[2,],rr[3,])
  }

  plot.new()
  plot.window(xlim=c(0,1.1),ylim=c(0,1),xaxs="i",yaxs="i")

  ras<-as.raster(t(cm)[nr:1,,drop=FALSE])
  rasterImage(ras,0,0,1,1,interpolate=FALSE)
  box()

  axis(1,at=c(0,0.5,1),labels=c("","0.5","1.0"),
       cex.axis=0.9,tcl=-0.20)
  axis(2,at=seq(0,1,0.1),
       labels=formatC(seq(0,1,0.1),format="f",digits=1),
       las=1,cex.axis=0.9,tcl=-0.15)

  wd<-0.15
  hd<-0.018
  rect(0,1-hd,wd,1,col=col1_max,border=1,lwd=0.8,xpd=NA)
  rect(0,0.5-hd/2,wd,0.5+hd/2,border=1,col="#FF0",lwd=0.8,xpd=NA)
  rect(0,0,wd,hd,border=1,col=col2_max,lwd=0.8,xpd=NA)

  rect(1,0,1+wd*0.8,1,col=self_max,border=1,lwd=0.8,xpd=NA)

  mtext(expression(x[0]/B),
        side=1,line=2.2,cex=0.9)
  mtext("(self-directed",
        side=1,line=3.25,cex=0.9)
  mtext("share)",
        side=1,line=4.3,cex=0.9)

  mtext("Relative tournament-1 share among tournament time",
        side=2,line=1.65,cex=0.9)

  text(0.15,0.025,"2",adj=c(0,0),cex=1.0,col="white",font=2)
  text(0.15,0.975,"1",adj=c(0,1),cex=1.0,col="white",font=2)
  text(0.96,0.50,"self",adj=c(1,0.5),cex=0.95,col="white",font=2)
  text(1.25,0.5,"selfmaxxers",srt=90,adj=c(0.5,0.5),cex=1,xpd=NA)
  text(-0.09,1.032,"1maxxers",cex=1,xpd=NA,col=col1_max)
  text(-0.09,-0.032,"2maxxers",cex=1,xpd=NA,col=col2_max)
}

#------------------------------------------------------------
# Prepare the SAME four scenarios as Figure 05
#------------------------------------------------------------

prep<-lapply(scenarios,function(sc){
  cache<-population_cache(sc$beta1,sc$beta2,sc$phi1,sc$phi2,B,nxpop)

  k0<-benchmark_scaler(cache,B)
  fixed<-population_state(k0[1],k0[2],cache,eta,B)
  eq<-responsive_state(cache,eta,B,k0)

  list(pars=sc,cache=cache,k0=k0,fixed=fixed,eq=eq)
})

cat("Population mean budget shares:\n")
for(j in seq_along(prep)){
  z<-prep[[j]]

  cat("\nScenario",j,"\n")
  cat("  fixed equal E: E0/B =",
      formatC(z$fixed$E0/B,digits=3,format="f"),
      " E1/B =",formatC(z$fixed$E1/B,digits=3,format="f"),
      " E2/B =",formatC(z$fixed$E2/B,digits=3,format="f"),
      " scaler =",formatC(z$k0,digits=3,format="f"),"\n")

  cat("  responsive:       E0/B =",
      formatC(z$eq$state$E0/B,digits=3,format="f"),
      " E1/B =",formatC(z$eq$state$E1/B,digits=3,format="f"),
      " E2/B =",formatC(z$eq$state$E2/B,digits=3,format="f"),
      " scaler1 =",formatC(z$eq$scaler[1],digits=3,format="f"),
      " scaler2 =",formatC(z$eq$scaler[2],digits=3,format="f"),
      " [",z$eq$status,"]\n")
}

#------------------------------------------------------------
# Figure
#------------------------------------------------------------

draw_asymmetric_optout_maps<-function(){
  alpha<-seq(alim[1],alim[2],length.out=nalpha)
  hdr_curves<-talent_hdr_curves(median,sd,probs=c(0.50,0.90,0.99))

  # EXACT Figure-07 layout.
  legw<-0.52
  widths<-c(1,1,1,1,legw)

  layout(
    matrix(c(1,2,3,4,9,
             5,6,7,8,9),nrow=2,byrow=TRUE),
    widths=widths,
    heights=c(1,1)
  )

  eta_at<-(sum(widths[1:4])+widths[5]/2)/sum(widths)

  # EXACT Figure-07 outer margins.
  par(oma=c(2,3.5,3,0.2))

  for(i in 1:2){
    for(j in 1:4){
      obj<-prep[[j]]
      sc<-obj$pars

      if(i==1){
        scaler1<-obj$k0[1]
        scaler2<-obj$k0[2]
        pa<-obj$fixed
      } else {
        scaler1<-obj$eq$scaler[1]
        scaler2<-obj$eq$scaler[2]
        pa<-obj$eq$state
      }

      sm<-strategy_map_optout(
        alpha,B,
        sc$beta1,sc$beta2,sc$phi1,sc$phi2,
        scaler1,scaler2,eta,nx
      )

      cm<-allocation_colours(sm$x0,sm$x1,sm$x2)

      # EXACT Figure-07 panel graphics.
      par(
        mar=c(2.75,2.75,0.5,0.2),
        mgp=c(1.75,0.48,0),
        tcl=-0.23,
        pty="s",
        xpd=FALSE
      )

      plot(
        NA,
        xlim=alim,ylim=alim,
        xaxs="i",yaxs="i",
        axes=FALSE,xlab="",ylab="",asp=1
      )

      # Figure-07 diversifier/self-directed colour logic.
      draw_colour_raster(alpha,cm)

      # Pure maxxers in the same solid Figure-07 colours.
      draw_maxxer_masks(alpha,sm)

      abline(0,1,col="#000",lty=2,lwd=0.9)

      # Reference HDR contours: 50%, 90%, 99%.
      for(hc in hdr_curves){
        lines(hc[1,],hc[2,],col="#EEE",lwd=0.9)
      }

      # Figure-05 logic: exact 50:50 stationary curve, only where globally optimal.
      hc<-half_curve(
        alpha,B,
        sc$beta1,sc$beta2,sc$phi1,sc$phi2,
        scaler1,scaler2,eta,alim,nx
      )

      ii<-which(is.finite(hc$x) & hc$global)

      if(length(ii)>1){
        groups<-split(ii,cumsum(c(1,diff(ii)!=1)))

        for(r in groups){
          if(length(r)>1){
            lines(hc$x[r],hc$y[r],col=coldiag,lwd=3)
          }
        }
      }

      # Equal positive investments with some outside time. Global checks
      # exclude discontinuous switches between unequal allocations.
      equal_partial<-optout_equal_partial(alpha,c(sc$beta1,sc$beta2),
        c(sc$phi1,sc$phi2),c(scaler1,scaler2),eta,B,alim)
      for(branch in 1:2) lines(equal_partial[,branch],alpha,
        col=coldiag,lwd=1,lty=2)

      box()

      # EXACT Figure-07 population bar.
      draw_population_bar(pa$E0,pa$E1,pa$E2,alim)

      if(i==2){
        axis(1,at=seq(0,2.5,0.5),cex.axis=1.2)
        mtext(
          expression(alpha[1]~~"(talent in tournament 1)"),
          side=1,line=2.4,cex=0.95
        )
      }

      if(j==1){
        axis(2,at=seq(0,2.5,0.5),las=1,cex.axis=1.2)
        mtext(
          expression(alpha[2]~~"(talent in tournament 2)"),
          side=2,line=2.4,cex=0.95
        )
      }

      # EXACT Figure-05 two-line coloured parameter header.
      if(i==1){
        tx1<-alim[1]+0.30*diff(alim)
        tx2<-alim[1]+0.70*diff(alim)

        mtext(
          bquote(beta[1]==.(formatC(sc$beta1,format="f",digits=2))),
          side=3,line=1.55,at=tx1,cex=1.12,col=col1_max
        )
        mtext(
          bquote(beta[2]==.(formatC(sc$beta2,format="f",digits=2))),
          side=3,line=1.55,at=tx2,cex=1.12,col=col2_max
        )

        mtext(
          bquote(varphi[1]==.(formatC(sc$phi1,format="f",digits=2))),
          side=3,line=0.15,at=tx1,cex=1.12,col=col1_max
        )
        mtext(
          bquote(varphi[2]==.(formatC(sc$phi2,format="f",digits=2))),
          side=3,line=0.15,at=tx2,cex=1.12,col=col2_max
        )
      }

      # Only show a warning in a responsive panel if the active solver did not
      # actually converge. This is diagnostic text, not a strategy label.
      if(i==2 && !obj$eq$status%in%c("active","inactive")){
        text(
          alim[1]+0.04*diff(alim),
          alim[2]-0.04*diff(alim),
          obj$eq$status,
          adj=c(0,1),
          cex=0.62,
          col="white"
        )
      }
    }
  }

  # Row labels follow Figures 05/07.
  mtext(
    expression("Fixed "*E[1]==E[2]),
    side=2,outer=TRUE,at=0.76,line=1.25,cex=1.25
  )

  mtext(
    expression("Responsive "*E[t]),
    side=2,outer=TRUE,at=0.26,line=1.25,cex=1.25
  )

  # eta placement follows Figure 07.
  mtext(
    bquote(eta==.(formatC(eta,format="f",digits=2))),
    side=3,outer=TRUE,at=eta_at,line=0.50,cex=1.5,font=2
  )

  # EXACT Figure-07 legend panel margins and legend.
  par(mar=c(4.8,4.2,0.95,1.6),pty="m",xpd=FALSE)
  draw_allocation_legend()
}

optout_save(draw_asymmetric_optout_maps,"Figure_asymmetric_tournaments_optout_HDR",eta)
