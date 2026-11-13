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

