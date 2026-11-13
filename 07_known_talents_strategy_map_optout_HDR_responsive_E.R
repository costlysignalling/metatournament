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
nx<-501

#Baseline talent distribution
median<-1
sd<-0.5
m<-log(median)
slog<-sqrt(log((1+sqrt(1+4*(sd/median)^2))/2))

# The responsive population is set by metatournament.optout.npop (default 8192).

nxpop<-501

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
#Global allocation with self-directed time
#------------------------------------------------------------

#For each row, return the running maximum and its first location.
#This lets us optimize V1(x1)+V2(x2) subject to x1+x2<=B without
#searching the full x1-x2 triangle separately for every talent pair.
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
#Responsive common scaler
#------------------------------------------------------------


# The shared solver uses continuous best responses and a deterministic QMC
# population; the same population supplies the equilibrium and its budget bar.
responsive_scaler<-function(beta,phi,eta,B,nx=501,...){
  eq<-optout_responsive(rep(beta,2),rep(phi,2),eta,B)
  st<-eq$state; k<-eq$scaler[1]
  c(scaler=k,scaler0=optout_fixed_scalers(beta,phi,B)$scaler[1],
    lambda=if(k>0) eta/k else Inf,E0=st$E0,E1=st$E1,E2=st$E2,
    M1=st$M1,M2=st$M2,Utotal=st$Utotal,Uidentity=B-(1-eta)*st$E0,
    collapsed=as.numeric(eq$status=='inactive'),npop=nrow(optout_population()),residual=eq$err)
}
#Analytic/numerical boundaries of the three pure maxxer regions.
#
#Self-maxxing is optimal iff each tournament separately fails to beat the
#self-directed return. This gives the scalar talent threshold
#
#  alpha_self = min_x {[(1+eta*x/scaler)^(1/phi)-1] / g(x)}.
#
#For a 1-maxxer, start from x1=B. If x2 units are moved to tournament 2,
#the best amount that can remain in tournament 1 is the prefix maximum of
#V1 over [0,B-x2]. Requiring every such deviation to be unprofitable gives,
#for each alpha1, the largest alpha2 compatible with full maxxing in 1.
#The 2-maxxer boundary is its mirror image.
maxxer_boundaries<-function(alpha_boundary,B,beta,phi,scaler,eta,nxb=4001){
  #Inactive-tournament boundary state: every talent combination is a self-maxxer.
  if(scaler<=0){
    return(list(alpha=alpha_boundary,self=Inf,
                max_other=rep(NA_real_,length(alpha_boundary))))
  }

  xx<-seq(0,B,length.out=nxb)
  gx<-g(xx,beta)

  #Complete self-maxxer threshold. x=0 is excluded because both numerator
  #and denominator vanish there.
  xpos<-xx[-1]
  self_alpha<-min(eta/(scaler*phi*beta),
    ((1+eta*xpos/scaler)^(1/phi)-1)/g(xpos,beta))

  other_max<-rep(NA_real_,length(alpha_boundary))
  cap<-nxb:1

  for(i in seq_along(alpha_boundary)){
    a1<-alpha_boundary[i]
    v1<-scaler*w(xx,a1,beta,phi)-eta*xx
    v1B<-v1[nxb]

    #If x1=B is not itself a global one-tournament optimum, there can be
    #no 1-maxxer region for this alpha1.
    if(v1B < max(v1)-1e-11) next

    p1<-cummax(v1)
    x2<-xx[-1]
    best1<-p1[cap[-1]]
    gap<-v1B-best1

    #Solve V2(x2,alpha2)=gap for alpha2 and impose the condition for every
    #possible x2. The most restrictive deviation determines the boundary.
    rhs<-1+(gap+eta*x2)/scaler
    aa<-rep(Inf,length(x2))
    ok<-rhs>0 & gx[-1]>0
    aa[ok]<-(rhs[ok]^(1/phi)-1)/gx[-1][ok]
    other_max[i]<-max(0,min(aa,na.rm=TRUE))
  }

  list(alpha=alpha_boundary,self=self_alpha,max_other=other_max)
}

#Classify grid points that are already known analytically to be pure maxxers.
#Any overlap of the two tournament-maxxer inequalities is a tie. Off the
#diagonal it is assigned to the higher-talent tournament; exact diagonal ties
#are left to the numerical optimizer.
maxxer_masks<-function(alpha,bounds){
  n<-length(alpha)
  A1<-matrix(alpha,n,n)
  A2<-t(A1)

  ok<-is.finite(bounds$max_other)
  b<-rep(NA_real_,n)
  if(sum(ok)>=2){
    inside<-alpha>=min(bounds$alpha[ok]) & alpha<=max(bounds$alpha[ok])
    b[inside]<-approx(bounds$alpha[ok],bounds$max_other[ok],
                      xout=alpha[inside],rule=1,ties="ordered")$y
  }
  B1<-matrix(b,n,n)

  self<-A1<=bounds$self & A2<=bounds$self
  one<-!self & is.finite(B1) & A2<=B1
  two<-!self & t(one)

  both<-one & two
  one[both & A1<=A2]<-FALSE
  two[both & A2<=A1]<-FALSE

  list(self=self,one=one,two=two)
}

#Returns globally optimal shares x0/B, x1/B and x2/B throughout talent space.
#Pure self/1/2-maxxer regions are filled directly from their boundaries; the
#discrete budget-grid optimization is performed only on the remaining points.
#The transformed objective is
#  V1(x1)+V2(x2),
#where Vt(x)=scaler*w(x,alpha)-eta*x and x1+x2<=B.
strategy_map_optout<-function(alpha,B,beta,phi,scaler,eta,bounds,nx=501){
  sm<-optout_talent_map(alpha,rep(beta,2),rep(phi,2),rep(scaler,2),eta,B)
  sm$bounds<-bounds;sm
}
#Draw the pure maxxer regions as vector polygons. Their edges therefore follow
#the calculated talent-space boundaries rather than the strategy-map raster.
draw_maxxer_regions<-function(bounds,alim){
  x<-bounds$alpha
  y<-pmin(bounds$max_other,alim[2])
  valid<-is.finite(y) & y>=alim[1] & x>=alim[1] & x<=alim[2]

  idx<-which(valid)
  if(length(idx)){
    groups<-split(idx,cumsum(c(1,diff(idx)!=1)))
    for(ii in groups){
      if(length(ii)<2) next
      xx<-x[ii]
      yy<-y[ii]
      px<-c(xx,rev(xx))
      py<-c(rep(alim[1],length(ii)),rev(yy))
      polygon(px,py,col=col1_max,border=NA)
      polygon(py,px,col=col2_max,border=NA)
    }
  }

  if(bounds$self>alim[1]){
    a<-min(bounds$self,alim[2])
    rect(alim[1],alim[1],a,a,col=self_max,border=NA)
  }
}

#------------------------------------------------------------
#Population mean allocation
#------------------------------------------------------------

#Population-average allocation under the baseline talent distribution.
#
# Bars show first-response or equilibrium means from the same deterministic population.
population_summary_optout<-function(beta,phi,scaler,eta,B,nx=801,npop=201){
  st<-optout_state(rep(scaler,2),rep(beta,2),rep(phi,2),eta,B)
  c(E0=st$E0,E1=st$E1,E2=st$E2,M1=st$M1,M2=st$M2,
    Utour=scaler*(st$M1+st$M2),Uself=eta*st$E0,Utotal=st$Utotal,
    tournament_sink=B-scaler*(st$M1+st$M2))
}
#------------------------------------------------------------
#Talent HDR contours
#------------------------------------------------------------

#Exact iso-density curves enclosing the requested probability mass for
#the bivariate lognormal baseline talent distribution.
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
#Colours
#------------------------------------------------------------

#When x0=0 this reproduces the old blue-white-red diversifier gradient.
#As x0/B increases, the colour is blended toward self_col. Complete
#self-maxxers and the two tournament maxxer types are then overlaid in solid
#colours, and the global 50:50 splitter line is highlighted in yellow.
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

#Draw a raster whose matrix coordinates follow image(x=alpha,y=alpha,z=...).
draw_colour_raster<-function(alpha,colmat){
  ras<-as.raster(t(colmat)[ncol(colmat):1,,drop=FALSE])
  rasterImage(ras,min(alpha),min(alpha),max(alpha),max(alpha),interpolate=FALSE)
}

#Stacked side bar showing population-average budget allocation.
#Tournament 1 is at the bottom, self-directed time in the middle,
#and tournament 2 at the top, matching the side bars in Figure 05.
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

  #Self-directed time at the bottom, followed by tournament 1 and tournament 2.
  rect(x0,y0,x1,y1,col=col0,border=NA,xpd=NA)
  rect(x0,y1,x1,y2,col=col1,border=NA,xpd=NA)
  rect(x0,y2,x1,y3,col=col2,border=NA,xpd=NA)
  rect(x0,y0,x1,y3,border="black",col=NA,lwd=0.7,xpd=NA)

  xm<-(x0+x1)/2

  #Labels are centered in their population-share rectangles. Very thin
  #segments are left unlabeled rather than forcing overlapping text.
  if(y1-y0>0.055*dy)
    text(xm,(y0+y1)/2,expression(E[0]),cex=0.67,xpd=NA,col="white",font=2)
  if(y2-y1>0.055*dy)
    text(xm,(y1+y2)/2,expression(E[1]),cex=0.67,xpd=NA,col="white",font=2)
  if(y3-y2>0.055*dy)
    text(xm,(y2+y3)/2,expression(E[2]),cex=0.67,xpd=NA,col="white",font=2)
}

#Two-dimensional legend: horizontal = self-directed share; vertical =
#relative tournament-1 share among the time that remains in tournaments.
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
  axis(2,at=seq(0,1,0.1),labels=formatC(seq(0,1,0.1),format="f",digits=1),
       las=1,cex.axis=0.9,tcl=-0.15)

  wd<-0.15
  hd<-0.018
  rect(0,1-hd,wd,1,col=col1_max,border=1,lwd=0.8,xpd=NA)
  rect(0,0.5-hd/2,wd,0.5+hd/2,border=1,col="#FF0",lwd=0.8,xpd=NA)
  rect(0,0,wd,hd,border=1,col=col2_max,lwd=0.8,xpd=NA)

  rect(1,0,1+wd*0.8,1,col=self_max,border=1,lwd=0.8,xpd=NA)

  mtext(expression(x[0]/B),
        side = 1, line = 2.2, cex = 0.9)
    mtext("(self-directed",
        side = 1, line = 3.25, cex = 0.9)
    mtext("share)",
        side = 1, line = 4.3, cex = 0.9)
  
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
#Prepare maps
#------------------------------------------------------------

alpha<-seq(alim[1],alim[2],length.out=nalpha)
alpha_boundary<-seq(alim[1],alim[2],length.out=1201)
hdr_curves<-talent_hdr_curves(median,sd,probs=c(0.50,0.90,0.99))

eq_state<-vector("list",length(betas)*length(phis))
dim(eq_state)<-c(length(betas),length(phis))

maps<-vector("list",length(betas)*length(phis))
dim(maps)<-c(length(betas),length(phis))

pop_summary<-vector("list",length(betas)*length(phis))
dim(pop_summary)<-c(length(betas),length(phis))

for(i in seq_along(betas)){
  for(j in seq_along(phis)){
    eq_state[[i,j]]<-responsive_scaler(
      betas[i],phis[j],eta,B,nxpop
    )
    scaler<-unname(eq_state[[i,j]]["scaler"])

    bounds<-maxxer_boundaries(alpha_boundary,B,betas[i],phis[j],
                              scaler,eta,nxb=4001)
    maps[[i,j]]<-strategy_map_optout(alpha,B,betas[i],phis[j],
                                     scaler,eta,bounds,nx)

    # The displayed bar uses the same QMC population as the equilibrium calculation.
    pop_summary[[i,j]]<-population_summary_optout(
      betas[i],phis[j],scaler,eta,B,nxpop
    )
  }
}

cat("\nResponsive population equilibria:\n")
for(i in seq_along(betas)){
  for(j in seq_along(phis)){
    eq<-eq_state[[i,j]]
    pa<-pop_summary[[i,j]]
    cat("beta =",betas[i]," phi =",phis[j],
        " scaler0 =",formatC(eq["scaler0"],digits=4,format="f"),
        " scaler* =",formatC(eq["scaler"],digits=4,format="f"),
        " lambda* =",formatC(eq["lambda"],digits=4,format="f"),
        " state =",ifelse(eq["collapsed"]>0.5,"collapsed","active"),
        " QMC actors =",formatC(eq["npop"],digits=0,format="f"),
        " E0/B =",formatC(pa["E0"]/B,digits=3,format="f"),
        " Ubar/B =",formatC(1-(1-eta)*pa["E0"]/B,digits=3,format="f"),
        "\n")
  }
}

#------------------------------------------------------------
#Figure
#------------------------------------------------------------

draw_strategy_map_optout<-function(){
  legw<-0.52
  widths<-c(1,1,1,1,legw)
  layout(matrix(c(1,2,3,4,9,
                  5,6,7,8,9),nrow=2,byrow=TRUE),
         widths=widths,heights=c(1,1))

  eta_at<-(sum(widths[1:4])+widths[5]/2)/sum(widths)

  par(oma=c(2,3.5,3,0.2))

  for(i in seq_along(betas)){
    for(j in seq_along(phis)){
      sm<-maps[[i,j]]
      cm<-allocation_colours(sm$x0,sm$x1,sm$x2)

      par(mar=c(2.75,2.75,0.5,0.2),
          mgp=c(1.75,0.48,0),tcl=-0.23,pty="s",xpd=FALSE)

      plot(NA,xlim=alim,ylim=alim,xaxs="i",yaxs="i",
           axes=FALSE,xlab="",ylab="",asp=1)
      draw_colour_raster(alpha,cm)

      #Pure maxxer regions are exact/smooth polygons, independent of the
      #talent-grid raster used for mixed allocations.
      draw_maxxer_regions(sm$bounds,alim)

      tol<-0.5/(nx-1)
      is_half<-(sm$x0<=tol & abs(sm$x1-0.5)<=tol & abs(sm$x2-0.5)<=tol)

      abline(0,1,col="#000",lty=2,lwd=0.9)

      #Reference HDR contours: 50%, 90%, 99%.
      for(hc in hdr_curves){
        lines(hc[1,],hc[2,],col="#EEE",lwd=0.9)
      }

      splitidx<-which(diag(is_half))
      if(length(splitidx)>1){
        groups<-split(splitidx,cumsum(c(1,diff(splitidx)!=1)))
        for(r in groups){
          if(length(r)>1){
            lines(alpha[r],alpha[r],col=coldiag,lwd=3)
          }
        }
      }

      box()

      #Population-average allocation under the baseline talent distribution.
      pa<-pop_summary[[i,j]]
      draw_population_bar(pa["E0"],pa["E1"],pa["E2"],alim)

      #The common scaler is the responsive population fixed point.
      # text(alim[2]-0.06*diff(alim),alim[2]-0.06*diff(alim),
      #      bquote(scaler^"*"==.(formatC(eq_state[[i,j]]["scaler"],
      #                                format="f",digits=2))),
      #      adj=c(1,1),cex=0.75)
      #This $\scaler$ tells us how much utility one unit of raw spoils weight is worth under the normalization inherited from the corresponding no-opt-out population.
      #$\eta$ is being compared to tournament utility after this panel-specific normalization.

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
        mtext(bquote(varphi==.(formatC(phis[j],format="f",digits=2))),
              side=3,line=0.95,cex=1.4)
      }
    }
  }

  mtext(expression(beta==0.4),side=2,outer=TRUE,at=0.76,line=1.25,cex=1.4)
  mtext(expression(beta==0.8),side=2,outer=TRUE,at=0.26,line=1.25,cex=1.4)
  mtext(bquote(eta==.(formatC(eta,format="f",digits=2))),
        side=3,outer=TRUE,at=eta_at,line=0.50,cex=1.5,font=2)
  mtext(expression("Responsive "*E[t]),
        side=3,outer=TRUE,at=-0.03,line=0.5,cex=1.15,adj=0,col="#082")

  par(mar=c(4.8,4.2,0.95,1.6),pty="m",xpd=FALSE)
  draw_allocation_legend()
}

optout_save(draw_strategy_map_optout,"Figure_known_talent_optout_strategy_map_HDR_responsive_E",eta)
