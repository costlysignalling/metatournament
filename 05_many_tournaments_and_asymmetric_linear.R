# Four-panel bridge from many identical tournaments to asymmetric tournaments.
# Run from this scripts folder. Requires only base R and the two supplied
# participation_many_summary.csv files; see REPRODUCING_PARTICIPATION_FIGURES.md
# for their calculation code. All strategy maps use analytic marginal-return
# conditions, with vectorized bisection for unequal training exponents.
# Panel D uses random allocations with E1=E2 and separately calculated M1,M2.
source("05_random_allocation_reference.R")
# Panel C is the exact alpha1=alpha2 slice. HDRs use a reference distribution
# in which the shared talent and alpha3 are independent lognormal variables
# with median 1 and original-scale SD 0.5; they are not masses in full 3D space.

B <- 50
talent_sdlog <- sqrt(log((1+sqrt(2))/2))
map_colours <- colorRampPalette(c('#315A9B','white','#B54343'))(801)
Ts <- 1:10
line_colours <- c('#444444',hcl.colors(9,'Dark 3'))

# Exact HDR contours for independent lognormal coordinate variables.
talent_hdr_curves <- function(probs=c(.5,.9,.99)) {
  s <- talent_sdlog
  th <- seq(0,2*pi,length.out=1001)
  lapply(qchisq(probs,df=2,ncp=2*s*s),function(q)
    exp(-s*s+s*sqrt(q)*rbind(cos(th),sin(th))))
}

# x is investment in tournament 1; both positive scalers are fixed.
linear_allocation <- function(a1,a2,beta1,beta2,B=50,scaler_ratio=1) {
  lo <- rep(0,length(a1)); hi <- rep(B,length(a1))
  for(i in 1:52) {
    mid <- (lo+hi)/2
    increase <- scaler_ratio*a1*beta1*(1+mid)^(beta1-1) >
      a2*beta2*(1+B-mid)^(beta2-1)
    lo[increase] <- mid[increase]; hi[!increase] <- mid[!increase]
  }
  x <- (lo+hi)/2
  x[scaler_ratio*a1*beta1<=a2*beta2*(1+B)^(beta2-1)] <- 0
  x[scaler_ratio*a1*beta1*(1+B)^(beta1-1)>=a2*beta2] <- B
  x[a1==0 & a2==0] <- NA_real_
  x
}

panel_letter <- function(letter) mtext(letter,side=3,adj=0,line=.85,font=2,cex=1.5)

draw_frequencies <- function(folder,letter,title,legend=FALSE) {
  d <- read.csv(file.path('Etude3_multiple_tournaments_analysis',folder,
    'participation_many_summary.csv'))
  d <- subset(d,beta==.8 & phi==1)
  d <- rbind(d[c('T','k','proportion')],data.frame(T=1,k=1,proportion=1))
  plot(NA,xlim=c(10.2,.8),ylim=c(-3,103),axes=FALSE,
    xaxs='i',yaxs='i',xlab='Tournaments entered (k)',ylab=if(letter=='A') 'Fraction of actors (%)' else '')
  abline(h=seq(0,100,20),col='grey85',lwd=.65)
  abline(v=1:10,col='grey94',lwd=.5)
  axis(1,at=10:1,cex.axis=.8,gap.axis=0)
  axis(2,at=seq(0,100,20),las=1,cex.axis=.85)
  for(t in Ts) {
    z <- d[d$T==t,]; z <- z[order(z$k,decreasing=TRUE),]
    lines(z$k,100*z$proportion,col=line_colours[t],lwd=1.25,
      type='o',pch=16,cex=.65)
  }
  box(); panel_letter(letter)
  mtext(title,side=3,line=.75,adj=1,cex=1.3)
  if(legend) legend('topleft',legend=paste('T =',Ts),col=line_colours,
    pch=16,lwd=1.1,ncol=2,bty='n',cex=.65,inset=c(.005,.015),
    x.intersp=.55,y.intersp=.85,seg.len=1.25)
}

draw_map <- function(z,alpha,letter,title,xlab,ylab,slopes,hdr=TRUE) {
  image(alpha,alpha,z,col=map_colours,zlim=c(0,1),xlim=c(0,2.5),ylim=c(0,2.5),
    axes=FALSE,xlab=xlab,ylab=ylab,xaxs='i',yaxs='i',useRaster=TRUE)
  image(alpha,alpha,ifelse(z==0,1,NA_real_),col='#08F',breaks=c(.5,1.5),
    add=TRUE,axes=FALSE,useRaster=TRUE)
  image(alpha,alpha,ifelse(z==1,1,NA_real_),col='#E06',breaks=c(.5,1.5),
    add=TRUE,axes=FALSE,useRaster=TRUE)
  if(hdr) for(curve in talent_hdr_curves()) lines(curve[1,],curve[2,],col='#EEE',lwd=.8)
  abline(0,1,col='black',lty=2,lwd=.75)
  abline(0,slopes[1],col='#08F',lwd=1.1)
  abline(0,slopes[3],col='#E06',lwd=1.1)
  abline(0,slopes[2],col='#FF0',lwd=3)
  axis(1,at=seq(0,2.5,.5),cex.axis=.8)
  axis(2,at=seq(0,2.5,.5),las=1,cex.axis=.85)
  box(); panel_letter(letter)
  mtext(title,side=3,line=.75,adj=1,cex=1.3)
}

draw_composite_slice <- function(alpha) {
  grid <- expand.grid(pair=alpha,third=alpha)
  ratio <- (grid$pair/grid$third)^(1/(1-.8))
  z <- (ratio*(B+1)-1)/(ratio+.5)
  z[is.infinite(ratio)] <- B
  z <- pmax(0,pmin(B,z))/B
  z[grid$pair==0 & grid$third==0] <- NA_real_
  # y/x slopes: third-only (blue), 50:50 pair/third, pair-only (pink).
  slopes <- c((1+B)^(.2),((1+B/2)/(1+B/4))^(.2),(1+B/2)^(-.2))
  draw_map(matrix(z,length(alpha)),alpha,'C',expression(beta==0.80),
    expression(alpha[1]==alpha[2]),expression(alpha[3]~"(talent in tournament 3)"),slopes)
}

# Change panel D's two training exponents here. Its population reference,
# allocation map, boundary rays, equal-allocation ray and title all follow.
draw_asymmetric <- function(alpha,beta1=.8,beta2=.4) {
  grid <- expand.grid(a1=alpha,a2=alpha)
  ref <- random_allocation_reference(beta1,beta2,1,1,B)
  z <- linear_allocation(grid$a1,grid$a2,beta1,beta2,B,scaler_ratio=ref$s)/B
  # alpha1/alpha2 at: 2maxxing boundary, equal investment, 1maxxing boundary.
  rho <- (beta2/(ref$s*beta1))*c((1+B)^(beta2-1),
    (1+B/2)^(beta2-beta1),(1+B)^(1-beta1))
  draw_map(matrix(z,length(alpha)),alpha,'D',bquote(beta[1]==.(formatC(beta1,format='f',digits=2))*','~~
           beta[2]==.(formatC(beta2,format='f',digits=2))),
           expression(alpha[1]~"(talent in tournament 1)"),expression(alpha[2]~"(talent in tournament 2)"),1/rho)
}

# Draw all four panels; panel C uses equal paired talents.
# The reduced device height is matched by reduced vertical margins, preserving
# panel sizes and typography while trimming whitespace above and below.
draw_bridge <- function() {
  par(mfrow=c(1,4),mar=c(3.1,3.2,2.375,.35),oma=c(.475,0,0,.4),
    mgp=c(1.7,.40,0),tcl=-.22,cex=1.25,pty='s')
  draw_frequencies('T10_phi125','A',expression(B==50),TRUE)
  draw_frequencies('proportional_budget','B',expression(B[T]==25*T))
  draw_composite_slice(seq(0,2.5,length.out=401))
  draw_asymmetric(seq(0,2.5,length.out=401))
}

png('Figure_many_tournaments_and_asymmetric_linear.png',
    width=13.2,height=3.95,units='in',res=600)
draw_bridge()
dev.off()
# To export as vector PDF instead, replace png() above with:
# pdf('Figure_many_tournaments_and_asymmetric_linear.pdf',width=13.2,height=3.95)
