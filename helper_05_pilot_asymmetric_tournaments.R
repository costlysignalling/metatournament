# Fixed E1=E2: uniform random allocations independent of talent, E1=E2=B/2.
# Each population mean spoils weight M_t is calculated separately.
# Responsive E_t updates the allocation distribution and both E_t and M_t.
source("05_random_allocation_reference.R")
source("05_linear_population_means.R")
pilot_asymmetric<-function(beta1,beta2,phi1,phi2,
                           mode=c("fixed","responsive"),
                           alim=c(0,3.5),
                           hdr_probs=c(0.50,0.90,0.99),
                           hdr_col="black"){
  mode<-match.arg(mode)
  alpha<-seq(alim[1],alim[2],length.out=nalpha)
  
  cache<-population_cache(beta1,beta2,phi1,phi2,B,nxpop)
  
  if(mode=="fixed"){
    st<-random_allocation_reference(beta1,beta2,phi1,phi2,B,m,slog)
    s_use<-st$s
    response<-if(phi1==1 && phi2==1) linear_population_means(beta1,beta2,s_use,B,m,slog) else population_state(s_use,cache,B)
    E1<-response$E1
    E2<-response$E2
  } else {
    eq<-equilibrium_s(cache,B)
    s_use<-unname(eq["s"])
    E1<-unname(eq["E1"])
    E2<-unname(eq["E2"])
  }
  
  sm<-strategy_map(alpha,B,beta1,beta2,phi1,phi2,s_use,nx,nroot)

  # Highest-density contours of the joint talent distribution.
  # With independent lognormal talents, equal-density contours are circles
  # in log(alpha) space centred at (m-slog^2, m-slog^2), then exponentiated.
  hdr_q<-qchisq(hdr_probs,df=2,ncp=2*slog^2)
  hdr_theta<-seq(0,2*pi,length.out=1001)
  hdr_centre<-m-slog^2
  hdr_curves<-lapply(hdr_q,function(q){
    rr<-slog*sqrt(q)
    rbind(
      exp(hdr_centre+rr*cos(hdr_theta)),
      exp(hdr_centre+rr*sin(hdr_theta))
    )
  })
  
  par(mar=c(4,4,3.3,1.3),
      mgp=c(2.4,0.7,0),tcl=-0.23,pty="s",xpd=FALSE)
  
  image(alpha,alpha,sm$z,col=cols,zlim=c(0,1),
        xlim=alim,ylim=alim,xaxs="i",yaxs="i",
        axes=FALSE,xlab="",ylab="",useRaster=TRUE)
  
  z2<-ifelse(sm$max2,1,NA_real_)
  z1<-ifelse(sm$max1,1,NA_real_)
  
  image(alpha,alpha,z2,col=col2_max,breaks=c(0.5,1.5),
        add=TRUE,axes=FALSE,useRaster=TRUE)
  image(alpha,alpha,z1,col=col1_max,breaks=c(0.5,1.5),
        add=TRUE,axes=FALSE,useRaster=TRUE)
  
  lines(sm$b2,alpha,col=adjustcolor(col2_max,0.9),lwd=1.05)
  lines(alpha,sm$b1,col=adjustcolor(col1_max,0.9),lwd=1.05)
  
  hc<-half_curve(alpha,B,beta1,beta2,phi1,phi2,s_use,alim,nx)
  
  ii<-which(is.finite(hc$x) & hc$global)
  if(length(ii)>1){
    groups<-split(ii,cumsum(c(1,diff(ii)!=1)))
    for(r in groups){
      if(length(r)>1){
        lines(hc$x[r],hc$y[r],col=coldiag,lwd=3)
      }
    }
  }

  # Talent-density / HDR contours. The outer 99% contour is emphasized.
  hdr_lty<-rep(2,length(hdr_probs))
  hdr_lwd<-rep(1.15,length(hdr_probs))
  hdr_lty[length(hdr_probs)]<-1
  hdr_lwd[length(hdr_probs)]<-1.8

  for(k in seq_along(hdr_curves)){
    z<-hdr_curves[[k]]
    lines(z[1,],z[2,],col=hdr_col,lty=hdr_lty[k],lwd=hdr_lwd[k])

    # Put a small percentage label near the upper-right part of each contour.
    inside<-which(z[1,]>=alim[1] & z[1,]<=alim[2] &
                  z[2,]>=alim[1] & z[2,]<=alim[2])
    if(length(inside)){
      kk<-inside[which.max(z[1,inside]+z[2,inside])]
      text(z[1,kk],z[2,kk],paste0(round(100*hdr_probs[k]),"%"),
           pos=2,offset=0.25,cex=.72,col=hdr_col,xpd=FALSE)
    }
  }
  
  box()
  draw_budget_bar(E1,E2,alim)
  
  axis(1,at=seq(alim[1],alim[2],0.5),cex.axis=1.2)
  axis(2,at=seq(alim[1],alim[2],0.5),las=1,cex.axis=1.2)
  
  mtext(expression(alpha[1]~~"(talent in tournament 1)"),
        side=1,line=2.4,cex=0.95)
  mtext(expression(alpha[2]~~"(talent in tournament 2)"),
        side=2,line=2.4,cex=0.95)
  
  tx1<-alim[1]+0.30*diff(alim)
  tx2<-alim[1]+0.70*diff(alim)
  
  mtext(bquote(beta[1]==.(formatC(beta1,format="f",digits=2))),
        side=3,line=1.55,at=tx1,cex=1.12,col=col1_max)
  mtext(bquote(beta[2]==.(formatC(beta2,format="f",digits=2))),
        side=3,line=1.55,at=tx2,cex=1.12,col=col2_max)
  
  mtext(bquote(varphi[1]==.(formatC(phi1,format="f",digits=2))),
        side=3,line=0.15,at=tx1,cex=1.12,col=col1_max)
  mtext(bquote(varphi[2]==.(formatC(phi2,format="f",digits=2))),
        side=3,line=0.15,at=tx2,cex=1.12,col=col2_max)
  
  if(mode=="fixed"){
    mtext(expression("Fixed "*E[1]==E[2]),side=2,line=4,cex=0.95)
  } else {
    mtext(expression("Responsive "*E[t]),side=2,line=4,cex=0.95)
  }
  
  abline(0,1,lty=2)
  
  invisible(list(
    s=s_use,E1=E1,E2=E2,map=sm,
    hdr_probs=hdr_probs,hdr_curves=hdr_curves
  ))
}


pilot_asymmetric(0.01,0.6,0.5,0.5,"fixed")
pilot_asymmetric(0.01,0.6,0.5,0.5,"responsive")

# example_combinations<-list(
#   c(beta1=0.5,beta2=0.2,phi1=1.5,phi2=1.5),
#   c(beta1=0.8,beta2=0.8,phi1=1.5,phi2=0.66),
#   c(beta1=0.8,beta2=0.5,phi1=0.66,phi2=1),
#   c(beta1=0.8,beta2=0.2,phi1=1.5,phi2=0.66)
# )

pilot_asymmetric(0.5,0.2,1.5,1.5,"fixed")
pilot_asymmetric(0.5,0.2,1.5,1.5,"responsive")

pilot_asymmetric(0.8,0.8,1.5,0.66,"fixed")
pilot_asymmetric(0.8,0.8,1.5,0.66,"responsive")

pilot_asymmetric(0.8,0.5,0.66,1,"fixed")
pilot_asymmetric(0.8,0.5,0.66,1,"responsive")

pilot_asymmetric(0.8,0.2,1.5,0.66,"fixed")
pilot_asymmetric(0.8,0.2,1.5,0.66,"responsive")




#Tady pokusy:
#Different phis
pilot_asymmetric(0.6,0.6,1,1.5,"fixed")
pilot_asymmetric(0.6,0.6,1,1.5,"responsive")

pilot_asymmetric(0.6,0.6,0.5,0.9,"fixed")
pilot_asymmetric(0.6,0.6,0.5,0.9,"responsive")

pilot_asymmetric(0.6,0.6,1,1.5,"fixed")
pilot_asymmetric(0.6,0.6,1,1.5,"responsive")

#Different betas
pilot_asymmetric(0.6,0.8,1,1,"fixed")
pilot_asymmetric(0.6,0.8,1,1,"responsive")


