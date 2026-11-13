talent_density<-function(median=1,sd=0.5,rho=0,
                         probs=c(0.50,0.90,0.99),contcols=rep("white",4),
                         offsx=c(0.09,0.14,0.16),offsy=rep(-0.15,3),
                         alim=NULL,n=501,
                         margcol="grey80",
                         statcols=c("dodgerblue3","darkorange2","firebrick2")){
  
  m<-log(median)
  s<-sqrt(log((1+sqrt(1+4*(sd/median)^2))/2))
  
  mu<-c(m,m)
  Sigma<-s^2*matrix(c(1,rho,rho,1),2,2)
  Sinv<-solve(Sigma)
  
  one<-c(1,1)
  centre<-mu-Sigma%*%one
  ncp<-drop(t(one)%*%Sigma%*%one)
  
  qs<-qchisq(probs,df=2,ncp=ncp)
  
  #Exact transformed density contours
  th<-seq(0,2*pi,length.out=1001)
  L<-t(chol(Sigma))
  
  curves<-lapply(qs,function(q){
    y<-drop(centre)+L%*%rbind(sqrt(q)*cos(th),sqrt(q)*sin(th))
    rbind(exp(y[1,]),exp(y[2,]))
  })
  
  if(is.null(alim)){
    mx<-max(vapply(curves,max,1))
    alim<-c(0,1.05*mx)
  }
  
  alpha<-seq(max(alim[1],1e-5),alim[2],length.out=n)
  
  a1<-matrix(alpha,n,n)
  a2<-matrix(alpha,n,n,byrow=TRUE)
  
  y1<-log(a1)
  y2<-log(a2)
  d1<-y1-m
  d2<-y2-m
  
  Q<-Sinv[1,1]*d1^2+
    2*Sinv[1,2]*d1*d2+
    Sinv[2,2]*d2^2
  
  dens<-exp(-Q/2)/(2*pi*sqrt(det(Sigma))*a1*a2)
  
  par(mar=c(3.5,3.5,2.5,2.5),
      mgp=c(2,0.65,0),tcl=-0.23,pty="s")
  
  cols<-hcl.colors(256,"Viridis")
  image(alpha,alpha,dens,col=cols,
        xlim=alim,ylim=alim,
        xaxs="i",yaxs="i",
        xlab=expression(alpha[1]~~"(talent in tournament 1)"),
        ylab=expression(alpha[2]~~"(talent in tournament 2)"),
        useRaster=TRUE)
  
  for(i in rev(seq_along(curves))){
    z<-curves[[i]]
    lines(z[1,],z[2,],lwd=1.2,col=contcols[i])
    k<-which.max(z[1,]+z[2,])
    text(z[1,k]+offsx[i],z[2,k]+offsy[i],
         paste0(round(100*probs[i]),"%"),
         pos=4,cex=0.85,xpd=NA,col=contcols[i])
  }
  
  points(exp(centre[1]),exp(centre[2]),pch=3,cex=1.1,lwd=1.3,
         col=contcols[length(probs)+1])
  
  #Marginal LN density, rescaled just for plotting
  mdens<-dlnorm(alpha,m,s)
  mdens<-mdens/max(mdens)*0.15*diff(alim)
  
  #Top marginal (alpha_1)
  polygon(c(alpha,rev(alpha)),
          c(alim[2]+mdens,rep(alim[2],length(alpha))),
          col=margcol,border=NA,xpd=NA)
  
  #Right marginal (alpha_2)
  polygon(c(alim[2]+mdens,rep(alim[2],length(alpha))),
          c(alpha,rev(alpha)),
          col=margcol,border=NA,xpd=NA)
  
  #Analytical summaries of LN(m,s)
  xmode<-exp(m-s^2)
  xmed <-exp(m)
  xmean<-exp(m+s^2/2)
  xstats<-c(xmode,xmed,xmean)
  
  #Markers: top marginal
  for(i in 1:3){
    h<-dlnorm(xstats[i],m,s)/max(dlnorm(alpha,m,s))*0.15*diff(alim)
    segments(xstats[i],alim[2],
             xstats[i],alim[2]+h,
             col=statcols[i],lwd=1.4,xpd=NA)
  }
  
  #Markers: right marginal
  for(i in 1:3){
    h<-dlnorm(xstats[i],m,s)/max(dlnorm(alpha,m,s))*0.15*diff(alim)
    segments(alim[2],xstats[i],
             alim[2]+h,xstats[i],
             col=statcols[i],lwd=1.4,xpd=NA)
  }
  
  #Labels inside marginal plots
  #Labels inside marginal plots, slightly before the mode
  #Label positions inside marginal densities
  #1st value = shift along the axis, relative to the mode (negative = before mode)
  #2nd value = how deep into the grey shape the label sits, in units of local height
  #Label positions: distance from mode, height within density
  lab1<-c(-0.05,0.3)  #alpha_1
  lab2<-c(-0.05,0.3)  #alpha_2
  
  maxdens<-max(dlnorm(alpha,m,s))
  labpos<-function(lab){
    p<-xmode+lab[1]*diff(alim)
    h<-dlnorm(p,m,s)/maxdens*0.15*diff(alim)
    c(p,alim[2]+lab[2]*h)
  }
  
  pos1<-labpos(lab1)
  pos2<-labpos(lab2)
  
  text(pos1[1],pos1[2],expression(alpha[1]),cex=0.8,xpd=NA)
  text(pos2[2],pos2[1],expression(alpha[2]),cex=0.8,srt=270,xpd=NA)
  
  #Legend top-right outside
  legend(alim[2]-0.17*diff(alim),alim[2]+0.2*diff(alim),
         legend=c("mode","median","mean"),
         col=statcols,lwd=1.6,bty="n",xpd=NA,cex=0.8)
  
  box()
  
  invisible(list(
    m=m,s=s,Sigma=Sigma,
    centre=drop(centre),
    probabilities=probs,
    q=qs,
    curves=curves,
    mode=xmode,median=xmed,mean=xmean
  ))
}

talent_density()

w<-4
png("FigureS02.png",width=w,height=w,units="in",res=600)
talent_density()
rect(0,0,2.5,2.5,lty=2,border="white")
dev.off()
