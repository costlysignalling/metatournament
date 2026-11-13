source('09_harvesting_model.R');source('08_finite_pool_model.R')
strategy_hdr<-function(){s<-boredom_slog();th<-seq(0,2*pi,length.out=1001);lapply(qchisq(c(.5,.9,.99),2,ncp=2*s*s),function(q)exp(-s*s+s*sqrt(q)*rbind(cos(th),sin(th))))}
strategy_polygon<-function(a,field,col){
 bands<-isoband::isobands(a,a,t(field),0,max(1,max(field)+1))
 if(length(bands[[1]]$x))plot(sf::st_sfc(isoband::iso_to_sfg(bands)[[1]]),add=TRUE,col=col,border=NA)
}
prepare_strategy_panel<-function(p,mode,a,finite=FALSE){
 actor<-as.matrix(expand.grid(a,a));alloc<-if(finite)function(z)finite_pool_alloc(p,mode,z)else function(z)harvest_alloc(p,mode,z)
 q<-matrix(alloc(actor),length(a));B<-p$B
 if(!finite){
  bd<-harvest_boundaries(p,mode,a);f1<-outer(a,bd[,2],'-');f2<- -outer(a,bd[,1],'-');eq<-bd[,3];eq[bd[,4]<.5]<-NA
  return(list(q=q,field1=f1,field2=f2,equal=list(list(x=eq,y=a))))
 }
 # Refine every transition along talent 1, including disconnected regions.
 rootfield<-function(target){
  mask<-q>=target;ij<-which(mask[-1,,drop=FALSE]!=mask[-nrow(mask),,drop=FALSE],arr.ind=TRUE)
  roots<-vector('list',length(a))
  if(nrow(ij)){
   lo<-a[ij[,1]];hi<-a[ij[,1]+1];y<-a[ij[,2]];side<-mask[cbind(ij[,1]+1,ij[,2])]
   for(k in 1:26){mid<-(lo+hi)/2;above<-alloc(cbind(mid,y))>=target;movehi<-above==side;hi[movehi]<-mid[movehi];lo[!movehi]<-mid[!movehi]}
   for(j in seq_along(a))roots[[j]]<-(lo+hi)[ij[,2]==j]/2
  }
  field<-matrix(0,length(a),length(a))
  for(j in seq_along(a))field[,j]<-(if(length(roots[[j]]))apply(abs(outer(a,roots[[j]],'-')),1,min)else rep(diff(range(a)),length(a)))*ifelse(mask[,j],1,-1)
  list(field=field,roots=roots)
 }
 one<-rootfield(B-1e-5);two<-rootfield(1e-5);half<-rootfield(B/2)
 equal<-list();nr<-max(lengths(half$roots))
 if(nr)for(k in 1:nr){
  xx<-vapply(half$roots,function(v)if(length(v)>=k)v[k]else NA_real_,0.);ii<-which(is.finite(xx))
  if(length(ii)){
   m<-if(mode=='first')p$reference else p$responsive$moment;bg<-(p$N-1)*c(m[1],B-m[1],m[2:3]);xbest<-alloc(cbind(xx[ii],a[ii]))
   value<-function(x){u<-boredom_w(x,xx[ii],p$beta[1],p$phi[1]);v<-boredom_w(B-x,a[ii],p$beta[2],p$phi[2]);(bg[1]+x)*u/(bg[3]+u)+(bg[2]+B-x)*v/(bg[4]+v)}
   xx[ii[value(B/2)<value(xbest)-1e-7*(1+B)]]<-NA
  }
  equal[[k]]<-list(x=xx,y=a)
 }
 list(q=q,field1=one$field,field2= -two$field,equal=equal)
}
draw_strategy_examples<-function(prep,stem,finite=FALSE,nalpha=361L,png=TRUE,pdf=TRUE,training=FALSE){
 for(pkg in c('sf','isoband'))if(!requireNamespace(pkg,quietly=TRUE))stop('Install ',pkg)
 a<-seq(0,2.5,length.out=nalpha);B<-50;panels<-lapply(prep,function(p)lapply(c('first','responsive'),function(m)if(training)prepare_training_comparison_panel(p,m,a)else prepare_strategy_panel(p,m,a,finite)))
 cols<-colorRampPalette(c('#315A9B','white','#B54343'))(801)
 draw<-function(){
  layout(matrix(c(1:4,9,5:8,9),2,byrow=TRUE),widths=c(1,1,1,1,.46),heights=c(1,1));par(oma=c(2,3.5,3,0))
  for(i in 1:2)for(j in 1:4){
   p<-prep[[j]];z<-panels[[j]][[i]];par(mar=c(2.75,2.75,.5,.2),mgp=c(1.75,.48,0),tcl=-.23,pty='s',xpd=FALSE)
   plot(NA,xlim=range(a),ylim=range(a),xaxs='i',yaxs='i',axes=FALSE,xlab='',ylab='',asp=1)
   cm<-matrix(cols[pmax(1,pmin(801,1+round(800*z$q/B)))],length(a))
   rasterImage(as.raster(t(cm)[ncol(cm):1,,drop=FALSE]),0,0,2.5,2.5,interpolate=TRUE)
   strategy_polygon(a,z$field2,'#08F');strategy_polygon(a,z$field1,'#E06')
   abline(0,1,lty=2,lwd=.9,col='black');for(h in strategy_hdr())lines(h[1,],h[2,],col='#EEE',lwd=.8)
   for(k in 1:2)for(cl in contourLines(a,a,if(k==1)z$field1 else z$field2,levels=0))lines(cl$x,cl$y,col=if(k==1)'#E06'else '#08F',lwd=1.05)
   for(eq in z$equal)lines(eq$x,eq$y,col='#FF0',lwd=2.8)
   box()
   energy<-if(finite){if(i==1)p$first$mapped[1]else p$responsive$mapped[1]}else{if(i==1)p$first[1]else p$responsive[1]}
   rect(2.54,0,2.69,2.5*energy/B,col='#E06',border=NA,xpd=NA);rect(2.54,2.5*energy/B,2.69,2.5,col='#08F',border=NA,xpd=NA);rect(2.54,0,2.69,2.5,xpd=NA,lwd=.7)
   if(i==2){axis(1,at=seq(0,2.5,.5),cex.axis=1.2);mtext(expression(alpha[1]~~'(talent in tournament 1)'),1,line=2.4,cex=.95)}
   if(j==1){axis(2,at=seq(0,2.5,.5),las=1,cex.axis=1.2);mtext(expression(alpha[2]~~'(talent in tournament 2)'),2,line=2.4,cex=.95)}
   if(i==1){
    mtext(LETTERS[j],3,line=1.8,at=0,adj=0,cex=1.25,font=2,xpd=NA)
    if(finite)mtext(bquote(N==.(p$N)),3,line=.15,at=0,adj=0,cex=1,xpd=NA)
    for(t in 1:2){at<-2.5*c(.48,.88)[t];color<-c('#E06','#08F')[t]
     fixed_type<-!finite&&!training&&(t==1||p$both_fixed)
     header<-if(fixed_type)bquote(P[.(t)]==alpha[.(t)])else bquote(beta[.(t)]==.(sprintf('%.2f',p$beta[t])))
     lhs<-if(fixed_type)bquote(P[.(t)])else bquote(beta[.(t)])
     rhs<-if(fixed_type)bquote(alpha[.(t)])else sprintf('%.2f',p$beta[t])
     # Align the equals signs, rather than the centres of unequal expressions.
     width<-function(expr)strwidth(expr,units='inches',cex=1.12)
     offset<-.5*(width(bquote(varphi[.(t)]))-width(sprintf('%.2f',p$phi[t]))-width(lhs)+width(rhs))
     header_at<-at+offset*diff(par('usr')[1:2])/par('pin')[1]
     mtext(header,3,line=1.55,at=header_at,cex=1.12,col=color)
     mtext(bquote(varphi[.(t)]==.(sprintf('%.2f',p$phi[t]))),3,line=.15,at=at,cex=1.12,col=color)
    }
   }
  }
  mtext(if(finite)expression('Fixed '*E[t])else expression('Fixed '*E[1]==E[2]),2,outer=TRUE,at=.76,line=1.25,cex=1.25)
  mtext(expression('Responsive '*E[t]),2,outer=TRUE,at=.26,line=1.25,cex=1.25)
  # Approved S14 strategy legend, reused verbatim (local colour names resolved).
  # Keep wording, ticks, font sizes, spacing and border style identical to S14.
  par(mar=c(2.5,0.55,0.5,1.60),pty="m",xpd=NA)
  plot.new()
  plot.window(xlim=c(-0.08,1.00),ylim=c(0,1),xaxs="i",yaxs="i")

  zz<-seq(0,1,length.out=801)

  for(k in 1:800){
    rect(0.34,zz[k],0.64,zz[k+1],col=cols[k],border=NA)
  }

  rect(0.34,0,0.64,1,border="black",lwd=0.8)

  wd<-0.006
  rect(0.34,0-wd,0.64,0+wd,col='#08F',lwd=0.8,xpd=NA,border=1)
  rect(0.34,0.5-wd,0.64,0.5+wd,col='#FF0',lwd=0.8,xpd=NA,border=1)
  rect(0.34,1-wd,0.64,1+wd,col='#E06',lwd=0.8,xpd=NA,border=1)

  ticks<-seq(0,1,0.25)
  segments(0.64,ticks,0.75,ticks,lwd=0.8)
  text(0.88,ticks,sprintf("%.2f",ticks),
       adj=c(0,0.5),cex=1.2,xpd=NA)

  text(0.10,0.5,"Relative time invested in tournament 1",
       srt=90,cex=1.7,xpd=NA)

  text(0.49,1.035,"1maxxers",
       cex=1.5,xpd=NA,col='#E06')

  text(0.49,-0.035,"2maxxers",
       cex=1.5,xpd=NA,col='#08F')

 }
 if(png){grDevices::png(paste0(stem,'.png'),width=3580/300,height=1740/300,units='in',res=600,type=if(capabilities('cairo'))'cairo'else getOption('bitmapType'));tryCatch(draw(),finally=dev.off())}
 if(pdf){cairo_pdf(paste0(stem,'.pdf'),width=3580/300,height=1740/300,family='Arial');tryCatch(draw(),finally=dev.off())}
 invisible(panels)
}
draw_all_harvesting_examples<-function(){
 d<-readRDS('09_harvesting_results/examples.rds')
 for(n in names(d)){message('Plot ',n);draw_strategy_examples(d[[n]],paste0('Figure_harvesting_examples_',n))}
}
