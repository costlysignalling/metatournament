# Match the selected-example layout of main-text Figure 9.
# For Figure 10 and its utility contours (S33), use draw_paper_boredom() in
# 08_prepare_paper_boredom.R. Redraws use saved data, without compiling.
# The original example set remains available; more_boredom now draws S36.
boredom_example_colours<-function(q,B,c){
 d<-pmax(c*B-q,0);z2<-pmin((1-c)*B,B-q);r<-q/pmax(q+z2,1e-15)
 vv<-colorRamp(c('#315A9B','white','#B54343'))(as.vector(r))/255
 dd<-as.vector(d)/B;vv<-vv*(1-dd)+(85/255)*dd
 out<-rgb(vv[,1],vv[,2],vv[,3]);dim(out)<-dim(q);out[q>=B-1e-5]<-'#E06';out
}
boredom_hdr<-function(){
 s<-sqrt(log((1+sqrt(2))/2));Sigma<-diag(s*s,2);centre<- -rowSums(Sigma);theta<-seq(0,2*pi,length.out=1001)
 lapply(qchisq(c(.5,.9,.99),df=2,ncp=2*s*s),function(q)exp(centre+t(chol(Sigma))%*%rbind(sqrt(q)*cos(theta),sqrt(q)*sin(theta))))
}
plot_boredom_examples<-function(compulsory=.5,N=30,png=TRUE,pdf=TRUE,variant=c('original','more_boredom'),data=NULL,contours=FALSE,paper=FALSE,stem=NULL){
 variant<-match.arg(variant);tag<-paste0('c_',format(compulsory,trim=TRUE),'_N_',N)
 if(is.null(data))data<-readRDS(file.path('08_boredom_results',paste0(if(variant=='original')'examples_'else 'more_boredom_',tag,'.rds')))
 B<-data$B;a<-data$alpha;alim<-range(a);hdr<-boredom_hdr()
 if(contours){ur<-range(unlist(lapply(data$prep,function(p)p$utility)));ulevels<-pretty(ur,12);vir<-hcl.colors(201,'Viridis');ucols<-vir[round((seq(0,1,length.out=201)^.65)*200)+1]}
 draw<-function(){
  # Identical panel geometry in the allocation and utility figures. Reserve
  # the full legend strip even when contours make the colour legend redundant.
  widths<-c(1,1,1,1,.52)
  layout(matrix(c(1:4,9,5:8,9),2,byrow=TRUE),widths=widths,heights=c(1,1))
  eta_at<-(4+.52/2)/sum(widths);par(oma=c(2,3.5,4.5,.2))
  for(i in 1:2)for(j in 1:4){
   p<-data$prep[[j]];q<-p$maps[[i]];cols<-boredom_example_colours(q,B,compulsory)
   par(mar=c(2.75,2.75,.5,.2),mgp=c(1.75,.48,0),tcl=-.23,pty='s',xpd=FALSE)
   plot(NA,xlim=alim,ylim=alim,xaxs='i',yaxs='i',axes=FALSE,xlab='',ylab='',asp=1)
   cols[q>=B-1e-5]<-'#B54343' # maxxer fill comes from a rooted polygon
   if(contours){u<-p$utility[[i]];cols<-matrix(ucols[pmax(1,pmin(201,1+floor(200*(u-ur[1])/diff(ur))))],nrow(u))}
   rasterImage(as.raster(t(cols)[ncol(cols):1,,drop=FALSE]),min(a),min(a),max(a),max(a),interpolate=TRUE)
   cl<-p$boundaries[[i]];lower<-cl[[2]]$x;upper<-cl[[3]]$x
   if(!contours){
    edge<-cl[[4]]$x;edge[!is.finite(edge)]<-max(a);edge[q[1,]>=B-1e-5]<-min(a)
    polygon(c(edge,max(a),max(a)),c(a,max(a),min(a)),col='#E06',border=NA)
    edge0<-cl[[1]]$x;edge0[!is.finite(edge0)]<-min(a);edge0[q[nrow(q),]<=1e-5]<-max(a)
    polygon(c(min(a),edge0,min(a)),c(min(a),a,max(a)),col=boredom_example_colours(matrix(0,1,1),B,compulsory)[1],border=NA)
   }
   upper[!is.finite(upper)&abs(q[nrow(q),]-B/2)<1e-4]<-max(a)
   valid<-which(is.finite(lower)&is.finite(upper)&upper>=lower)
   if(!contours&&length(valid)>1)for(ii in split(valid,cumsum(c(1,diff(valid)!=1))))if(length(ii)>1)
    polygon(c(lower[ii],rev(upper[ii])),c(a[ii],rev(a[ii])),col='#FF0',border=NA)
   abline(0,1,col=if(contours)'#888888'else '#000',lty=2,lwd=.9)
   if(contours)contour(a,a,u,levels=ulevels,add=TRUE,drawlabels=TRUE,labcex=.85,lwd=.7,col='#444444')
   for(h in hdr)lines(h[1,],h[2,],col='#EEE',lwd=.9)
   lines(cl[[1]]$x,cl[[1]]$y,col='black',lwd=1.2,lty=2)
   lines(cl[[4]]$x,cl[[4]]$y,col='#E06',lwd=1.4)
   if(abs(compulsory-.5)>1e-8)lines(cl[[2]]$x,cl[[2]]$y,col='#FF0',lwd=3)
   box()
   m<-if(i==1)p$first$mapped else p$responsive$mapped;dy<-diff(alim)
   x0<-max(a)+.015*dy;x1<-max(a)+.075*dy;yy<-min(a)+c(0,m[4],m[4]+m[1],B)*dy/B
   for(k in 1:3)rect(x0,yy[k],x1,yy[k+1],col=c('black','#E06','#08F')[k],border=NA,xpd=NA)
   rect(x0,min(a),x1,max(a),border='black',lwd=.7,xpd=NA)
   labs<-list(expression(bar(d)),expression(bar(z)[1]),expression(bar(z)[2]))
   for(k in 1:3)if(yy[k+1]-yy[k]>.055*dy)text((x0+x1)/2,mean(yy[k:(k+1)]),labs[[k]],cex=.67,xpd=NA,col='white',font=2)
   if(i==2){axis(1,at=seq(0,2.5,.5),cex.axis=1.2);mtext(expression(alpha[1]~~'(talent in tournament 1)'),1,line=2.4,cex=.95)}
   if(j==1){axis(2,at=seq(0,2.5,.5),las=1,cex.axis=1.2);mtext(expression(alpha[2]~~'(talent in tournament 2)'),2,line=2.4,cex=.95)}
   if(i==1){
    if(paper)mtext(LETTERS[j],3,line=1.8,at=min(a),adj=0,cex=1.25,font=2,xpd=NA)
    mtext(bquote(N==.(p$N)),3,line=.15,at=min(a),adj=0,cex=1.0,xpd=NA)
   }
   if(i==1)for(t in 1:2){
    at<-min(a)+c(.48,.88)[t]*diff(alim)
    mtext(bquote(beta[.(t)]==.(sprintf('%.2f',p$beta[t]))),3,line=1.55,at=at,cex=1.12,col=c('#E06','#08F')[t])
    mtext(bquote(varphi[.(t)]==.(sprintf('%.2f',p$phi[t]))),3,line=.15,at=at,cex=1.12,col=c('#E06','#08F')[t])
   }
  }
  mtext(expression('Fixed '*E[t]),2,outer=TRUE,at=.76,line=1.25,cex=1.25)
  mtext(expression('Responsive '*E[t]),2,outer=TRUE,at=.26,line=1.25,cex=1.25)
  mtext(bquote(c==.(sprintf('%.2f',compulsory))),3,outer=TRUE,at=eta_at,line=1.0,cex=1.12)
  if(contours){par(mar=c(4.8,4.2,.95,1.6),pty='m');plot.new()}else{
  par(mar=c(4.8,4.2,.95,1.6),pty='m',xpd=FALSE)
  nr<-ns<-181;d<-seq(0,1,length.out=ns);r<-seq(0,1,length.out=nr)
  gr<-expand.grid(d=d,r=r);vv<-colorRamp(c('#315A9B','white','#B54343'))(gr$r)/255*(1-gr$d)+(85/255)*gr$d
  cm<-matrix(rgb(vv[,1],vv[,2],vv[,3]),ns,nr)
  plot.new();plot.window(xlim=c(0,1.1),ylim=c(0,1),xaxs='i',yaxs='i')
  rasterImage(as.raster(t(cm)[nr:1,,drop=FALSE]),0,0,1,1,interpolate=FALSE);box()
  axis(1,at=c(0,.5,1),labels=c('','0.5','1.0'),cex.axis=.9,tcl=-.2)
  axis(2,at=seq(0,1,.1),labels=sprintf('%.1f',seq(0,1,.1)),las=1,cex.axis=.9,tcl=-.15)
  rect(0,.982,.15,1,col='#E06',lwd=.8,xpd=NA);rect(0,.491,.15,.509,col='#FF0',lwd=.8,xpd=NA)
  rect(0,0,.15,.018,col='#08F',lwd=.8,xpd=NA);rect(1,0,1.12,1,col='black',lwd=.8,xpd=NA)
  mtext(expression(d/B),1,line=2.2,cex=.9);mtext('(bored time',1,line=3.25,cex=.9);mtext('share)',1,line=4.3,cex=.9)
  mtext('Relative tournament-1 share among productive time',2,line=1.65,cex=.9)
  text(.15,.025,'2',adj=c(0,0),cex=1,col='white',font=2);text(.15,.975,'1',adj=c(0,1),cex=1,col='white',font=2)
  text(1.25,.5,'bored time',srt=90,cex=1,xpd=NA);text(-.09,1.032,'1maxxers',cex=1,xpd=NA,col='#E06')
  text(-.09,-.035,'2maxxers',cex=1,xpd=NA,col='#08F')
  }
 }
 if(is.null(stem))stem<-paste0(if(variant=='original')'Figure_local_boredom_examples_'else 'Figure_local_boredom_more_examples_',tag)
 if(png){grDevices::png(paste0(stem,'.png'),width=3640/300,height=6.6,units='in',res=600);tryCatch(draw(),finally=dev.off())}
 if(pdf){cairo_pdf(paste0(stem,'.pdf'),width=3640/300,height=6.6,family='Arial');tryCatch(draw(),finally=dev.off())}
 if(interactive())draw();invisible(data)
}
if(sys.nframe()==0L)plot_boredom_examples()
