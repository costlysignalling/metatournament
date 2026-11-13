# Bored time alone, using the same computed populations as the invested-time map.
# A common scale spans the calculated population-mean range over both conditions.
# Redraw from saved results; no model recompilation is required.
plot_boredom_level<-function(compulsory=.5,N=30,png=TRUE,pdf=TRUE){
 tag<-paste0('c_',format(compulsory,trim=TRUE),'_N_',N)
 refined_file<-file.path('08_boredom_results',paste0('level_refined_',tag,'.rds'))
 if(!file.exists(refined_file))stop('Run refine_boredom_level() in 08_refine_boredom_level.R before redrawing this figure.')
 d<-readRDS(refined_file)
 stopifnot(all(vapply(d$panels,function(p)all(p$status==1),TRUE)))
 # s0 is mean bored time / B, not the share of actors who are bored.
 fields<-lapply(d$panels,function(p)list(p$first$s0,p$responsive$s0))
 lim<-range(unlist(fields));stopifnot(all(is.finite(lim)),diff(lim)>0)
 message(sprintf('Calculated mean bored-time range: %.8f%% to %.8f%% of B',100*lim[1],100*lim[2]))
 xs<-seq(min(d$beta),max(d$beta),length.out=601)
 interp<-function(z){
  v<-vapply(seq_along(d$beta),function(j)splinefun(d$beta,z[,j],method='monoH.FC')(xs),numeric(length(xs)))
  out<-t(vapply(seq_along(xs),function(i)splinefun(d$beta,v[i,],method='monoH.FC')(xs),numeric(length(xs))))
  out[]<-pmax(lim[1],pmin(lim[2],out));out
 }
 cols<-gray(seq(.995,.045,length.out=501))
 draw<-function(){
  layout(rbind(c(1,2,3,0,10,11,12),c(4,5,6,0,13,14,15),c(7,8,9,0,16,17,18),rep(19,7)),
   widths=c(1,1,1,.32,1,1,1),heights=c(1,1,1,.80))
  par(oma=c(.45,3.35,6.55,.45))
  for(mode in 1:2)for(r in 1:3)for(c in 1:3){
   p<-d$panels[[c,r]];z<-interp(if(mode==1)p$first$s0 else p$responsive$s0)
   par(mar=c(1.5,1.5,.5,.5),mgp=c(1.35,.36,0),tcl=-.24,xaxs='i',yaxs='i',pty='s',xpd=FALSE)
   image(xs,xs,z,col=cols,zlim=lim,axes=FALSE,xlab='',ylab='',useRaster=TRUE)
   abline(0,1,lty=3,lwd=.7,col='#888888')
   if(r==3){axis(1,at=c(.2,.4,.6,.8),cex.axis=1.02);mtext(expression(beta[1]),1,line=1.9,cex=1.062)}
   if(c==1){axis(2,at=c(.2,.4,.6,.8),cex.axis=1.02);mtext(expression(beta[2]),2,line=1.4,cex=1.062,las=1)
    text(par('usr')[1]-.36,mean(par('usr')[3:4]),bquote(varphi[2]==.(sprintf('%.2f',p$phi2))),srt=90,xpd=NA,cex=1.77)}
   if(r==1)mtext(bquote(varphi[1]==.(sprintf('%.2f',p$phi1))),3,line=.22,cex=1.18)
   box(lwd=.8)
  }
  par(mar=c(.1,.1,.1,.1),pty='m',xpd=NA);plot.new();plot.window(xlim=c(0,1),ylim=c(0,1),xaxs='i',yaxs='i')
  rasterImage(as.raster(matrix(cols,nrow=1)),.24,.30,.76,.55);rect(.24,.30,.76,.55,lwd=.8)
  at<-seq(lim[1],lim[2],length.out=5);xx<-.24+.52*(at-lim[1])/diff(lim)
  segments(xx,.30,xx,.24);text(xx,.14,sprintf('%.3f',100*at),cex=1.1)
  text(.5,.75,'Population-average bored time (% of total budget)',cex=1.45)
  mtext(expression('Fixed '*E[t]),3,outer=TRUE,at=.225,line=3.95,cex=1.55)
  mtext('Bored time after the first talent-informed response',3,outer=TRUE,at=.225,line=2.25,cex=1.1)
  mtext(expression('Responsive '*E[t]),3,outer=TRUE,at=.775,line=3.95,cex=1.55)
  mtext('Bored time at the responsive density',3,outer=TRUE,at=.775,line=2.25,cex=1.1)
  mtext(bquote(c==.(sprintf('%.2f',compulsory))~~N==.(N)),3,outer=TRUE,at=.5,line=4.3,cex=1.55)
 }
 stem<-paste0('Figure_local_boredom_level_',tag)
 if(png){grDevices::png(paste0(stem,'.png'),width=3500/300,height=7.6,units='in',res=600);tryCatch(draw(),finally=dev.off())}
 if(pdf){cairo_pdf(paste0(stem,'.pdf'),width=3500/300,height=7.6,family='Arial');tryCatch(draw(),finally=dev.off())}
 if(interactive())draw()
 invisible(list(range=lim,data=d))
}
if(sys.nframe()==0L)plot_boredom_level()
