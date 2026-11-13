# Run from the manuscript scripts folder. This compact figure follows
# panel A of Figure_fixed_performance_harvesting, now showing spoils at phi=1.25.
plot_harvesting_paper_explanation<-function(phi=1.25,B=50,alpha=c(.8,1.5,1.4),
 investment=c(.6,.3,.5)*B,png=TRUE,pdf=TRUE){
 stopifnot(phi>0,B>0,length(alpha)==3,length(investment)==3,all(alpha>0),
           all(investment>=0),all(investment<=B))
 color<-c('#D09','#0A0','#08F');height<-expm1(phi*log1p(alpha))
 draw<-function(){
  par(mar=c(3,3.4,.2,.45),mgp=c(1.8,.5,0),tcl=-.25,bty='l',
      xaxs='i',yaxs='i',cex=1.05,pty='m')
  plot(NA,xlim=c(0,B),ylim=c(0,max(height)*1.3),xlab=expression('Investment '*x),
       ylab=bquote('Spoils weight '~(1+alpha)^.(phi)-1))
  for(i in seq_along(alpha)){
   rect(0,0,investment[i],height[i],col=adjustcolor(color[i],alpha.f=.12),border=NA)
   segments(0,height[i],investment[i],height[i],col=color[i],lwd=2.2)
   abline(h=height[i],col=adjustcolor(color[i],.6),lty=2,lwd=.7)
  }
 }
 stem<-'Figure_fixed_performance_harvesting_paper'
 if(png){grDevices::png(paste0(stem,'.png'),width=4.2,height=2.9,units='in',res=600,
  type=if(capabilities('cairo'))'cairo'else getOption('bitmapType'));tryCatch(draw(),finally=dev.off())}
 if(pdf){cairo_pdf(paste0(stem,'.pdf'),width=4.2,height=2.9,family='Arial');tryCatch(draw(),finally=dev.off())}
 if(interactive())draw()
 invisible(list(talent=alpha,spoils_rate=height,claim=investment*height))
}

# Exact transformed-talent CDF and quantile. Original talent has median 1, SD .5.
harvesting_slog<-function()sqrt(log((1+sqrt(2))/2))
q_harvesting_talent<-function(p,phi,sdlog=harvesting_slog()){
 stopifnot(phi>0,sdlog>0,all(p>=0&p<=1))
 expm1(phi*log1p(qlnorm(p,meanlog=0,sdlog=sdlog)))
}
p_harvesting_talent<-function(z,phi,sdlog=harvesting_slog()){
 stopifnot(phi>0,sdlog>0)
 plnorm(expm1(log1p(pmax(z,0))/phi),meanlog=0,sdlog=sdlog)
}
if(sys.nframe()==0L)plot_harvesting_paper_explanation()
