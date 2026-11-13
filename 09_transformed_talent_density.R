# Exact change-of-variable density for the equivalent harvesting talent.
# Run from the manuscript scripts folder; reuses the CDF/quantile helpers.
source('09_harvesting_explanations.R')
d_harvesting_talent<-function(z,phi,sdlog=harvesting_slog(),log=FALSE){
 stopifnot(length(phi)==1L,is.finite(phi),phi>0,sdlog>0)
 out<-rep(-Inf,length(z));ok<-is.finite(z)&z>0
 alpha<-expm1(log1p(z[ok])/phi)
 out[ok]<-dlnorm(alpha,meanlog=0,sdlog=sdlog,log=TRUE)-base::log(phi)+(1/phi-1)*log1p(z[ok])
 out[is.na(z)]<-NA_real_
 if(log)out else exp(out)
}
plot_transformed_talent_density<-function(phis=c(1.25,2),png=TRUE,pdf=TRUE){
 stopifnot(length(phis)>=1,all(is.finite(phis)),all(phis>0))
 slog<-harvesting_slog();blue<-'#08F'
 limits<-vapply(phis,function(p)max(qlnorm(.995,0,slog),q_harvesting_talent(.995,p)),0.)
 fields<-Map(function(phi,limit){
  x<-sort(unique(c(0,seq(0,limit,length.out=2501),exp(seq(log(1e-5),log(limit),length.out=2501)))))
  list(x=x,original=dlnorm(x,0,slog),transformed=d_harvesting_talent(x,phi))
 },phis,limits)
 ymax<-1.08*max(vapply(fields,function(p)max(p$original,p$transformed),0.))
 draw<-function(){
  par(mfrow=c(1,length(phis)),mar=c(3.3,3.7,1.8,.6),mgp=c(2,.65,0),
      tcl=-.25,bty='l',xaxs='i',yaxs='i',cex=1)
  for(j in seq_along(phis)){
   d<-fields[[j]]
   plot(NA,xlim=c(0,limits[j]),ylim=c(0,ymax),xlab='Talent value',
        ylab=if(j==1)'Probability density'else '',las=1)
   polygon(c(d$x,rev(d$x)),c(d$transformed,rep(0,length(d$x))),col=adjustcolor(blue,.08),border=NA)
   lines(d$x,d$original,col='black',lwd=2)
   lines(d$x,d$transformed,col=blue,lwd=2.2)
   mtext(bquote(varphi==.(sprintf('%.2f',phis[j]))),3,line=.4,cex=1.15)
   if(j==length(phis))legend('topright',legend=c(expression('Original '*alpha),expression('Transformed '*Z)),
                 col=c('black',blue),lwd=c(2,2.2),bty='n',cex=.95,x.intersp=.65,y.intersp=.9)
  }
 }
 stem<-'Figure_harvesting_transformed_talent_density';width<-3.75*length(phis)
 if(png){grDevices::png(paste0(stem,'.png'),width=width,height=3.2,units='in',res=600,
  type=if(capabilities('cairo'))'cairo'else getOption('bitmapType'));tryCatch(draw(),finally=dev.off())}
 if(pdf){cairo_pdf(paste0(stem,'.pdf'),width=width,height=3.2,family='Arial');tryCatch(draw(),finally=dev.off())}
 if(interactive())draw();invisible(fields)
}
if(sys.nframe()==0L)plot_transformed_talent_density()
