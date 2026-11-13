# Fixed-performance rectangles, continuously accumulated or interspersed.
# Each actor has exactly the same area alpha_i * x_i in the two panels.
plot_fixed_performance_harvesting<-function(B=50,alpha=c(.8,1.5,1.4),investment=c(.6,.3,.5)*B,png=TRUE,pdf=TRUE){
 stopifnot(length(alpha)==3,length(investment)==3,B>0,all(alpha>0),all(investment>=0),all(investment<=B))
 colors<-c('#D09','#0A0','#08F')
 slots<-local({had<-exists('.Random.seed',.GlobalEnv,inherits=FALSE);if(had)old<-get('.Random.seed',.GlobalEnv);on.exit(if(had)assign('.Random.seed',old,.GlobalEnv)else rm('.Random.seed',envir=.GlobalEnv));set.seed(2601);lapply(investment,function(v)sort(sample.int(100,round(100*v/B))))})
 draw<-function(){par(mfrow=c(1,2),oma=c(.1,.1,.1,.1),mar=c(4,4.4,3,1.5),mgp=c(2,.65,0),tcl=-.25,bty='l',xaxs='i',yaxs='i',cex=1.05)
  for(panel in 1:2){plot(NA,xlim=c(0,B),ylim=c(0,max(alpha)*1.3),xlab=if(panel==1)expression('Investment '*x)else 'Clock time',ylab=expression('Fixed performance '*P==alpha))
   for(i in 1:3){
    if(panel==1){rect(0,0,investment[i],alpha[i],col=adjustcolor(colors[i],alpha.f=.12),border=NA);segments(0,alpha[i],investment[i],alpha[i],col=colors[i],lwd=2.2)}else{
     for(k in slots[[i]]){rect((k-1)*B/100,0,k*B/100,alpha[i],col=adjustcolor(colors[i],alpha.f=.12),border=NA);segments((k-1)*B/100,alpha[i],k*B/100,alpha[i],col=colors[i],lwd=2.2)}
    }
    abline(h=alpha[i],col=adjustcolor(colors[i],.6),lty=2,lwd=.7)
   }
   mtext(LETTERS[panel],3,line=1.5,at=0,adj=0,cex=1.35,font=2)
   mtext(c('Accumulated investment','Interchangeable time')[panel],3,line=1.5,cex=1.12)
  }
 }
 stem<-'Figure_fixed_performance_harvesting'
 if(png){grDevices::png(paste0(stem,'.png'),width=9,height=3.5,units='in',res=600);tryCatch(draw(),finally=dev.off())}
 if(pdf){cairo_pdf(paste0(stem,'.pdf'),width=9,height=3.5,family='Arial');tryCatch(draw(),finally=dev.off())}
 if(interactive())draw();invisible(slots)
}
if(sys.nframe()==0L)plot_fixed_performance_harvesting()
