# Fractions entering k of T identical tournaments, T=2,...,10.
# Data are calculated in C++ and called directly from R through Rcpp; Python
# and reticulate are not required. The C++ source, R calculation/validation
# scripts, raw integration results and plotting summaries are all supplied in
# Etude3_multiple_tournaments_analysis (see REPRODUCING_PARTICIPATION_FIGURES.md).
# Set the working directory to this scripts folder before running this file.
# By default this file only redraws the figure using the saved CSV summary.
# To recalculate EVERY probability from source, uncomment BOTH lines below.
# Requires Rcpp and a C++17 toolchain (Rtools on Windows); this can take hours.
# source("reproduce_participation_figures.R")
# reproduce_participation_figures(recompute=TRUE, figures="proportional")
# The helper also redraws the figure after computation. Recalculation replaces
# the saved raw data, summaries, convergence tables and PNG for this figure.
b<-25  #Total budget is B_T=b*T; T=2 therefore retains B=50.
betas<-c(.4,.8)
phis<-c(.66,1,1.25,1.5)
Ts<-2:10
data_file<-'participation_many_summary.csv'
if(!file.exists(data_file))data_file<-file.path('Etude3_multiple_tournaments_analysis','proportional_budget',data_file)
result<-read.csv(data_file)
cols<-hcl.colors(length(Ts),'Dark 3')
pchs<-rep(16,length(Ts))

draw_participation<-function(){
  par(mfrow=c(2,4),oma=c(3.0,3.8,3.5,.5),
      mar=c(2.7,2.2,1.5,.15),mgp=c(1.6,.45,0),tcl=-.2)
  for(i in seq_along(betas))for(j in seq_along(phis)){
    beta<-betas[i];phi<-phis[j]
    plot(NA,xlim=c(10.2,.8),ylim=c(-3,103),
         xaxs='i',yaxs='i',axes=FALSE,xlab='',ylab='')
    abline(h=seq(0,100,20),col='grey85',lwd=.65)
    abline(v=1:10,col='grey93',lwd=.55)
    axis(1,at=10:1,cex.axis=1.05)
    axis(2,at=seq(0,100,20),las=1,cex.axis=1.15)
    for(t in seq_along(Ts)){
      part<-subset(result,beta==betas[i]&phi==phis[j]&T==Ts[t])
      part<-part[order(part$k,decreasing=TRUE),]
      lines(part$k,100*part$proportion,col=cols[t],lwd=1.5,
            type='o',pch=pchs[t],cex=.85)
    }
    box()
    if(i==1) mtext(bquote(varphi==.(formatC(phi,format="f",digits=2))),
                    side=3,line=.6,cex=1.25)
  }
  mtext('Number of tournaments entered (k)',side=1,outer=TRUE,line=1.3,cex=1.05)
  mtext('Fraction of actors (%)',side=2,outer=TRUE,line=2.4,cex=1.05)
  mtext(bquote(beta==.(formatC(betas[1],format="f",digits=2))),side=2,outer=TRUE,line=.6,at=.75,cex=1.2)
  mtext(bquote(beta==.(formatC(betas[2],format="f",digits=2))),side=2,outer=TRUE,line=.6,at=.25,cex=1.2)
  par(fig=c(0,1,0,1),oma=c(0,0,0,0),mar=c(0,0,0,0),new=TRUE)
  plot.new()
  plot.window(c(0,1),c(0,1),xaxs='i',yaxs='i')
  legend(.5,.985,legend=paste0('T = ',Ts),col=cols,pch=pchs,lwd=1.5,
         horiz=TRUE,xjust=.5,yjust=1,bty='n',cex=1.15,x.intersp=.6)
}
png('Figure_known_talents_multiple_tournaments_proportional_budget.png',
    width=11.2,height=6.5,units='in',res=300)
draw_participation()
dev.off()
#if(interactive()) draw_participation()
