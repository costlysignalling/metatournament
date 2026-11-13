# Parameter concordance guide: no numerical population calculations are needed.
# Run from the scripts folder to regenerate the PNG (600 dpi) and PDF.
draw_parameter_concordance <- function() {
  phi <- c(0.66, 1, 1.5)
  green <- "#65B880"; red <- "#E58683"; neutral <- "#EEEEEE"
  layout(rbind(matrix(1:9,nrow=3,byrow=TRUE),c(10,10,10)),
         heights=c(1,1,1,.30))
  par(oma=c(.1,2,1.5,.15),cex=1,family="sans")
  for (i in 1:3) for (j in 1:3) {
    par(mar=c(1.5,1.5,.35,.2),mgp=c(1,.2,0),tcl=-.14,
        xaxs="i",yaxs="i",pty="s",xpd=FALSE,cex=1)
    plot.new(); plot.window(c(0,1),c(0,1))
    if (i==j) rect(0,0,1,1,col=neutral,border=NA) else {
      polygon(c(0,1,1),c(0,0,1),col=if(j>i) green else red,border=NA)
      polygon(c(0,0,1),c(0,1,1),col=if(j>i) red else green,border=NA)
    }
    abline(0,1,lty=3,lwd=0.7); box(lwd=0.8)
    axis(1, at=c(0,1), cex.axis=.72, mgp=c(1,0,0), padj=-0.2, lwd=.8, lwd.ticks=.8)
        axis(2, at=c(0,1), las=1, cex.axis=.72, mgp=c(1,.2,0), lwd=.8, lwd.ticks=.8)
    mtext(expression(beta[1]),1,line=.3,cex=1.0,xpd=NA)
    mtext(expression(beta[2]),2,line=.3,cex=1.0,las=1,xpd=NA)
    if (i==1) mtext(bquote(varphi[1]==.(formatC(phi[j],format="f",digits=2))),
                    3,line=.45,cex=1.10,xpd=NA)
    if (j==1) mtext(bquote(varphi[2]==.(formatC(phi[i],format="f",digits=2))),
                    2,line=1.85,cex=1.10,xpd=NA)
  }
  par(mar=rep(0,4),pty="m",xpd=NA,cex=1)
  plot.new(); plot.window(c(0,1),c(0,1),xaxs="i",yaxs="i")
  legend(.5,.45,legend=expression(Concordant,Discordant,"Equal "*varphi),
         fill=c(green,red,neutral),border=NA,horiz=TRUE,xjust=0.5,yjust=0.5,
         bty="n",cex=.90,x.intersp=.5,seg.len=.9)
}
png("Figure_parameter_concordance_guide.png",width=4.0,height=4,units="in",res=600,pointsize=12)
draw_parameter_concordance(); dev.off()
pdf("Figure_parameter_concordance_guide.pdf",width=4.0,height=4,pointsize=12,useDingbats=FALSE)
draw_parameter_concordance(); dev.off()
