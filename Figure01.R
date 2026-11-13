w<-2.6

widths<-c(3,0.8,3,0.8,0.4,1.2-0.4)
heights<-c(4.15,1)

png("Figure01.png",
    w=(w/widths[1])*sum(widths),
    h=(w/heights[1])*sum(heights),
    units="in",res=600)

cols<-c("#D09","#0C0","#08D")
col1_div<-"#B54343"
col2_div<-"#315A9B"

layout(matrix(c(1,2,3,4,5,5,
                6,0,7,8,8,0),nrow=2,byrow=TRUE),
       widths=widths,
       heights=heights)

B<-50
l<-401

beta1<-0.6
beta2<-0.4

r<-c(0.6,0.3,0.5)
tals1<-c(1.1,0.9,1.7)
tals2<-c(0.8,1.5,1.4)

phi1<-0.5
phi2<-2

x1<-r*B
x2<-(1-r)*B

P1<-tals1*((1+x1)^beta1-1)
P2<-tals2*((1+x2)^beta2-1)

W1<-(1+P1)^phi1-1
W2<-(1+P2)^phi2-1

rs1<-W1/sum(W1)
rs2<-W2/sum(W2)

E1<-sum(x1)
E2<-sum(x2)

U1<-E1*rs1
U2<-E2*rs2
U<-U1+U2

alpha_lab1<-as.expression(lapply(1:3,function(i) bquote(alpha[.(i)]==.(tals1[i]))))
alpha_lab2<-as.expression(lapply(1:3,function(i) bquote(alpha[.(i)]==.(tals2[i]))))

plot_energy_strip<-function(xv,xlab_expr,ylab_expr){
  par(mar=c(2.3,3.5,0,1),mgp=c(1.4,0.38,0))
  
  mids<-barplot(rev(xv),
                horiz=TRUE,
                col=adjustcolor(rev(cols),alpha.f=0.35),
                border=rev(cols),
                axes=FALSE,
                xlim=c(0,B),
                space=0.35,
                names.arg=rep("",3),
                xlab=xlab_expr)
  
  axis(1,at=seq(0,B,10),labels=seq(0,B,10),tck=-0.04)
  box(bty="l")
  
  text(x=-0.04*B,y=mids,labels=3:1,xpd=NA,adj=1,
       col=rev(cols),cex=0.9)
  
  mtext(ylab_expr,2,line=1.9,cex=0.7)
}

eps<-4
halfw<-0.35

# panel 1
par(mar=c(3.5,3.5,2,1),mgp=c(2,0.65,0))
for(i in 1:3){
  x<-seq(0,x1[i],l=l)
  y<-tals1[i]*((1+x)^beta1-1)
  if(i==1){
    plot(x,y,xaxs="i",yaxs="i",type="l",lwd=2,col=cols[i],
         ylim=c(0,12),xlim=c(0,B),bty="n",
         ylab=expression(alpha((1+x)^beta-1)),xlab=expression(x[1]))
  } else {
    lines(x,y,lwd=2,col=cols[i])
  }
  abline(h=P1[i],lty=2,col=cols[i])
  for(i in 1:3){
    text(x=50,y=P1[i]+0.35,
         labels=bquote(P[.(i)]),
         col=cols[i],cex=0.9,pos=2)
  }
}
title(main="Tournament 1",adj=0,col.main=col1_div,cex.main=1.4)
mtext(bquote(beta==.(beta1)),3,line=0.2,adj=1,col=col1_div,cex=0.85)
mtext(bquote(varphi==.(phi1)),3,line=0.2,at=64,
      col=col1_div,xpd=NA,cex=0.85)
legend("topleft",legend=alpha_lab1,text.col=cols,bty="n",cex=1,inset=c(-0.05,0.035))

# panel 2
par(mar=c(3.5,3.5,2,0.8),mgp=c(2,0.65,0))
plot(NULL,xaxs="i",yaxs="i",ylim=c(0,12/sum(P1)),xlim=c(0,1),bty="n",
     axes=TRUE,xaxt="n",
     ylab=expression("Relative tournament spoils" %prop% (1+P)^varphi-1),
     xlab="")
abline(h=P1/sum(P1),col=cols,lty=2)
abline(h=rs1,col=cols,lty=1,lwd=2)

# panel 3
par(mar=c(3.5,3.5,2,1),mgp=c(2,0.65,0))
for(i in 1:3){
  x<-seq(0,x2[i],l=l)
  y<-tals2[i]*((1+x)^beta2-1)
  if(i==1){
    plot(x,y,xaxs="i",yaxs="i",type="l",lwd=2,col=cols[i],
         ylim=c(0,12),xlim=c(0,B),bty="n",
         ylab=expression(alpha((1+x)^beta-1)),xlab=expression(x[2]))
  } else {
    lines(x,y,lwd=2,col=cols[i])
  }
  abline(h=P2[i],lty=2,col=cols[i])
  for(i in 1:3){
    text(x=50,y=P2[i]+0.35,
         labels=bquote(P[.(i)]),
         col=cols[i],cex=0.9,pos=2)
  }
}
title(main="Tournament 2",adj=0,col.main=col2_div,cex.main=1.4)
mtext(bquote(beta==.(beta2)),3,line=0.2,adj=1,col=col2_div,cex=0.85)
mtext(bquote(varphi==.(formatC(phi2,format="f",digits=1))),3,line=0.2,at=64,
      col=col2_div,xpd=NA,cex=0.85)
legend("topleft",legend=alpha_lab2,text.col=cols,bty="n",cex=1,inset=c(-0.05,0.035))

# panel 4
par(mar=c(3.5,3.5,2,0.8),mgp=c(2,0.65,0))
plot(NULL,xaxs="i",yaxs="i",ylim=c(0,12/sum(P2)),xlim=c(0,1),bty="n",
     axes=TRUE,xaxt="n",
     ylab=expression("Relative tournament spoils" %prop% (1+P)^varphi-1),
     xlab="")
abline(h=P2/sum(P2),col=cols,lty=2)
abline(h=rs2,col=cols,lty=1,lwd=2)

# panel 5
par(mar=c(3.5,3,2,1),mgp=c(2,0.65,0))
Ut_mat<-rbind(U1,U2)

bp<-barplot(Ut_mat,col=NA,border=NA,
            names.arg=rep("",3),
            ylab="Total utility",
            xlab="Actor",
            ylim=c(0,1.05*max(U)))

for(j in 1:3){
  cumulative<-c(0,cumsum(Ut_mat[,j]))
  for(i in 1:2){
    y0<-cumulative[i]
    y1<-cumulative[i+1]
    rect(bp[j]-halfw,y0,bp[j]+halfw,y1,
         col=c(col1_div,col2_div)[i],border=cols[j],lwd=1.5)
  }
}

text(bp,par("usr")[3]-0.04*diff(par("usr")[3:4]),
     labels=1:3,xpd=NA,col=cols)

# panel 6
plot_energy_strip(x1,expression(x[1]),bquote(E[1]==.(E1)))

# panel 7
plot_energy_strip(x2,expression(x[2]),bquote(E[2]==.(E2)))

# panel 8
par(mar=c(2.2,3.8,0,0),mgp=c(1.6,0.65,0))
Et_mat<-cbind(x1,x2)

bp<-barplot(Et_mat,col=NA,border=NA,
            names.arg=c("T1","T2"),las=1,
            ylim=c(0,1.08*max(E1,E2)),ylab="")

for(j in 1:2){
  cumulative<-c(0,cumsum(Et_mat[,j]))
  for(i in 1:3){
    y0<-cumulative[i]
    y1<-cumulative[i+1]
    if(i>1) y0<-y0+eps/2
    if(i<3) y1<-y1-eps/2
    rect(bp[j]-halfw,y0,bp[j]+halfw,y1,
         col=c(col1_div,col2_div)[j],border=cols[i],lwd=1.5,xpd=T)
  }
}

mtext(expression(E[t]),2,line=2.2,las=1,cex=0.7,xpd=NA)

dev.off()
