
g<-function(x,beta){
  (1+x)^beta-1
}

B<-10
x<-seq(1,B,0.01)
beta<-0.6
beta1<-0.6

w<-5

cols<-c("#D09","#088","#06B")
colg<-"#555"

beta1<-0.6
beta2<-0.4
beta3<-0.8

w<-5
png("FigureS01.png",w=w,h=w,units="in",res=600)

par(pty="s",mar=c(3.5,3.5,2,2),mgp=c(2,0.65,0))

x2<-seq(0,5,0.01)

plot(x2,x2^beta,xaxs="i",yaxs="i",type="l",lwd=2,col=cols[1],asp=1,
     ylim=c(0,5),xlim=c(0,5),bty="n",
     ylab=expression(x^beta),xlab=expression(x))

lines(x2,x2^beta2,lwd=2,col=cols[2])
lines(x2,x2^beta3,lwd=2,col=cols[3])

# curve labels
add_beta_label <- function(x, beta, col, offset=0.12){
  
  y <- x^beta
  
  # local slope of x^beta
  slope <- beta * x^(beta - 1)
  angle <- atan(slope) * 180 / pi
  
  text(x, y + offset,
       labels=bquote(beta == .(beta)),
       col=col,
       srt=angle,
       adj=c(0,0),
       xpd=NA)
}

add_beta_label(4.15, beta3, cols[3])  # beta = 0.8
add_beta_label(4.25, beta1, cols[1])  # beta = 0.6
add_beta_label(4.35, beta2, cols[2])  # beta = 0.4

abline(v=1,col="#CCC",lty=2)
abline(h=1,col="#CCC",lty=2)

segments(1,1,1,5,col=colg)
segments(1,1,5,1,col=colg)

ticks<-function(x,y,l=0.15,labels=NULL){
  if(length(x)>1 & length(y)==1){
    segments(x,y,x,y-l,xpd=NA,col=colg)
    text(x,y-l,labels,pos=1,xpd=NA,col=colg)
  }
  if(length(x)==1 & length(y)>1){
    segments(x-l,y,x,y,xpd=NA,col=colg)  
    text(x-l,y,labels,pos=2,xpd=NA,col=colg)
  }
}

ticks(1:5,1,labels=0:4)
ticks(1,1:5,labels=0:4)

mtext(expression((1+x)^beta-1),side=2,line=-2,adj=0.6,col=colg)
mtext(expression(x),side=1,line=-2,adj=0.6,col=colg)

dev.off()
