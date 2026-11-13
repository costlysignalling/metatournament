#Run this after setting up functions in 04_known_talents_strategy_map.R
pilot_map<-function(beta,phi,alim=c(0,2.5)){
  alpha<-seq(alim[1],alim[2],length.out=nalpha)
  sm<-strategy_map(alpha,B,beta,phi,nx,nroot,noverlap)
  z<-sm$z
  bound<-sm$bound
  
  par(mar=c(4,4,3,1),mgp=c(2.4,0.7,0),pty="s")
  
  image(alpha,alpha,z,col=cols,zlim=c(0,1),
        xlim=alim,ylim=alim,xaxs="i",yaxs="i",
        xlab=expression(alpha[1]~~"(talent in tournament 1)"),
        ylab=expression(alpha[2]~~"(talent in tournament 2)"),
        useRaster=TRUE)
  
  polygon(c(bound,alim[1],alim[1]),
          c(alpha,max(alpha),alim[1]),
          col=col2_max,border=NA)
  
  polygon(c(alpha,max(alpha),alim[1]),
          c(bound,alim[1],alim[1]),
          col=col1_max,border=NA)
  
  splitidx<-which(sm$split_global)
  if(length(splitidx)>1){
    groups<-split(splitidx,cumsum(c(1,diff(splitidx)!=1)))
    for(r in groups){
      if(length(r)>1){
        lines(alpha[r],alpha[r],col=coldiag,lwd=3)
      }
    }
  }
  
  box()
  title(bquote(beta==.(beta)~~","~~varphi==.(phi)))
}

pilot_map(0.8,1)

abline(a=0,b=(1+50)^(0.8-1),col=3,lwd=1)

(1+50)^(0.8-1)

al1<-seq(0,3,0.1)
al2<-al1*(1+50)^(0.8-1)
points(al1,al2)

al2<-seq(0,3,0.1)
al1<-al2*(1/((1+50)^(0.8-1)))
lines(al1,al2)


al2<-seq(0,3,0.1)
al1<-al2*(1+50)^(0.8-1)
lines(al1,al2,col=2,lwd=2)


al1<-al2*(1+50)^(0.4-1)

1/((1+50)^(0.4-1))

1/((1+50)^(0.8-1))
