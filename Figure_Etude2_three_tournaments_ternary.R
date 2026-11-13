# Etude 2: expected utility for three identical tournaments.
# Population allocations are uniform on the simplex (Dirichlet(1,1,1)),
# independently of talent. Each marginal share is Beta(1,2), not Uniform(0,1).
# A focal actor's shares can vary while the population normalization stays fixed.
stopifnot(requireNamespace('viridisLite',quietly=TRUE))
B <- 50; beta <- 0.4; phi <- 2.34
m <- 0
s <- sqrt(log((1+sqrt(2))/2)) # median 1; original-scale SD 0.5

gh <- function(n) {
  J <- matrix(0,n,n)
  J[cbind(seq_len(n-1),2:n)] <- sqrt(seq_len(n-1))
  eig <- eigen(J+t(J),symmetric=TRUE)
  list(a=exp(m+s*eig$values),w=eig$vectors[1,]^2)
}
rule <- gh(80)
q <- function(x,quadrature=rule) vapply(x,function(xx)
  sum(quadrature$w*expm1(phi*log1p(quadrature$a*expm1(beta*log1p(xx))))),numeric(1))
mean_q <- function(x) q(x,rule)
M <- integrate(function(r) 2*(1-r)*mean_q(B*r),0,1,
               rel.tol=1e-11,subdivisions=1000)$value
utility <- function(r1,r2,r3) B/(3*M)*(mean_q(B*r1)+mean_q(B*r2)+mean_q(B*r3))
check_x <- seq(0,B,length.out=101)
stopifnot(max(abs(q(check_x,gh(160))-mean_q(check_x)))<1e-10*max(mean_q(check_x)))
stopifnot(abs(3*(B/(3*M))*M-B)<1e-10)

# Map r1 to bottom-left, r2 to bottom-right, r3 to top.
height <- sqrt(3)/2
n <- 240L
triangles <- vector('list',n*n)
index <- 0L
for(i in 0:(n-1)) for(j in 0:(n-1-i)) {
  index <- index+1L
  triangles[[index]] <- rbind(c(i,j),c(i+1,j),c(i,j+1))/n
  if(i+j<=n-2L) {
    index <- index+1L
    triangles[[index]] <- rbind(c(i+1,j),c(i+1,j+1),c(i,j+1))/n
  }
}
stopifnot(index==n*n)
centers <- t(vapply(triangles,colMeans,numeric(2)))
values <- utility(1-rowSums(centers),centers[,1],centers[,2])
edge_value <- utility(.5,.5,0)
center_value <- utility(1/3,1/3,1/3)
specialist_value <- utility(1,0,0)
limits <- range(c(values,edge_value,center_value,specialist_value))
palette <- viridisLite::viridis(256)
colors <- palette[pmax(1,pmin(256,1+floor(255*(values-limits[1])/diff(limits))))]
vertices <- lapply(triangles,function(t) cbind(t[,1]+t[,2]/2,height*t[,2]))
xx <- unlist(lapply(vertices,function(v) c(v[,1],NA_real_)))
yy <- unlist(lapply(vertices,function(v) c(v[,2],NA_real_)))

draw <- function() {
  par(mar=c(1,1,1,1),family='sans',xpd=NA)
  plot.new(); plot.window(xlim=c(-.13,1.29),ylim=c(-.12,1.02),asp=1)
  text(.5,1.035,expression(beta==0.4~' | '~varphi==2.34~' | '~B==50),cex=1.0)
  polygon(xx,yy,col=colors,border=NA)
  # Share grid: each family holds one share constant.
  for(r in seq(.2,.8,by=.2)) {
    segments(r,0,(1+r)/2,height*(1-r),col='#FFFFFF65',lwd=.65)
    segments(1-r,0,(1-r)/2,height*(1-r),col='#FFFFFF65',lwd=.65)
    segments(r/2,height*r,1-r/2,height*r,col='#FFFFFF65',lwd=.65)
    text(1-r,-.035,paste0(round(100*r),'%'),cex=.72,col='#555555')
    text(r/2-.045,height*r,paste0(round(100*r),'%'),cex=.72,col='#555555',srt=60)
    text((1+r)/2+.045,height*(1-r),paste0(round(100*r),'%'),cex=.72,col='#555555',srt=-60)
  }
  polygon(c(0,1,.5),c(0,0,height),border='#333333',lwd=1.2,col=NA)
  text(-.035,.02,'Tournament 1\n100%',adj=c(1,.5),cex=.88)
  text(1.035,.02,'Tournament 2\n100%',adj=c(0,.5),cex=.88)
  text(.5,height+.065,'Tournament 3: 100%',cex=.88)
  points(c(.5,.25,.75),c(0,height/2,height/2),pch=21,bg='white',col='#222222',cex=1.15,lwd=1.2)
  points(.5,height/3,pch=23,bg='white',col='#222222',cex=1.15,lwd=1.2)
  # Continuous viridis key, same scale as the triangle.
  keyx <- c(1.14,1.185); keyy <- c(.13,.74)
  breaks <- seq(keyy[1],keyy[2],length.out=257)
  rect(keyx[1],breaks[-257],keyx[2],breaks[-1],col=palette,border=NA)
  rect(keyx[1],keyy[1],keyx[2],keyy[2],border='#555555',col=NA)
  ticks <- pretty(limits,n=5); ticks<-ticks[ticks>=limits[1] & ticks<=limits[2]]
  for(v in ticks) {
    y<-keyy[1]+diff(keyy)*(v-limits[1])/diff(limits)
    segments(keyx[2],y,keyx[2]+.009,y,col='#555555')
    text(keyx[2]+.017,y,formatC(v,format='f',digits=2),adj=0,cex=.8)
  }
  text(mean(keyx),keyy[2]+.065,'Expected\nutility',cex=.9)
  legend(.05,-.145,legend=c('50:50 across two tournaments','Equal three-way split'),
         pch=c(21,23),pt.bg='white',col='#222222',horiz=TRUE,bty='n',cex=.79,xjust=0,yjust=0)
}
w<-6
png('Figure_Etude2_three_tournaments_ternary.png',width=w,height=(w/9)*6.5,units="in",res=600)
draw(); dev.off()
pdf('Figure_Etude2_three_tournaments_ternary.pdf',width=w,height=(w/9)*6.5,useDingbats=FALSE)
draw(); dev.off()
cat(sprintf('Uniform-population normalization M=%.12f\n',M))
cat(sprintf('Specialist: %.8f; equal three-way: %.8f; two-way split: %.8f\n',
            specialist_value,center_value,edge_value))
cat(sprintf('Color scale: %.8f to %.8f\n',limits[1],limits[2]))
