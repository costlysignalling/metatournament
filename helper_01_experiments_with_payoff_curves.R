#first diminishing returns function, pretty straighforward, power law for 0<b<1
skill1<-function(x,a,b){a*x^b}

b<-0.4
col=2
lwd=1.5
curve(skill1(x,a=1.0,b=b),from=0,100,ylim=c(0,15),n=10001,col=col,lwd=lwd)
curve(skill1(x,a=1.5,b=b),add=T,n=10001,col=col,lwd=lwd)
curve(skill1(x,a=0.5,b=b),add=T,n=10001,col=col,lwd=lwd)

b<-0.6
col=3
lwd=1.5
curve(skill1(x,a=1.0,b=b),add=T,n=10001,col=col,lwd=lwd)
curve(skill1(x,a=1.5,b=b),add=T,n=10001,col=col,lwd=lwd)
curve(skill1(x,a=0.5,b=b),add=T,n=10001,col=col,lwd=lwd)
#They also cross so why not to go with the one that does +1 in there

#It has a little awkward property S′(e)→∞ as e→0.
skill2<-function(x,a,b,S0){S0+a*((1+x)^b-1)}

#This just adds constant
curve(skill2(x,a=1,b=0.4,S0=0),from=0,100,ylim=c(0,15),n=10001)
curve(skill2(x,a=1,b=0.4,S0=1),n=10001,add=T)

#So I will do it without s0
skill<-function(x,a,b){a*((1+x)^b-1)}

curve(skill(x,a=1,b=0.4),from=0,100,ylim=c(0,15),n=10001)
curve(skill(x,a=1.5,b=0.4),add=T,n=10001)
curve(skill(x,a=0.5,b=0.4),add=T,n=10001)

#This is nice
b<-0.4
col=2
lwd=1.5
curve(skill(x,a=1.0,b=b),from=0,100,ylim=c(0,15),n=10001,col=col,lwd=lwd)
curve(skill(x,a=1.5,b=b),add=T,n=10001,col=col,lwd=lwd)
curve(skill(x,a=0.5,b=b),add=T,n=10001,col=col,lwd=lwd)

b<-0.6
col=3
lwd=1.5
curve(skill(x,a=1.0,b=b),add=T,n=10001,col=col,lwd=lwd)
curve(skill(x,a=1.5,b=b),add=T,n=10001,col=col,lwd=lwd)
curve(skill(x,a=0.5,b=b),add=T,n=10001,col=col,lwd=lwd)

#Since in the first experiment at least, we set the maximum investment time to 1, we need to check how these things look on that scale
b<-0.4
col=2
lwd=1.5
curve(skill(x,a=1.0,b=b),from=0,1,ylim=c(0,1),n=10001,col=col,lwd=lwd)
curve(skill(x,a=1.5,b=b),add=T,n=10001,col=col,lwd=lwd)
curve(skill(x,a=0.5,b=b),add=T,n=10001,col=col,lwd=lwd)

b<-0.6
col=3
lwd=1.5
curve(skill(x,a=1.0,b=b),add=T,n=10001,col=col,lwd=lwd)
curve(skill(x,a=1.5,b=b),add=T,n=10001,col=col,lwd=lwd)
curve(skill(x,a=0.5,b=b),add=T,n=10001,col=col,lwd=lwd)

#Ok, here ii does not have at first sight so nice diminishing returns profile (mathematically these things will be equivalent, nevertheless what will be important is the power coefficient and the total time, maybe we will eventually reduce it to a single parameter, but for now, just going with "50" years budget is fine)

#The function we use is just equivalent to y=x^b when we subtract 1 from both x and y. This makes it go nicely through the 0,0, point with finite derivations. This is what it is visually
b<-0.6
x<-seq(-2,4,l=1001)
y<-(1+x)^b-1

plot(x,y,type="l",asp=1)

x<-seq(0,5,l=1001)
y<-x^b

lines(x,y,col=4)
abline(h=0,lty=2)
abline(v=0,lty=2)
abline(h=-1,lty=2,col=2)
abline(v=-1,lty=2,col=2)

lines(x-1,y-1,lty=1,col=2)

#Alternatives
b<-0.75
x<-seq(-2,4,l=1001)
y<-(1+x)^b-1
lines(x,y)

x<-seq(0,5,l=1001)
y<-x^b
lines(x,y,col=4)
lines(x-1,y-1,lty=1,col=2)





#Is b inherent to the tournament? What it affects?

#And now I am looking at the utility function, which is just "the total mass of people I beat as recognized by others", whatever it means. I will figure it out

#Also figure out the correlations between tournaments



