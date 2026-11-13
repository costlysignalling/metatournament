#Two games two matching talents

library(viridis)

#In the first etude, I build a system with just two tournaments and a random talent
#for each. Players are randomly distributed in talent space. We compare the utility
#of specialists and generalists.

#Define tournaments, now both have the same training exponent
tours<-data.frame(beta=c(0.4,0.4))
ntour<-nrow(tours)

#I draw random participants with talents.
#I compared lognormal with gamma distributions and other alternatives earlier.
n<-4000
median<-1
sd<-0.5
m<-log(median)
s<-sqrt(log((1+sqrt(1+4*(sd/median)^2))/2))

#How good are rare people? 1:1000, 1:10^6, 1:world population
qlnorm(1-1/1000,m,s)
qlnorm(1-1/1000000,m,s)
qlnorm(1-1/(9*10^9),m,s)

#These qlnorm values are one-exceedance quantiles, not expected sample maxima:
#each has upper-tail probability 1/N, so N independent people give one expected
#exceedance. The chance that at least one exceeds the quantile is
#1-(1-1/N)^N, approaching 1-exp(-1), or about 63.2%.
#
#For comparison, the maximum A_N has cumulative probability
#Pr(A_N <= a) = pnorm((log(a)-m)/s)^N.
#If U is unit-exponential, exp(-U/N) has the distribution of the maximum of
#N independent uniform variables. Thus the expected maximum is the integral
#of exp(m+s*qnorm(exp(-u/N)))*exp(-u) over u>0.
#The function below sets u=exp(t) for numerical integration. Computing the
#upper-tail probability with -expm1(-u/N) avoids rounding exp(-u/N) to 1.
#A wider integration window can be used to check convergence.
expected_max_talent<-function(N,m=0,
    s=sqrt(log((1+sqrt(2))/2)),lower=-40,upper=5){
  stopifnot(length(N)==1,N>=1,s>0)
  integrand<-function(t){
    u<-exp(t)
    z<-qnorm(-expm1(-u/N),lower.tail=FALSE)
    exp(m+s*z+t-u)
  }
  integrate(integrand,lower,upper,rel.tol=1e-10,
            subdivisions=1000)$value
}
sapply(c(1000,1000000,9*10^9),expected_max_talent,m=m,s=s)
#Expected maxima: about 4.13054, 8.29593, and 16.33851 times the median.
#Both calculations support the manuscript's rounded illustrations of 4, 8, 16.


#Besides talent, randomly draw how much of their time people allocate to each
#tournament. I work with a "50 years total budget" to make diminishing returns clear.
budget<-50

rel1<-runif(n,0,1)
people<-list(tal=cbind(rlnorm(n,m,s),
                       rlnorm(n,m,s)),
             rel=cbind(rel1,
                       1-rel1),
             tour=cbind(rel1*budget,
                        (1-rel1)*budget))

str(people)

#Training first determines performance.
g<-function(x,beta){
  (1+x)^beta-1
}

performance<-function(x,alpha,beta){
  alpha*g(x,beta)
}

people$perf<-sapply(1:ntour,function(t){
  performance(people$tour[,t],people$tal[,t],tours$beta[t])
})
people$total<-rowSums(people$perf)

#I will plot a few players
nplot<-10
grain<-0.1

cols<-turbo(nplot)
lwd1<-2
lwd2<-1.8

layout(matrix(c(1,2,3),nrow=1),widths=c(3,3,1))
par(mar=c(3.5,3.5,1.5,0.5),mgp=c(2,0.65,0))

for(t in 1:ntour){
  plot(NULL,xlim=c(0,50),ylim=c(0,10),
       xlab="Time invested",ylab="Performance",
       xaxs="i",yaxs="i")
  title(paste("Tournament",t),adj=0)

  for(i in 1:nplot){
    times<-c(seq(0,people$tour[i,t],by=grain),people$tour[i,t])
    lines(times,
          performance(times,people$tal[i,t],tours$beta[t]),
          col=cols[i],lwd=lwd1)
  }

  abline(h=people$perf[1:nplot,t],lty=2,col=cols,lwd=lwd2)
}

plot(NULL,xlim=c(0,1),ylim=c(0,10),
     xlab="",ylab="Total performance",
     xaxt="n",xaxs="i",yaxs="i")
abline(h=people$total[1:nplot],lty=1,col=cols,lwd=lwd1)
title("Sum",adj=0)


#Tournament utility
#Performance is translated into a nonnegative spoils weight only after performance
#has been produced:
#
# W(P,phi) = (1+P)^phi - 1.
#
#This is the safe version of the performance-to-spoils transformation. Its derivative
#at P=0 is finite for every phi>0. In this first model phi=1, so W(P,1)=P exactly.

spoils_weight<-function(P,phi){
  (1+P)^phi-1
}

tours$E<-colSums(people$tour)
sum(tours$E)

tours$phi<-c(1,1)

people$spoils<-sapply(1:ntour,function(t){
  w<-spoils_weight(people$perf[,t],tours$phi[t])
  tours$E[t]*w/sum(w)
})

#Time-allocation entropy and total utility
people$ent<- -rowSums(ifelse(people$rel>0,
                             people$rel*log(people$rel),0))
people$utility<-rowSums(people$spoils)


#Analytical expectation for phi=1
#Because W(P,1)=P=alpha*g(x), the talent mean cancels exactly in the symmetric
#population expectation, as in the original baseline derivation.

mean_g<-function(B,beta){
  ((1+B)^(beta+1)-1)/(B*(beta+1))-1
}

expected_utility_r<-function(r,B,beta){
  mu_g<-mean_g(B,beta)
  B/(2*mu_g)*(g(B*r,beta)+g(B*(1-r),beta))
}

expected_utility_t<-function(t,B,beta){
  mu_g<-mean_g(B,beta)
  B/(2*mu_g)*(g(t,beta)+g(B-t,beta))
}


#Functions connected to entropy
entropy2<-function(r){
  ifelse(r==0 | r==1,0,-r*log(r)-(1-r)*log(1-r))
}

r_from_entropy<-function(h){
  if(h<=0) return(0)
  if(h>=log(2)) return(0.5)
  uniroot(function(r){entropy2(r)-h},c(0,0.5))$root
}

expected_utility_H<-function(h,B,beta){
  r<-vapply(h,r_from_entropy,1)
  expected_utility_r(r,B,beta)
}


#Check whether generalists on average beat specialists
layout(matrix(c(1,2),nrow=1),widths=c(1,1))
par(mar=c(3.5,3.5,1.5,0.5),mgp=c(2,0.65,0))

plot(people$ent,people$utility,
     pch=19,col="#00808020",
     ylim=c(0,max(people$utility)),
     xlab="Time allocation entropy",
     ylab="Total utility")

m<-lm(people$utility~people$ent)
abline(m,lwd=1,col=1)
text(range(people$ent),0,
     c("specialists","generalists"),
     pos=c(4,2),cex=0.8,offset=0)

hh<-seq(0,log(2),length.out=501)
lines(hh,
      expected_utility_H(hh,B=budget,beta=tours$beta[1]),
      lwd=lwd1,col=2)


#Same conclusion as a function of time invested in tournament 1
plot(people$tour[,1],people$utility,
     pch=19,col="#00808020",
     ylim=c(0,max(people$utility)),
     xlab="Time invested in tournament 1",
     ylab="Total utility")

x<-people$tour[,1]
m<-lm(people$utility~x+I(x^2))
xx<-seq(0,budget,length.out=501)
lines(xx,predict(m,newdata=data.frame(x=xx)),lwd=1,col=1)

text(c(0,budget),0,c("spec 2","spec 1"),pos=c(4,2),cex=0.8)
text(budget/2,0,"generalists",cex=0.8)

lines(xx,
      expected_utility_t(xx,B=budget,beta=tours$beta[1]),
      lwd=lwd1,col=2)

#The red curve is the analytical expectation, not a fitted regression.
