#I was considering gamma and lognormal

n<-10000
mu<-1
sd<-0.4
E<-(mu/sd)^2
gamma<-mu/sd^2
v<-rgamma(n,shape=E,rate=gamma)

mean(v)
sd(v)

plot(density(v))

#Here I compare lognormal with gamma distribution, what I do not like about gamma is how it apporaches exponential distribution for small sd at mu=1
mu<-1
sd<-1

#parameters gamma 
E<-(mu/sd)^2
gamma<-mu/sd^2

#parameters log-normal
m<-log((mu^2)/sqrt(mu^2+sd^2)) #both work, verified with 
s<-sqrt(log(1+(sd^2)/(mu^2)))

curve(dgamma(x,E,gamma),from=0,to=8,n=10001)
curve(dlnorm(x,m,s),col=2,n=10001,add=T)

#I take lognormal, but I will parametrize it via median and sd, that feels right, I might change this later
median<-1
sd<-1
m<-log(median)
s<-sqrt(log((1+sqrt(1+4*(sd/median)^2))/2))

v<-rlnorm(n,m,s)
mean(v)
median(v)
sd(v)

curve(dlnorm(x,m,s),col=2,n=10001,from=0,to=8)

max(v)
plot(density(v)) #there are few insanely tallented in that distribution

#Here is a function that calculates the expected maximum given parameters of the distribution and n
expected_max_lnorm <- function(n, meanlog = 0, sdlog = 1) {
  integrand <- function(z) {
    exp(log(n) + meanlog + sdlog * z + dnorm(z, log = TRUE) +
          (n - 1) * pnorm(z, log.p = TRUE))
  }
  integrate(integrand, lower = -Inf, upper = Inf, rel.tol = 1e-10)$value
}

expected_max_lnorm(n,m,s)
expected_max_lnorm(10^9,m,s) #One billion people
expected_max_lnorm(8*10^9,m,s) #Eight billion people

#Probably too much? I would try to restrict the sd a little
median<-1
sd<-0.5
m<-log(median)
s<-sqrt(log((1+sqrt(1+4*(sd/median)^2))/2))
curve(dlnorm(x,m,s),col=2,n=10001,from=0,to=8)
expected_max_lnorm(10^9,m,s) #One billion people
expected_max_lnorm(8*10^9,m,s) #Eight billion people

qlnorm(c(0.005,0.01,0.1,0.25,0.5,0.75,0.9,0.99,0.995),m, s) #quantiles 98% of people are between 0.36 and 2.74, the implied global maximum is quite large, but we do not need to worry about it right now

#Eventually, it might be a good idea to also try Gumbel distribution to compare it with the behavior of cumulative Henrich-like model
