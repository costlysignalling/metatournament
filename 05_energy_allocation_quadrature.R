# Shared deterministic quadrature rules for the energy-allocation figure.
energy_gh <- function(n,sdlog) {
 J<-matrix(0,n,n);i<-1:(n-1);J[cbind(i,i+1)]<-sqrt(i/2)
 e<-eigen(J+t(J),symmetric=TRUE)
 list(a=exp(sqrt(2)*sdlog*e$values),w=e$vectors[1,]^2)
}
energy_gl <- function(n) {
 J<-matrix(0,n,n);i<-1:(n-1);J[cbind(i,i+1)]<-i/sqrt(4*i*i-1)
 e<-eigen(J+t(J),symmetric=TRUE)
 list(x=e$values,w=2*e$vectors[1,]^2)
}
energy_feedback_class <- function(Efirst,Eequilibrium,B=50,tol=1e-6) {
 d0<-Efirst-B/2;d1<-Eequilibrium-B/2
 out<-rep(0L,length(d0));active<-abs(d0)>tol
 same<-d0*d1>=0;stronger<-abs(d1)>abs(d0)
 out[active & same & stronger]<-1L
 out[active & same & !stronger]<-2L
 out[active & !same & !stronger]<-3L
 out[active & !same & stronger]<-4L
 out
}
