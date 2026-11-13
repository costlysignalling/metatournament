Rcpp::sourceCpp('05_equal_talent_conditional_refinement.cpp',cacheDir='conditional_cpp_cache')
sd<-sqrt(log((1+sqrt(2))/2))
gh<-function(n){J<-matrix(0,n,n);J[cbind(1:(n-1),2:n)]<-sqrt((1:(n-1))/2);e<-eigen(J+t(J),symmetric=TRUE);list(a=exp(sqrt(2)*sd*e$values),w=e$vectors[1,]^2)}
gl<-function(n){J<-matrix(0,n,n);i<-1:(n-1);J[cbind(i,i+1)]<-i/sqrt(4*i*i-1);e<-eigen(J+t(J),symmetric=TRUE);list(x=e$values,w=2*e$vectors[1,]^2)}
x<-readRDS('05_equal_talent_geometry_inputs.rds')
pars<-rbind(c(.5588,.9606),c(.5882,.9018),c(.6078,.8528),c(.6372,.8136),
            c(.6568,.7646),c(.7058,.7156),c(.7352,.794),c(.9704,.99))
initial<-sapply(seq_len(nrow(pars)),function(k)x$ratios[[3,3]][which.min(abs(x$betas-pars[k,1])),which.min(abs(x$betas-pars[k,2]))])
run<-function(n1,n2,nx){a<-gh(n1);b<-gl(n2);print(system.time(z<-conditional_equal_talent_nash_cpp(pars[,1],pars[,2],initial,a$a,a$w,b$x,b$w,50,sd,nx=nx,threads=3)));z}
a<-run(96,96,81);b<-run(192,192,161)
error<-abs(log(a$ratio/b$ratio))
print(data.frame(pars,relative_log_ratio_change=error))
stopifnot(all(a$status==1),all(b$status==1),max(error)<1e-5)
# At this almost wholly specialist point, the independently integrated
# endpoint control-variate calculation provides a further cross-check.
stopifnot(abs(log(b$ratio[nrow(pars)]/initial[nrow(pars)]))<1e-8)
saveRDS(list(pars=pars,a=a,b=b,error=error),'05_equal_talent_conditional_resolution_check.rds')
cat('Switch-split conditional integrals passed doubled-order and finer-grid checks.\n')
