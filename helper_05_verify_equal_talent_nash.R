Rcpp::sourceCpp('05_equal_talent_nash_refinement.cpp',cacheDir='equal_talent_rcpp_cache_v2')
sdlog<-sqrt(log((1+sqrt(2))/2))
radical<-function(n,b) {out<-numeric(length(n));f<-1/b;while(any(n>0)) {
  out<-out+f*(n%%b);n<-n%/%b;f<-f/b};out}
gh<-function(n) {J<-matrix(0,n,n);J[cbind(1:(n-1),2:n)]<-sqrt((1:(n-1))/2)
  e<-eigen(J+t(J),symmetric=TRUE);list(a=exp(sqrt(2)*sdlog*e$values),w=e$vectors[1,]^2)}
inputs<-readRDS('05_equal_talent_geometry_inputs.rds')
pars<-rbind(c(.7058,.7156),c(.7058,.7842),c(.7254,.8332),c(.745,.8332),
            c(.794,.8038),c(.8136,.9018),c(.9508,.99),c(.9704,.99))
initial<-vapply(seq_len(nrow(pars)),function(i) {
  j<-which.min(abs(inputs$betas-pars[i,1]));k<-which.min(abs(inputs$betas-pars[i,2]))
  inputs$ratios[[3,3]][j,k]
},1)
run<-function(n,ngh,nx,fine=FALSE) {
  q<-gh(ngh)
  result <- vector('list',2)
  groups <- list(which(pars[,1]<.73),which(pars[,1]>=.73))
  for(k in 1:2) {
    ids <- groups[[k]]
    points <- if(k==1) {if(fine) 131072L else 65536L} else n
    grid <- if(k==1) 161L else nx
    ii<-64+seq_len(points)
    result[[k]] <- refine_equal_talent_nash_cpp(pars[ids,1],pars[ids,2],initial[ids],
      radical(ii,2),radical(ii,3),q$a,q$w,50,sdlog,nx=grid,threads=2)
  }
  out <- lapply(names(result[[1]]),function(name) {
    value <- numeric(nrow(pars))
    for(k in 1:2) value[groups[[k]]] <- result[[k]][[name]]
    value
  })
  names(out)<-names(result[[1]])
  out
}
a<-run(16384,192,81);b<-run(65536,256,161,fine=TRUE)
stopifnot(all(a$status %in% c(1,2)),all(b$status %in% c(1,2)))
error<-abs(log(a$ratio/b$ratio))
print(data.frame(pars,relative_log_ratio_change=error,
                 coarse_status=a$status,fine_status=b$status))
stopifnot(max(error[pars[,1]>=.73])<1e-5,
          max(error[pars[,1]<.73])<1e-4)
# Independent adaptive integration of the continuous endpoint-only expectations.
q<-gh(128);p<-refine_equal_talent_nash_cpp(.9704,.99,initial[nrow(pars)],numeric(),numeric(),
 q$a,q$w,50,sdlog,threads=1)$ratio
g1<-51^.9704-1;g2<-51^.99-1
quad<-function(f) integrate(function(z) f(exp(sdlog*z))*dnorm(z),-12,12,rel.tol=1e-10)$value
w1<-function(a) (1+a*g1)^1.5-1;w2<-function(a) (1+a*g2)^1.5-1
t1<-function(a) ((1+w2(a)/p)^(1/1.5)-1)/g1
t2<-function(a) ((1+p*w1(a))^(1/1.5)-1)/g2
E1<-50*quad(function(a) pnorm(log(t1(a))/sdlog,lower.tail=FALSE))
M1<-quad(function(a) w1(a)*pnorm(log(t2(a))/sdlog))
M2<-quad(function(a) w2(a)*pnorm(log(t1(a))/sdlog))
stopifnot(abs(log((E1/M1)/((50-E1)/M2)/p))<1e-8)
saveRDS(list(parameters=pars,coarse=a,fine=b,error=error),
        '05_equal_talent_nash_resolution_check.rds')
cat('Refined Nash ratios passed resolution and independent adaptive-integration checks.\n')
