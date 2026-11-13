for(e in parse('06B_predict_many_tournaments.R'))
  if(is.call(e) && identical(e[[1]],as.name('<-')))eval(e)
z <- readRDS('06B_many_results/three_tournament_population.rds')
reference128 <- many_random_reference(z$beta,z$phi,z$B,nquad=128L)
stopifnot(max(abs(log(reference128$M/z$fixed$M)))<1e-8)
A <- many_lognormal_population(1024,3)
fit <- many_allocate(A,z$beta,z$phi,z$B,tail(z$checks,1)[[1]]$responsive$scaler)
stopifnot(all(fit$within_tolerance),
  max(abs(rowSums(fit$allocation)-z$B))<1e-8,min(fit$allocation)>-1e-9)
# Independent analytic marginal-return check on every active/excluded activity.
X <- fit$allocation;C <- fit$scaler;D <- X
for(t in 1:3)D[,t] <- C[t]*z$phi[t]*A[,t]*z$beta[t]*(1+X[,t])^(z$beta[t]-1)*
  (1+A[,t]*((1+X[,t])^z$beta[t]-1))^(z$phi[t]-1)
for(i in seq_len(nrow(X))) {
  active <- X[i,]>1e-6;lambda <- mean(D[i,active])
  stopifnot(max(abs(D[i,active]-lambda))<1e-7)
  if(any(!active))stopifnot(max(D[i,!active])<lambda+1e-7)
}
# Equal linear-spoils tournaments have a known exact active-set solution.
for(T in 2:10) {
  alpha <- exp(seq(-.7,.7,length.out=T));beta <- .6;v <- alpha^(1/(1-beta));o <- order(v,decreasing=TRUE)
  k <- max(which((50+seq_len(T))*v[o]>cumsum(v[o])))
  expected <- numeric(T);expected[o[seq_len(k)]] <- (50+k)*v[o[seq_len(k)]]/sum(v[o[seq_len(k)]])-1
  ans <- many_allocate(alpha,rep(beta,T),rep(1,T),50)
  stopifnot(max(abs(ans$allocation-expected))<2e-6)
}
# Independently enumerate a triangular 500-step simplex in base R.
a <- c(.457,.457,.457);n <- 500L;best <- -Inf
w <- function(x) (1+.457*((1+x)^.5-1))^3-1
for(i in 0:n) {
  j <- 0:(n-i);value <- w(i/n)+w(j/n)+w(1-(i+j)/n)
  best <- max(best,max(value))
}
ans <- many_allocate(a,rep(.5,3),rep(3,3),1,tolerance=1e-6,max_grid=2048L)
stopifnot(best<=ans$utility+ans$utility_gap_bound+1e-10,
          ans$utility>.68230,all(ans$allocation>0),all(ans$within_tolerance))
cat('Final validation: marginal conditions, budgets, exact T=2..10 solutions, and independent nonconcave simplex enumeration passed.\n')
