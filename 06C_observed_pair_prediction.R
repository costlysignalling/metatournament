# Exact conditional prediction with unobserved outside tournaments.
# Use only the two observed talents and their combined budget b.
source('06B_predict_many_tournaments.R')

predict_observed_pair <- function(alpha,b,beta,phi,weights=NULL,
  reference_scaler=NULL,start=c(1,1),...) {
  alpha <- as.matrix(alpha)
  if(ncol(alpha)!=2 || length(beta)!=2 || length(phi)!=2)
    stop('Supply exactly two observed tournaments.')
  if(!is.null(reference_scaler)) {
    # First response to an independently known initial population reference.
    return(many_allocate(alpha,beta,phi,b,reference_scaler,...))
  }
  # Responsive conditional prediction. The measured b is held fixed for each
  # actor; no observed within-pair splits are used to obtain this solution.
  # Passing max_grid/allocation_tolerance is possible through ... here.
  fit <- many_responsive(alpha,beta,phi,b,start=start,weights=weights,...)
  if(!fit$converged)stop('Conditional consistency equation did not converge.')
  fit
}

pair_regret <- function(alpha,b,observed_x1,beta,phi,scaler) {
  prediction <- many_allocate(alpha,beta,phi,b,scaler)
  if(any(observed_x1<0 | observed_x1>b))stop('Observed allocations must lie in [0,b].')
  observed <- cbind(observed_x1,b-observed_x1)
  W <- observed
  for(t in 1:2)W[,t] <- expm1(phi[t]*log1p(alpha[,t]*expm1(beta[t]*log1p(observed[,t]))))
  value <- as.vector(W%*%prediction$scaler)
  # For a grid solution, true regret lies in [lower, upper].
  data.frame(lower=pmax(0,prediction$utility-value),
    upper=pmax(0,prediction$utility+prediction$utility_gap_bound-value))
}

# Reader example with a supplied CSV (replace names as necessary):
# dat <- read.csv('observed_pair.csv')
# alpha <- as.matrix(dat[c('alpha1','alpha2')])
# fit <- predict_observed_pair(alpha,dat$pair_budget,
#   beta=c(.6,.4),phi=c(.66,1),weights=dat$sampling_weight)
# dat$predicted_x1 <- fit$allocation[,1]
# dat$predicted_x2 <- fit$allocation[,2]
# write.csv(dat,'observed_pair_predictions.csv',row.names=FALSE)
# many_pie(list(allocation=fit$allocation[1,,drop=FALSE]),'Actor 1: conditional prediction')
#
# To test an observed responsive population directly, calculate its scalers:
# observed_X <- cbind(dat$x1,dat$pair_budget-dat$x1)
# measured <- many_means(alpha,observed_X,c(.6,.4),c(.66,1),dat$sampling_weight)
# regret <- pair_regret(alpha,dat$pair_budget,dat$x1,c(.6,.4),c(.66,1),measured$scaler)
# For a first-response test use the INITIAL measured population's scalers,
# not these post-response observed means.
