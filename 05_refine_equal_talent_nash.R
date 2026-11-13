# Remove discrete winner-count noise throughout the beta1,beta2>=.70 corner
# for phi1=phi2=1.5, including the old refinement boundary at beta=.8.
# Endpoint-only expectations are smooth 1D integrals. A stratified QMC
# residual accounts for genuine diversification, rather than assuming it away.
if(!requireNamespace('Rcpp',quietly=TRUE)) stop('Rcpp is needed to rebuild this refinement.')
dir.create('equal_talent_rcpp_cache_v2',showWarnings=FALSE)
Rcpp::sourceCpp('05_equal_talent_nash_refinement.cpp',cacheDir='equal_talent_rcpp_cache_v2')
inputs <- readRDS('05_equal_talent_geometry_inputs.rds')
stopifnot(inputs$B==50,inputs$talent$meanlog==0,
          inputs$talent$sd_original==.5,inputs$talent$median==1,
          identical(inputs$phis,c(.66,1,1.5)))
radical <- function(n,base) {
  ans<-numeric(length(n));f<-1/base
  while(any(n>0)) {ans<-ans+f*(n%%base);n<-n%/%base;f<-f/base};ans
}
gh <- function(n=128L) {
  J<-matrix(0,n,n);J[cbind(1:(n-1),2:n)]<-sqrt((1:(n-1))/2);J<-J+t(J)
  e<-eigen(J,symmetric=TRUE)
  list(a=exp(sqrt(2)*sqrt(log((1+sqrt(2))/2))*e$values),w=e$vectors[1,]^2)
}
n_conditional <- 16384L
n_conditional_boundary <- 65536L
ii<-64L+seq_len(n_conditional_boundary);u<-radical(ii,2);v<-radical(ii,3)
q<-gh(192L)
grid<-expand.grid(i=seq_along(inputs$betas),j=seq_along(inputs$betas))
grid<-subset(grid,i<j & inputs$betas[i]>=.70 & inputs$betas[j]>=.70)
# One checkpoint per beta row, so an interrupted rebuild can resume.
checkpoint_dir <- 'equal_talent_nash_continuous_checkpoints'
dir.create(checkpoint_dir,showWarnings=FALSE)
chunks <- split(seq_len(nrow(grid)),grid$i)
results <- vector('list',length(chunks))
for(k in seq_along(chunks)) {
  ids <- chunks[[k]]
  boundary <- inputs$betas[grid$i[ids[1]]]<.73
  row_n <- if(boundary) n_conditional_boundary else n_conditional
  row_nx <- if(boundary) 161L else 81L
  file <- file.path(checkpoint_dir,sprintf('row_%03d_n%d_gh192_nx%d.rds',grid$i[ids[1]],row_n,row_nx))
  if(file.exists(file)) {
    results[[k]] <- readRDS(file)
  } else {
    results[[k]] <- refine_equal_talent_nash_cpp(inputs$betas[grid$i[ids]],inputs$betas[grid$j[ids]],
      inputs$ratios[[3,3]][cbind(grid$i[ids],grid$j[ids])],u[seq_len(row_n)],v[seq_len(row_n)],q$a,q$w,
      inputs$B,sqrt(log((1+sqrt(2))/2)),nx=row_nx,threads=4)
    if(any(!results[[k]]$status %in% c(1,2))) {
      saveRDS(list(indices=grid[ids,],results=results[[k]]),'failed_row.rds')
      print(data.frame(beta1=inputs$betas[grid$i[ids]],beta2=inputs$betas[grid$j[ids]],
            status=results[[k]]$status,residual=results[[k]]$residual))
      stop('Numerical check failed; see failed_row.rds')
    }
    saveRDS(results[[k]],file)
  }
  message('Continuous Nash refinement row ',k,'/',length(chunks),' complete.')
}
# Restore the original expand.grid order after the row-wise checkpoints.
z <- lapply(names(results[[1]]),function(name) {
  out <- numeric(nrow(grid))
  for(k in seq_along(chunks)) out[chunks[[k]]] <- results[[k]][[name]]
  out
})
names(z) <- names(results[[1]])
stopifnot(all(z$status %in% c(1,2)),
          all(z$log_bracket_width[z$status==2]<1e-10))
rho<-inputs$ratios[[3,3]]
rho[cbind(grid$i,grid$j)]<-z$ratio
rho[cbind(grid$j,grid$i)]<-1/z$ratio
diag(rho)<-1
result<-list(betas=inputs$betas,rho=rho,indices=grid,results=z,
             n_conditional=c(main=n_conditional,boundary=n_conditional_boundary),
             n_hermite=192L,nx=c(main=81L,boundary=161L),
             high_accuracy_boundary=.73,
             method='Continuous endpoint moments + stratified diversification residual',
             beta_min=.70,root_solver='Bracketed Illinois; smooth roots or sign-changing brackets < 1e-10')
saveRDS(result,'05_equal_talent_nash_refinement.rds')
inputs$ratios[[3,3]]<-rho
inputs$source$high_beta_refinement<-result[c('method','n_conditional','n_hermite','nx','beta_min','root_solver','high_accuracy_boundary')]
saveRDS(inputs,'05_equal_talent_geometry_inputs.rds')
