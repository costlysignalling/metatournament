# Refine the phi1=phi2=1.5 percentile-contour neighbourhood deterministically.
# Conditional on alpha1, locate the global strategy switch across B/2 and
# integrate the two alpha2 intervals separately. This prevents a quadrature
# point from jumping between separated optima as the parameters vary.
# Outer integral: Gauss-Hermite; conditional log-talent integrals: Gauss-Legendre.
# All global utility comparisons retain endpoints and refined interior maxima.
Rcpp::sourceCpp('05_equal_talent_conditional_refinement.cpp',cacheDir='conditional_cpp_cache')
sd <- sqrt(log((1+sqrt(2))/2))
gh <- function(n) {
 J<-matrix(0,n,n);i<-1:(n-1);J[cbind(i,i+1)]<-sqrt(i/2)
 e<-eigen(J+t(J),symmetric=TRUE)
 list(a=exp(sqrt(2)*sd*e$values),w=e$vectors[1,]^2)
}
gl <- function(n) {
 J<-matrix(0,n,n);i<-1:(n-1);J[cbind(i,i+1)]<-i/sqrt(4*i*i-1)
 e<-eigen(J+t(J),symmetric=TRUE)
 list(x=e$values,w=2*e$vectors[1,]^2)
}
inputs <- readRDS('05_equal_talent_geometry_inputs.rds')
stopifnot(inputs$B==50,identical(inputs$phis,c(.66,1,1.5)))
b<-inputs$betas
# Keep the selected cells fixed between resumed runs. The band is deliberately
# much wider than the contour displacement from the original integration error.
selection_file <- '05_equal_talent_conditional_selection.rds'
if(file.exists(selection_file)) {
 selected <- readRDS(selection_file)
 stopifnot(identical(selected$betas,b))
 grid <- selected$grid
} else {
 grid<-expand.grid(i=seq_along(b),j=seq_along(b));grid<-subset(grid,i<j)
 alpha<-uniroot(function(total) integrate(function(a)
  dlnorm(a,0,sd)*plnorm(total-a,0,sd),0,total,rel.tol=1e-11)$value-.5,c(2,5))$root/2
 d<-b[grid$i]-b[grid$j];L<-1+inputs$B/2
 z<-log(inputs$ratios[[3,3]][cbind(grid$i,grid$j)])+log(b[grid$i]/b[grid$j])+d*log(L)+
  .5*(log1p(alpha*expm1(b[grid$i]*log(L)))-log1p(alpha*expm1(b[grid$j]*log(L))))
 grid<-grid[b[grid$i]>=.45 & b[grid$j]>=.7 & abs(z/d)<.05,,drop=FALSE]
 saveRDS(list(betas=b,grid=grid),selection_file)
}
q1<-gh(96);q2<-gl(96)
checkpoint_dir<-'equal_talent_conditional_checkpoints_v1'
dir.create(checkpoint_dir,showWarnings=FALSE)
chunks<-split(seq_len(nrow(grid)),grid$i)
results<-vector('list',length(chunks))
for(k in seq_along(chunks)) {
 ids<-chunks[[k]]
 file<-file.path(checkpoint_dir,sprintf('row_%03d_gh96_gl96_nx81.rds',grid$i[ids[1]]))
 record<-if(file.exists(file)) readRDS(file) else NULL
 if(!is.null(record)&&identical(record$j,grid$j[ids])) {
  results[[k]]<-record$result
 } else {
  results[[k]]<-conditional_equal_talent_nash_cpp(b[grid$i[ids]],b[grid$j[ids]],
   inputs$ratios[[3,3]][cbind(grid$i[ids],grid$j[ids])],q1$a,q1$w,q2$x,q2$w,
   inputs$B,sd,nx=81,threads=4)
  if(any(results[[k]]$status!=1)) {
   saveRDS(list(grid=grid[ids,],results=results[[k]]),'conditional_failed_row.rds')
   stop('Conditional integration did not converge; see conditional_failed_row.rds')
  }
  saveRDS(list(j=grid$j[ids],result=results[[k]]),file)
 }
 message('Switch-split conditional refinement row ',k,'/',length(chunks),' complete.')
}
z<-lapply(names(results[[1]]),function(name) {
 value<-numeric(nrow(grid));for(k in seq_along(chunks))value[chunks[[k]]]<-results[[k]][[name]];value
})
names(z)<-names(results[[1]])
rho<-inputs$ratios[[3,3]]
rho[cbind(grid$i,grid$j)]<-z$ratio;rho[cbind(grid$j,grid$i)]<-1/z$ratio;diag(rho)<-1
result<-list(betas=b,rho=rho,indices=grid,results=z,n_hermite=96L,n_legendre=96L,nx=81L,
 method='Deterministic conditional integration split at global strategy switch',
 standard_normal_limits=c(-10,10),selection_file=selection_file)
saveRDS(result,'05_equal_talent_conditional_refinement.rds')
inputs$ratios[[3,3]]<-rho
inputs$source$conditional_refinement<-result[c('method','n_hermite','n_legendre','nx','standard_normal_limits')]
saveRDS(inputs,'05_equal_talent_geometry_inputs.rds')
