// Model evaluations for direct transition-curve tracing.
#include "06_coexistence_hdr_scores.cpp"
// [[Rcpp::plugins(cpp11)]]
// [[Rcpp::plugins(openmp)]]
// [[Rcpp::export]]
NumericVector coex_random_moments_cpp(NumericVector beta,NumericVector phi,
 NumericVector aa,NumericVector aw,NumericVector xx,NumericVector ww,double B,int threads=3){
 int n=beta.size();NumericVector out(n);
#ifdef _OPENMP
 omp_set_num_threads(threads);
#pragma omp parallel for schedule(dynamic)
#endif
 for(int k=0;k<n;++k){double s=0;
  for(int i=0;i<aa.size();++i)for(int j=0;j<xx.size();++j)
   s+=aw[i]*ww[j]*wval(B*(1+xx[j])/2,aa[i],beta[k],phi[k])/2;
  out[k]=s;
 };return out;
}
// [[Rcpp::export]]
NumericVector coex_selected_scores_cpp(NumericVector beta1,NumericVector beta2,
 NumericVector phi1,NumericVector phi2,NumericVector rho,IntegerVector field,
 NumericVector xt,NumericVector xl,NumericVector xr,NumericVector theta,
 double hc,double hr,double B,int threads=3){
 int n=beta1.size();NumericVector out(n);
#ifdef _OPENMP
 omp_set_num_threads(threads);
#pragma omp parallel for schedule(dynamic)
#endif
 for(int k=0;k<n;++k){double sm,sd;
  if(field[k]<=2)scores_one(beta1[k],beta2[k],phi1[k],phi2[k],rho[k],xt,xl,xr,theta,hc,hr,B,sm,sd);
  else scores_one(beta2[k],beta1[k],phi2[k],phi1[k],1/rho[k],xt,xl,xr,theta,hc,hr,B,sm,sd);
  out[k]=(field[k]==1 || field[k]==4)?sm:sd;
 };return out;
}
