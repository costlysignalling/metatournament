// Smooth conditional integration of rare bored time, using the existing
// globally searched individual best response. Includes the same Local class.
#include "08_boredom_core.cpp"
// [[Rcpp::export]]
Rcpp::List boredom_conditional_cpp(Rcpp::NumericVector background,
 double beta1,double beta2,double phi1,double phi2,double B,double compulsory,
 Rcpp::NumericVector a2,Rcpp::NumericVector outer_weights,
 Rcpp::NumericVector nodes,Rcpp::NumericVector inner_weights,double slog,
 int grid_half=80,int root_iterations=20){
 Local f(B,compulsory,beta1,beta2,phi1,phi2,grid_half);
 double total=0.;int failures=0;
 auto best=[&](double a,double b){return f.br(a,b,background[0],background[1],background[2],background[3]).first;};
 for(int j=0;j<a2.size();j++){
  auto boundary=[&](double target){
   // Integrate in probability space, covering all but <1e-12 of the talent tail.
   double lo=0.,hi=1.-1e-12;
   if(best(std::exp(slog*R::qnorm(hi,0,1,1,0)),a2[j])<target){failures++;return hi;}
   for(int it=0;it<root_iterations;it++){
    double mid=(lo+hi)/2.;
    if(best(std::exp(slog*R::qnorm(mid,0,1,1,0)),a2[j])>=target)hi=mid;else lo=mid;
   }return (lo+hi)/2.;
  };
  // Most of the rare-region integral is a constant cB below the q=0 boundary.
  // Integrate only the remaining partial-boredom strip numerically.
  double u0=boundary(1e-7),uc=boundary(compulsory*B-1e-7);
  if(uc<u0-2e-6){failures++;std::swap(u0,uc);}
  double val=compulsory*B*u0;
  if(uc-u0>1e-8)for(int k=0;k<nodes.size();k++){
   double u=u0+(uc-u0)*(nodes[k]+1.)/2.;
   double q=best(std::exp(slog*R::qnorm(u,0,1,1,0)),a2[j]);
   val+=(uc-u0)*inner_weights[k]/2.*std::max(compulsory*B-q,0.);
  }
  total+=outer_weights[j]*val;
 }
 return Rcpp::List::create(Rcpp::_["mean_bored"]=total,Rcpp::_["failures"]=failures);
}
