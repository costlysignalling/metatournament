# Compiled continuous maximization for utility maps. Uses the same global
# objective as the population solver, refining every grid-detected local peak.
asymmetric_utility_values <- local({
  engine <- new.env(parent=globalenv())
  function(a1,a2,B,beta1,beta2,phi1,phi2,s,nx) {
    if(!exists('asymmetric_utilities_cpp',engine,inherits=FALSE)) {
      code <- readLines('05_energy_allocation_population_phispace.cpp',warn=FALSE)
      extra <- '
// [[Rcpp::export]]
NumericVector asymmetric_utilities_cpp(NumericVector a1,NumericVector a2,
 double B,double b1,double b2,double p1,double p2,double s,int nx,int threads) {
 int n=a1.size(); NumericVector out(n);
 ConditionalCell cell{b1,b2,p1,p2,B,1,s,nx};
#ifdef _OPENMP
 omp_set_num_threads(threads);
#pragma omp parallel for schedule(dynamic)
#endif
 for(int i=0;i<n;++i) {
  double aa=a1[i],bb=a2[i];
  if(p1<=1 && p2<=1) {
   double x=cell.optimum(0,B,aa,bb,nx);out[i]=cell.U(x,aa,bb);continue;
  }
  double step=B/(nx-1),prev=cell.U(0,aa,bb),cur=cell.U(step,aa,bb);
  double best=std::max(prev,cell.U(B,aa,bb));
  for(int k=1;k<nx-1;++k) {
   double next=cell.U((k+1)*step,aa,bb);
   if(cur>=prev && cur>=next) {
    double x=cell.golden((k-1)*step,(k+1)*step,aa,bb);
    best=std::max(best,cell.U(x,aa,bb));
   }
   best=std::max(best,cur);prev=cur;cur=next;
  }
  double x=cell.golden(0,step,aa,bb);best=std::max(best,cell.U(x,aa,bb));
  x=cell.golden(B-step,B,aa,bb);best=std::max(best,cell.U(x,aa,bb));
  out[i]=best;
 }
 return out;
}'
      Rcpp::sourceCpp(code=paste(c(code,extra),collapse='\n'),env=engine,
                      cacheDir='05_asymmetric_utility_cpp_cache')
    }
    if(length(a1)!=length(a2)) stop('Talent vectors must have equal lengths.')
    threads <- parallel::detectCores()
    if(is.na(threads)) threads <- 1L
    engine$asymmetric_utilities_cpp(a1,a2,B,beta1,beta2,phi1,phi2,s,
                                   as.integer(nx),max(1L,min(4L,threads-2L)))
  }
})
