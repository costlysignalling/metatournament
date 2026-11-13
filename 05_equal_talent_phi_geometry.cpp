// Included after 05_energy_allocation_population_phispace.cpp by the R loader.
// The log marginal ratio has at most one turning point, hence at most two roots.
// [[Rcpp::export]]
Rcpp::DataFrame equal_talent_phi_geometry_cpp(Rcpp::NumericVector p1,
 Rcpp::NumericVector p2,Rcpp::NumericVector rr,double b1,double b2,double B,
 bool check_global=false,int nx=401,int threads=4) {
 int n=p1.size();Rcpp::IntegerVector count(n),type(n);
 Rcpp::NumericVector first(n,NA_REAL),last(n,NA_REAL),lowout(n),highout(n);
 Rcpp::LogicalVector gf(n,NA_LOGICAL),gl(n,NA_LOGICAL);
 double L=1+B/2,G1=std::expm1(b1*std::log(L)),G2=std::expm1(b2*std::log(L));
 auto lpexp=[](double x){return x>0?x+std::log1p(std::exp(-x)):std::log1p(std::exp(x));};
#ifdef _OPENMP
 omp_set_num_threads(threads);
#pragma omp parallel for schedule(dynamic)
#endif
 for(int i=0;i<n;++i){
  double u=p1[i],v=p2[i],r=rr[i],tol=1e-11;
  double c=std::log(r)+std::log(u*b1/(v*b2))+(b1-b2)*std::log(L);
  double h=u==v?c+(u-1)*std::log(G1/G2):(u>v?INFINITY:-INFINITY);
  double A=(u-1)*G1-(v-1)*G2;
  lowout[i]=c;highout[i]=h;
  if(u==v && std::abs(c)<tol && (u==1 || G1==G2)) {type[i]=0;continue;}
  auto sign=[&](double x){return x>tol?1:(x< -tol?-1:0);};
  int sl=sign(c),sh=sign(h);if(sl==0)sl=sign(A);if(sh==0)sh=sl;
  if(sl==0)sl=sh;
  auto f=[&](double z){return c+(u-1)*lpexp(z+std::log(G1))-(v-1)*lpexp(z+std::log(G2));};
  std::vector<double> cuts{-1024.0};
  if(u!=v){double a=-A/((u-v)*G1*G2);if(std::isfinite(a)&&a>0)cuts.push_back(std::log(a));}
  cuts.push_back(1024.0);std::vector<double> roots;
  for(size_t k=0;k+1<cuts.size();++k){
   double lo=cuts[k],hi=cuts[k+1];int s0=k==0?sl:sign(f(lo));
   int s1=k+2==cuts.size()?sh:sign(f(hi));
   if(s0*s1>=0)continue;
   for(int it=0;it<65;++it){double mid=(lo+hi)/2;double fm=f(mid);
    if((fm>0?1:-1)==s0)lo=mid;else hi=mid;}
   roots.push_back((lo+hi)/2);
  }
  count[i]=roots.size();type[i]=roots.empty()?(sl>0?1:2):(sh>0?3:4);
  if(roots.empty())continue;
  first[i]=std::exp(roots.front());last[i]=std::exp(roots.back());
  if(check_global){
   ConditionalCell cell{b1,b2,u,v,B,1,r,nx};
   auto optimal=[&](double a) -> int {if(!std::isfinite(a))return NA_LOGICAL;
    double x=cell.optimum(0,B,a,a,nx),best=cell.U(x,a,a),mid=cell.U(B/2,a,a);
    return best<=mid+1e-8*std::max(1.0,std::abs(mid))?TRUE:FALSE;};
   gf[i]=optimal(first[i]);gl[i]=optimal(last[i]);
  }
 }
 return Rcpp::DataFrame::create(Rcpp::Named("phi1")=p1,Rcpp::Named("phi2")=p2,
 Rcpp::Named("ratio")=rr,Rcpp::Named("type")=type,Rcpp::Named("n_crossings")=count,
 Rcpp::Named("logR_low")=lowout,Rcpp::Named("logR_high")=highout,
 Rcpp::Named("alpha_first")=first,Rcpp::Named("alpha_cross")=last,
 Rcpp::Named("global_first")=gf,Rcpp::Named("global_equal_split")=gl);
}

