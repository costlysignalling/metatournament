#include <Rcpp.h>
#include <vector>
#include <cmath>
#include <algorithm>
using namespace Rcpp;
// [[Rcpp::plugins(cpp11)]]

struct Return {
  double a,b,p,c;
  double w(double x) const {return std::expm1(p*std::log1p(a*std::expm1(b*std::log1p(x))));}
  double value(double x) const {return c*w(x);}
  double slope(double x) const {
    if(c==0 || a==0) return 0;
    return c*p*a*b*std::exp((b-1)*std::log1p(x)+(p-1)*std::log1p(a*std::expm1(b*std::log1p(x))));
  }
};

template<class F> double bisect(F f,double lo,double hi) {
  double fl=f(lo);
  for(int i=0;i<65;++i) {
    double mid=(lo+hi)/2, fm=f(mid);
    if((fl>0)==(fm>0)) {lo=mid;fl=fm;} else hi=mid;
  }
  return (lo+hi)/2;
}

// The sign of w'' changes at most once. Split there, then find ALL
// stationary points of c*w(x)-eta*x on its monotone-slope intervals.
std::vector<double> single_candidates(const Return& r,double B,double eta) {
  std::vector<double> cuts{0,B}, out{0,B};
  double den=r.a*(r.b*r.p-1);
  if(den!=0) {
    double y=-(r.b-1)*(1-r.a)/den;
    if(y>0) {
      double x=std::expm1(std::log(y)/r.b);
      if(x>0 && x<B) cuts.push_back(x);
    }
  }
  std::sort(cuts.begin(),cuts.end());
  auto f=[&](double x){return r.slope(x)-eta;};
  for(size_t i=1;i<cuts.size();++i) {
    double lo=cuts[i-1],hi=cuts[i],fl=f(lo),fh=f(hi);
    if(fl*fh<0) out.push_back(bisect(f,lo,hi));
    if(fl==0) out.push_back(lo);
    if(fh==0) out.push_back(hi);
  }
  return out;
}

// [[Rcpp::export]]
NumericMatrix optout_alloc_cpp(NumericMatrix talent,NumericVector beta,
    NumericVector phi,NumericVector scaler,double B,double eta,int nscan=129) {
  int n=talent.nrow();
  NumericMatrix ans(n,4);
  colnames(ans)=CharacterVector::create("x0","x1","x2","utility");
  // Fixed scan plus logarithmic spacing near the two endpoints. Every
  // bracketed boundary stationary point is refined; endpoints are included.
  std::vector<double> grid;
  for(int k=0;k<nscan;++k) grid.push_back(B*k/(nscan-1));
  for(int k=1;k<=30;++k) {
    double x=B*std::pow(0.6,k);
    grid.push_back(x);grid.push_back(B-x);
  }
  std::sort(grid.begin(),grid.end());
  grid.erase(std::unique(grid.begin(),grid.end()),grid.end());
  for(int i=0;i<n;++i) {
    Return r1{talent(i,0),beta[0],phi[0],scaler[0]};
    Return r2{talent(i,1),beta[1],phi[1],scaler[1]};
    double best=eta*B,bx=0,by=0;
    auto consider=[&](double x,double y) {
      if(x+y>B+1e-10) return;
      double u=r1.value(x)+r2.value(y)+eta*std::max(0.0,B-x-y);
      if(u>best+1e-12*std::max(1.0,std::fabs(best))) {best=u;bx=x;by=y;}
    };
    auto s1=single_candidates(r1,B,eta),s2=single_candidates(r2,B,eta);
    for(double x:s1) for(double y:s2) consider(x,y);
    if(scaler[0]>0 || scaler[1]>0) {
      auto f=[&](double x){return r1.slope(x)-r2.slope(B-x);};
      double lo=grid[0],fl=f(lo);
      consider(0,B);consider(B,0);
      for(size_t k=1;k<grid.size();++k) {
        double hi=grid[k],fh=f(hi);
        // Each sampled point is also a global-search safeguard.
        consider(hi,B-hi);
        if(fl*fh<0) {double x=bisect(f,lo,hi);consider(x,B-x);}
        if(fl==0) consider(lo,B-lo);
        lo=hi;fl=fh;
      }
    }
    ans(i,0)=std::max(0.0,B-bx-by);ans(i,1)=bx;ans(i,2)=by;ans(i,3)=best;
  }
  return ans;
}
