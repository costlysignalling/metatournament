// Deterministic population integration for phi1=phi2=1.5.
// At each alpha1, split the alpha2 integral at the globally optimal switch
// across B/2. Integrating each side separately removes quadrature-node jumps.
#include <Rcpp.h>
#include <cmath>
#include <algorithm>
#ifdef _OPENMP
#include <omp.h>
#endif
using namespace Rcpp;
// [[Rcpp::plugins(cpp11)]]
// [[Rcpp::plugins(openmp)]]
struct ConditionalState {double e,m1,m2;};
struct ConditionalCell {
 double b1,b2,B,sd,rho;
 int nx;
 double w(double x,double a,double b) const {
  return std::expm1(1.5*std::log1p(a*std::expm1(b*std::log1p(x))));
 }
 double U(double x,double a1,double a2) const {return rho*w(x,a1,b1)+w(B-x,a2,b2);}
 double D(double x,double a1,double a2) const {
  return rho*1.5*a1*b1*std::pow(1+x,b1-1)*std::sqrt(1+a1*std::expm1(b1*std::log1p(x)))-
   1.5*a2*b2*std::pow(1+B-x,b2-1)*std::sqrt(1+a2*std::expm1(b2*std::log1p(B-x)));
 }
 double golden(double lo,double hi,double a1,double a2) const {
  double gr=.6180339887498948482,c=hi-gr*(hi-lo),d=lo+gr*(hi-lo);
  double fc=U(c,a1,a2),fd=U(d,a1,a2);
  for(int i=0;i<32;++i) {
   if(fc<fd){lo=c;c=d;fc=fd;d=lo+gr*(hi-lo);fd=U(d,a1,a2);}
   else {hi=d;d=c;fd=fc;c=hi-gr*(hi-lo);fc=U(c,a1,a2);}
  }
  return (lo+hi)/2;
 }
 double optimum(double lo,double hi,double a1,double a2,int n) const {
  double step=(hi-lo)/(n-1),bestx=lo,best=U(lo,a1,a2);
  int bestk=0;
  for(int k=1;k<n;++k) {double x=lo+step*k,v=U(x,a1,a2);if(v>best){best=v;bestx=x;bestk=k;}}
  if(bestk>0 && bestk<n-1) {
   double x=golden(lo+step*(bestk-1),lo+step*(bestk+1),a1,a2),v=U(x,a1,a2);
   if(v>best){best=v;bestx=x;}
  }
  if(D(lo,a1,a2)>0) {double x=golden(lo,lo+step,a1,a2),v=U(x,a1,a2);if(v>best){best=v;bestx=x;}}
  if(D(hi,a1,a2)<0) {double x=golden(hi-step,hi,a1,a2),v=U(x,a1,a2);if(v>best){best=v;bestx=x;}}
  return bestx;
 }
 double side_difference(double z,double a1) const {
  double a2=std::exp(sd*z);
  double left=optimum(0,B/2,a1,a2,21),right=optimum(B/2,B,a1,a2,21);
  return U(right,a1,a2)-U(left,a1,a2);
 }
 double switch_z(double a1) const {
  double lo=-10,hi=10;
  if(side_difference(lo,a1)<=0) return lo;
  if(side_difference(hi,a1)>=0) return hi;
  for(int i=0;i<37;++i){double mid=(lo+hi)/2;if(side_difference(mid,a1)>0)lo=mid;else hi=mid;}
  return (lo+hi)/2;
 }
 ConditionalState state(const NumericVector& aa,const NumericVector& aw,
                     const NumericVector& nodes,const NumericVector& weights) const {
  double E=0,M1=0,M2=0;
  for(int i=0;i<aa.size();++i) {
   if(std::abs(std::log(aa[i])/sd)>10) continue;
   double cut=switch_z(aa[i]);
   for(int side=0;side<2;++side) {
    double lo=side==0?-10:cut,hi=side==0?cut:10;
    if(hi-lo<1e-12)continue;
    for(int j=0;j<nodes.size();++j) {
     double z=(lo+hi)/2+(hi-lo)*nodes[j]/2,a2=std::exp(sd*z);
     double pw=aw[i]*weights[j]*(hi-lo)/2*std::exp(-z*z/2)/std::sqrt(2*M_PI);
     // Global comparison is retained even within each split interval.
     double x=optimum(0,B,aa[i],a2,nx);
     E+=pw*x;M1+=pw*w(x,aa[i],b1);M2+=pw*w(B-x,a2,b2);
    }
   }
  }
  return ConditionalState{E,M1,M2};
 }
};
// [[Rcpp::export]]
List conditional_equal_talent_nash_cpp(NumericVector beta1,NumericVector beta2,
 NumericVector initial,NumericVector aa,NumericVector aw,NumericVector nodes,
 NumericVector weights,double B,double sd,int nx=81,int threads=4) {
 int n=beta1.size();NumericVector ratio(n),residual(n);IntegerVector status(n);
#ifdef _OPENMP
 omp_set_num_threads(threads);
#pragma omp parallel for schedule(dynamic)
#endif
 for(int i=0;i<n;++i) {
  ConditionalCell cell{beta1[i],beta2[i],B,sd,initial[i],nx};
  auto f=[&](double z) {cell.rho=std::exp(z);ConditionalState st=cell.state(aa,aw,nodes,weights);
   return std::log(st.e)+std::log(st.m2)-std::log(B-st.e)-std::log(st.m1)-z;};
  double centre=std::log(initial[i]),lo=centre-.2,hi=centre+.2,flo=f(lo),fhi=f(hi);
  if(flo*fhi>0){lo=centre-2;hi=centre+2;flo=f(lo);fhi=f(hi);}
  if(!std::isfinite(flo)||!std::isfinite(fhi)||flo*fhi>0){status[i]=0;continue;}
  double root=centre;int previous=0;
  for(int k=0;k<60;++k){
   double c=(lo*fhi-hi*flo)/(fhi-flo);if(!std::isfinite(c)||c<=lo||c>=hi)c=(lo+hi)/2;
   double fc=f(c);root=c;if(std::abs(fc)<1e-10||hi-lo<1e-10)break;
   if(flo*fc>0){lo=c;flo=fc;if(previous==1)fhi*=.5;previous=1;}
   else{hi=c;fhi=fc;if(previous==-1)flo*=.5;previous=-1;}
  }
  ratio[i]=std::exp(root);residual[i]=f(root);status[i]=std::abs(residual[i])<1e-8?1:0;
 }
 return List::create(_["ratio"]=ratio,_["residual"]=residual,_["status"]=status);
}
