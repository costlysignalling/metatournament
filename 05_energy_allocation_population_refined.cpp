// Higher-accuracy population energies, splitting at ALL specialist departures
// as well as the switch across B/2. The integrand always uses the global optimum.
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
 double b1,b2,p1,p2,B,sd,rho;
 int nx;
 double w(double x,double a,double b,double p) const {
  return std::expm1(p*std::log1p(a*std::expm1(b*std::log1p(x))));
 }
 double U(double x,double a1,double a2) const {return rho*w(x,a1,b1,p1)+w(B-x,a2,b2,p2);}
 double D(double x,double a1,double a2) const {
  return rho*p1*a1*b1*std::pow(1+x,b1-1)*std::pow(1+a1*std::expm1(b1*std::log1p(x)),p1-1)-
   p2*a2*b2*std::pow(1+B-x,b2-1)*std::pow(1+a2*std::expm1(b2*std::log1p(B-x)),p2-1);
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
  if(p1<=1 && p2<=1) {
   if(D(lo,a1,a2)<=0)return lo;
   if(D(hi,a1,a2)>=0)return hi;
   for(int k=0;k<30;++k){double mid=(lo+hi)/2;if(D(mid,a1,a2)>0)lo=mid;else hi=mid;}
   return (lo+hi)/2;
  }
  if(b1*p1>1 && b2*p2>1 && a1>=(1-b1)/(b1*(p1-1)) && a2>=(1-b2)/(b2*(p2-1)))
   return U(lo,a1,a2)>U(hi,a1,a2)?lo:hi;
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
  if(p1<=1 && p2<=1)return D(B/2,a1,a2);
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
 // Inverse of the monotone global best response. Besides the B/2
 // switch, split at departure from each specialist endpoint; a jump
 // from an endpoint to a diversifier need not cross B/2.
 double endpoint_cut(double a1,double target) const {
  double lo=-10,hi=10;
  if(optimum(0,B,a1,std::exp(sd*lo),41)<=target)return lo;
  if(optimum(0,B,a1,std::exp(sd*hi),41)>=target)return hi;
  for(int i=0;i<37;++i){double mid=(lo+hi)/2;
   if(optimum(0,B,a1,std::exp(sd*mid),41)>target)lo=mid;else hi=mid;}
  return (lo+hi)/2;
 }
 // In the concave case the derivative condition defines a unique talent
 // cutoff. Fubini's theorem then gives smooth energy and weight integrals,
 // without solving an allocation problem at every quadrature point.
 double inverse_log_marginal(double logq,double b,double p,double x) const {
  double G=std::expm1(b*std::log1p(x));
  double target=logq-std::log(p*b)-(b-1)*std::log1p(x),z=target;
  for(int k=0;k<14;++k){
   double ag=std::exp(z)*G;
   double error=z+(p-1)*std::log1p(ag)-target;
   double step=error*(1+ag)/(1+p*ag);z-=step;
   if(std::abs(step)<1e-12)break;
  }
  return z;
 }
 ConditionalState state(const NumericVector& aa,const NumericVector& aw,
                     const NumericVector& nodes,const NumericVector& weights) const {
  double E=0,M1=0,M2=0;
  if(p1==1 && p2==1) {
   double sigma=std::sqrt(2.0)*sd,mu=std::exp(sd*sd/2);
   for(int j=0;j<nodes.size();++j){
    double x=B*(nodes[j]+1)/2,pw=B*weights[j]/2;
    double cut=std::log(b2/b1/rho)+(b2-1)*std::log1p(B-x)+(1-b1)*std::log1p(x);
    E+=pw*R::pnorm(cut/sigma,0,1,false,false);
    M1+=pw*mu*b1*std::pow(1+x,b1-1)*R::pnorm((cut-sd*sd)/sigma,0,1,false,false);
    M2+=pw*mu*b2*std::pow(1+B-x,b2-1)*R::pnorm((cut+sd*sd)/sigma,0,1,true,false);
   }
   return ConditionalState{E,M1,M2};
  }
  if(p1<=1 && p2<=1) {
   for(int j=0;j<nodes.size();++j){
    double x=B*(nodes[j]+1)/2,pw=B*weights[j]/2;
    for(int i=0;i<aa.size();++i){
     double a=aa[i];
     double l1=std::log(p1*a*b1)+(b1-1)*std::log1p(x)+
       (p1-1)*std::log1p(a*std::expm1(b1*std::log1p(x)));
     double l2=std::log(p2*a*b2)+(b2-1)*std::log1p(B-x)+
       (p2-1)*std::log1p(a*std::expm1(b2*std::log1p(B-x)));
     double P2=R::pnorm(inverse_log_marginal(std::log(rho)+l1,b2,p2,B-x)/sd,0,1,true,false);
     double P1=R::pnorm(inverse_log_marginal(l2-std::log(rho),b1,p1,x)/sd,0,1,true,false);
     E+=pw*aw[i]*P2;M1+=pw*aw[i]*std::exp(l1)*P2;M2+=pw*aw[i]*std::exp(l2)*P1;
    }
   }
   return ConditionalState{E,M1,M2};
  }
  for(int i=0;i<aa.size();++i) {
   if(std::abs(std::log(aa[i])/sd)>10) continue;
   std::vector<double> cuts{-10,endpoint_cut(aa[i],B*(1-1e-5)),
     switch_z(aa[i]),endpoint_cut(aa[i],B*1e-5),10};
   std::sort(cuts.begin(),cuts.end());
   cuts.erase(std::unique(cuts.begin(),cuts.end(),[](double a,double b){return std::abs(a-b)<1e-7;}),cuts.end());
   for(size_t side=0;side+1<cuts.size();++side) {
    double lo=cuts[side],hi=cuts[side+1];
    if(hi-lo<1e-12)continue;
    for(int j=0;j<nodes.size();++j) {
     double z=(lo+hi)/2+(hi-lo)*nodes[j]/2,a2=std::exp(sd*z);
     double pw=aw[i]*weights[j]*(hi-lo)/2*std::exp(-z*z/2)/std::sqrt(2*M_PI);
     // Global comparison is retained even within each split interval.
     double x=optimum(0,B,aa[i],a2,nx);
     E+=pw*x;M1+=pw*w(x,aa[i],b1,p1);M2+=pw*w(B-x,a2,b2,p2);
    }
   }
  }
  return ConditionalState{E,M1,M2};
 }
};

// [[Rcpp::export]]
List energy_population_cpp(NumericVector beta1,NumericVector beta2,NumericVector phi1,
 NumericVector phi2,NumericVector fixed_ratio,NumericVector initial,
 NumericVector aa,NumericVector aw,NumericVector nodes,NumericVector weights,
 double B,double sd,int nx=61,int threads=4) {
 int n=beta1.size();NumericVector first(n),equilibrium(n),ratios(n),residual(n);
 IntegerVector status(n),iterations(n);
#ifdef _OPENMP
 omp_set_num_threads(threads);
#pragma omp parallel for schedule(dynamic)
#endif
 for(int i=0;i<n;++i){
  ConditionalCell cell{beta1[i],beta2[i],phi1[i],phi2[i],B,sd,fixed_ratio[i],nx};
  first[i]=cell.state(aa,aw,nodes,weights).e;
  ConditionalState last{0,0,0};double lastz=INFINITY;
  auto f=[&](double z){cell.rho=std::exp(z);
   if(z!=lastz){last=cell.state(aa,aw,nodes,weights);lastz=z;}
   ConditionalState st=last;
   return std::log(st.e)+std::log(st.m2)-std::log(B-st.e)-std::log(st.m1)-z;};
  double root=std::log(initial[i]),fc=f(root);
  // The supplied equilibrium ratios are already close. Two safeguarded
  // secant/Newton corrections usually avoid a full bracket search.
  if(std::abs(fc)>1e-7){
   double z0=root,f0=fc,z1=z0+.005,f1=f(z1);
   for(int k=0;k<3 && std::abs(f0)>1e-7;++k){
    double slope=(f1-f0)/(z1-z0);
    if(!std::isfinite(slope)||slope>=-.05)break;
    double candidate=z0-f0/slope;
    if(!std::isfinite(candidate)||std::abs(candidate-root)>.2)break;
    double next=f(candidate);z1=z0;f1=f0;z0=candidate;f0=next;
    root=z0;fc=f0;iterations[i]++;
   }
  }
  if(std::abs(fc)>1e-7){
   double lo=root-.15,hi=root+.15,flo=f(lo),fhi=f(hi);
   if(flo*fhi>0){lo=root-2;hi=root+2;flo=f(lo);fhi=f(hi);}
   if(!std::isfinite(flo)||!std::isfinite(fhi)||flo*fhi>0){status[i]=0;continue;}
   int previous=0;
   for(int k=0;k<60;++k){
    double c=(lo*fhi-hi*flo)/(fhi-flo);if(!std::isfinite(c)||c<=lo||c>=hi)c=(lo+hi)/2;
    fc=f(c);root=c;iterations[i]=k+1;
    if(std::abs(fc)<1e-8||hi-lo<1e-10)break;
    if(flo*fc>0){lo=c;flo=fc;if(previous==1)fhi*=.5;previous=1;}
    else{hi=c;fhi=fc;if(previous==-1)flo*=.5;previous=-1;}
   }
  }
  cell.rho=std::exp(root);f(root);ConditionalState st=last;
  equilibrium[i]=st.e;ratios[i]=cell.rho;
  residual[i]=std::log(st.e)+std::log(st.m2)-std::log(B-st.e)-std::log(st.m1)-root;
  status[i]=std::isfinite(residual[i]) && std::abs(residual[i])<2e-6?1:0;
 }
 return List::create(_["E1_first"]=first,_["E1_equilibrium"]=equilibrium,
  _["ratio_equilibrium"]=ratios,_["residual"]=residual,_["status"]=status,_["iterations"]=iterations);
}
