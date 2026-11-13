#include <Rcpp.h>
#include <algorithm>
#include <cmath>
#include <vector>
using namespace Rcpp;
struct CampaignLinear {
 NumericVector c;NumericMatrix moments;double B,b1,b2;std::vector<double> y[2],s[2];
 void setup(){int n=c.size();double h=c[1]-c[0];
  for(int t=0;t<2;t++){
   y[t].resize(n);s[t].resize(n);std::vector<double>d(n-1);
   for(int j=0;j<n;j++)y[t][j]=std::log(moments(j,t==0?0:3)/moments(j,t==0?1:2));
   for(int j=0;j<n-1;j++)d[j]=(y[t][j+1]-y[t][j])/h;
   for(int j=1;j<n-1;j++)s[t][j]=d[j-1]*d[j]>0?2*d[j-1]*d[j]/(d[j-1]+d[j]):0;
   auto edge=[](double a,double b){double m=(3*a-b)/2;if(m*a<=0)return 0.;return a*b<0&&std::abs(m)>3*std::abs(a)?3*a:m;};
   s[t][0]=edge(d[0],d[1]);s[t][n-1]=edge(d[n-2],d[n-3]);
  }
 }
 double logL(double z,int t)const{int n=c.size();double h=c[1]-c[0],u=(z-c[0])/h;
  int j=std::max(0,std::min(n-2,(int)std::floor(u)));double f=std::max(0.,std::min(1.,u-j)),f2=f*f,f3=f2*f;
  return (2*f3-3*f2+1)*y[t][j]+(f3-2*f2+f)*h*s[t][j]+(-2*f3+3*f2)*y[t][j+1]+(f3-f2)*h*s[t][j+1];
 }
 std::pair<double,double> training(double z,double a1,double a2)const{
  double R=B-std::abs(z);if(R<=0)return {0,0};double l1=logL(z,0),l2=logL(z,1),q=std::log(a1/a2)+l1-l2,x;
  auto der=[&](double xx){return q+std::log(b1/b2)+(b1-1)*std::log1p(xx)-(b2-1)*std::log1p(R-xx);};
  if(der(0)<=0)x=0;else if(der(R)>=0)x=R;else if(b1==b2){double p=1/(1+std::exp(-q/(1-b1)));x=std::max(0.,std::min(R,(R+2)*p-1));}
  else{double lo=0,hi=R;for(int j=0;j<38;j++){double mid=(lo+hi)/2;if(der(mid)>0)lo=mid;else hi=mid;}x=(lo+hi)/2;}
  double v=std::exp(l1)*a1*std::expm1(b1*std::log1p(x))+std::exp(l2)*a2*std::expm1(b2*std::log1p(R-x));
  return {x,v};
 }
};
// [[Rcpp::export]]
NumericMatrix campaign_linear_focal_cpp(NumericMatrix talent,NumericMatrix moments,NumericVector effort,
 NumericVector beta,double B=50){
 CampaignLinear model{effort,moments,B,beta[0],beta[1]};model.setup();NumericMatrix out(talent.nrow(),7);
 for(int i=0;i<talent.nrow();i++){
  double a1=talent(i,0),a2=talent(i,1);auto base=model.training(0,a1,a2),best=base;double chosen=0;
  auto evaluate=[&](double c){auto z=model.training(c,a1,a2);if(z.second>best.second+1e-11*std::max(1.,best.second)){best=z;chosen=c;}return z.second;};
  std::vector<double> v;for(double c:effort)v.push_back(evaluate(c));
  for(int k=1;k<effort.size()-1;k++)if(v[k]>=v[k-1]&&v[k]>=v[k+1]){
   double lo=effort[k-1],hi=effort[k+1],r=.618033988749895,x=hi-r*(hi-lo),y=lo+r*(hi-lo),fx=evaluate(x),fy=evaluate(y);
   for(int it=0;it<38;it++)if(fx>fy){hi=y;y=x;fy=fx;x=hi-r*(hi-lo);fx=evaluate(x);}else{lo=x;x=y;fx=fy;y=lo+r*(hi-lo);fy=evaluate(y);}
  }
  if(best.second-base.second<=1e-8*(1+base.second)){best=base;chosen=0;}
  out(i,0)=best.first;out(i,1)=B-std::abs(chosen)-best.first;out(i,2)=chosen;out(i,3)=best.second;
  out(i,4)=base.first;out(i,5)=base.second;out(i,6)=best.second-base.second;
 }
 return out;
}
