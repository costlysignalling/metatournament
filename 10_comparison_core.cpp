#include <Rcpp.h>
#include <array>
#include <functional>
#include <algorithm>
#include <cmath>
// [[Rcpp::plugins(openmp)]]
#ifdef _OPENMP
#include <omp.h>
#endif
using namespace Rcpp;
struct PriceTraining {
 double b1,b2,p1,p2;
 double w(double x,double a,double b,double p)const{
  return std::expm1(p*std::log1p(a*std::expm1(b*std::log1p(x))));
 }
 double wp(double x,double a,double b,double p)const{
  return p*a*b*std::pow(1+x,b-1)*std::pow(1+a*std::expm1(b*std::log1p(x)),p-1);
 }
 double choose(double a1,double a2,double R,double H,int nx=32)const{
  if(R<=0)return 0;
  if(a1<=0&&a2<=0)return 0;
  if(p1==1&&p2==1&&b1==b2){
   double q=(H+std::log(a1/a2))/(1-b1);
   return std::max(0.,std::min(R,(R+2)/(1+std::exp(-q))-1));
  }
  double c1=std::exp(H/2),c2=std::exp(-H/2);
  auto value=[&](double x){return c1*w(x,a1,b1,p1)+c2*w(R-x,a2,b2,p2);};
  auto derivative=[&](double x){return c1*wp(x,a1,b1,p1)-c2*wp(R-x,a2,b2,p2);};
  if(p1<=1&&p2<=1){
   if(derivative(0)<=0)return 0;if(derivative(R)>=0)return R;
   double lo=0,hi=R;
   for(int j=0;j<32;j++){double mid=(lo+hi)/2;if(derivative(mid)>0)lo=mid;else hi=mid;}
   return (lo+hi)/2;
  }
  std::vector<double> v(nx+1);int best=0;
  for(int j=0;j<=nx;j++){v[j]=value(R*j/nx);if(v[j]>v[best])best=j;}
  double bx=R*best/nx,bv=v[best];std::vector<std::pair<double,double>> brackets;
  for(int j=1;j<nx;j++)if(v[j]>=v[j-1]&&v[j]>=v[j+1])brackets.push_back({R*(j-1)/nx,R*(j+1)/nx});
  if(derivative(0)>0&&derivative(R/nx)<0)brackets.push_back({0,R/nx});
  if(derivative(R)<0&&derivative(R*(nx-1)/nx)>0)brackets.push_back({R*(nx-1)/nx,R});
  for(auto bracket:brackets){
   double lo=bracket.first,hi=bracket.second,r=.618033988749895;
   double x=hi-r*(hi-lo),y=lo+r*(hi-lo),vx=value(x),vy=value(y);
   for(int j=0;j<45;j++)if(vx>vy){hi=y;y=x;vy=vx;x=hi-r*(hi-lo);vx=value(x);}else{lo=x;x=y;vx=vy;y=lo+r*(hi-lo);vy=value(y);}
   if(vx>bv){bx=x;bv=vx;}if(vy>bv){bx=y;bv=vy;}
  }
  return bx;
 }
};
using Eight=std::array<double,8>;
Eight plus8(Eight a,Eight b){for(int t=0;t<8;t++)a[t]+=b[t];return a;}
Eight times8(Eight a,double b){for(double &v:a)v*=b;return a;}
Eight adaptive8(const std::function<Eight(double)> &f,double lo,double hi,double tol,int depth=0){
 static double node[]={.9914553711208126,.9491079123427585,.8648644233597691,.7415311855993945,.5860872354676911,.4058451513773972,.2077849550078985,0};
 static double kw[]={.02293532201052922,.06309209262997855,.1047900103222502,.1406532597155259,.1690047266392679,.1903505780647854,.2044329400752989,.2094821410847278};
 static double gw[]={.1294849661688697,.2797053914892767,.3818300505051189,.4179591836734694};
 double mid=(lo+hi)/2,h=(hi-lo)/2;Eight c=f(mid),k=times8(c,kw[7]),g=times8(c,gw[3]);
 Eight left0={},left1={},right0={},right1={};
 for(int j=0;j<7;j++){
  Eight left=f(mid-h*node[j]),right=f(mid+h*node[j]);
  if(j==0){left0=left;right0=right;}if(j==1){left1=left;right1=right;}
  Eight v=plus8(left,right);k=plus8(k,times8(v,kw[j]));if(j%2)g=plus8(g,times8(v,gw[(j-1)/2]));
 }
 k=times8(k,h);g=times8(g,h);bool good=true;
 // Probe endpoint slivers too. Two nested Gaussian rules can otherwise
 // both miss a narrow strategy region starting at a cohort boundary.
 Eight fl=f(lo+1e-8*(hi-lo)),fr=f(hi-1e-8*(hi-lo));
 double extrap=(1-node[0])/(node[0]-node[1]);
 for(int t=0;t<8;t++){
  double error=std::max(std::abs(k[t]-g[t]),
   std::max(std::abs(fl[t]-(left0[t]+extrap*(left0[t]-left1[t]))),
    std::abs(fr[t]-(right0[t]+extrap*(right0[t]-right1[t]))))*h*(1-node[0]));
  if(error>tol*(std::abs(k[t])+1e-4))good=false;
 }
 if(good||depth>=24)return k;
 return plus8(adaptive8(f,lo,mid,tol,depth+1),adaptive8(f,mid,hi,tol,depth+1));
}
// [[Rcpp::export]]
NumericMatrix comparison_cohort_moments_cpp(NumericVector node,NumericVector weight,
 double R1,double R2,double H,double beta1,double beta2,double phi1,double phi2,
 double slog,double tolerance=1e-5,int nx=32,bool split_outer=false){
 PriceTraining z{beta1,beta2,phi1,phi2};Eight total={};
 auto conditional=[&](double m){
  auto f=[&](double r){
   double a1=std::exp(slog*(m+r)/std::sqrt(2.)),a2=std::exp(slog*(m-r)/std::sqrt(2.));
   double x=z.choose(a1,a2,R1,H,nx),y=z.choose(a2,a1,R2,H,nx);
   Eight v={x,R1-x,z.w(x,a1,beta1,phi1),z.w(R1-x,a2,beta2,phi2),
    y,R2-y,z.w(y,a2,beta1,phi1),z.w(R2-y,a1,beta2,phi2)};
   return times8(v,std::exp(-r*r/2)/std::sqrt(2*M_PI));
  };
  Eight ans={};
  int intervals=std::ceil(std::max(12.,std::min(38.,std::abs(H)/(slog*std::sqrt(2.)*std::min(phi1,phi2))+6))/2);
  for(int k=0;k<intervals;k++)ans=plus8(ans,adaptive8(f,2*k,2*k+2,tolerance));
  return ans;
 };
 if(split_outer&&node.size()){
  // Split the common-talent integral where a strategy boundary crosses
  // the cohort diagonal. This removes the kinks that a whole-line normal
  // rule can turn into artificial campaign-payoff peaks.
  std::vector<double> cuts={-10,-6,-4,-2,0,2,4,6,10};
  for(int group=0;group<2;group++){
   double R=group==0?R1:R2;
   auto allocation=[&](double m){
    double d=group==0?1e-7:-1e-7;
    double a1=std::exp(slog*(m+d)/std::sqrt(2.)),a2=std::exp(slog*(m-d)/std::sqrt(2.));
    return z.choose(a1,a2,R,H,nx)/R;
   };
   for(double target:{1e-8,.5,1-1e-8}){
    double left=-10;bool previous=allocation(left)>=target;
    for(int k=1;k<=40;k++){
     double right=-10+.5*k;bool current=allocation(right)>=target;
     if(current!=previous){
      double lo=left,hi=right;
      for(int it=0;it<26;it++){double mid=(lo+hi)/2;
       if((allocation(mid)>=target)==current)hi=mid;else lo=mid;}
      cuts.push_back((lo+hi)/2);
     }
     left=right;previous=current;
    }
   }
  }
  std::sort(cuts.begin(),cuts.end());
  cuts.erase(std::unique(cuts.begin(),cuts.end(),[](double a,double b){return std::abs(a-b)<1e-7;}),cuts.end());
  std::vector<double> outer_node,outer_weight;
  for(int k=0;k<(int)cuts.size()-1;k++){
   double center=(cuts[k]+cuts[k+1])/2,half=(cuts[k+1]-cuts[k])/2;
   for(int j=0;j<node.size();j++){
    double m=center+half*node[j];
    outer_node.push_back(m);
    outer_weight.push_back(half*weight[j]*std::exp(-m*m/2)/std::sqrt(2*M_PI));
   }
  }
  std::vector<Eight> value(outer_node.size());
  #ifdef _OPENMP
  #pragma omp parallel for num_threads(4) schedule(dynamic)
  #endif
  for(int j=0;j<(int)outer_node.size();j++)
   value[j]=times8(conditional(outer_node[j]),outer_weight[j]);
  // Sum in the original order for reproducible serial/parallel results.
  for(auto v:value)total=plus8(total,v);
 }else if(node.size()){
  std::vector<double> nd=as<std::vector<double>>(node),wt=as<std::vector<double>>(weight);
  std::vector<Eight> value(nd.size());
  #ifdef _OPENMP
  #pragma omp parallel for num_threads(4) schedule(dynamic)
  #endif
  for(int j=0;j<(int)nd.size();j++)value[j]=times8(conditional(nd[j]),wt[j]);
  for(auto v:value)total=plus8(total,v);
 }else{
  auto outer=[&](double m){return times8(conditional(m),std::exp(-m*m/2)/std::sqrt(2*M_PI));};
  const double edge[]={-10,-4,-2,0,2,4,10};
  for(int j=0;j<6;j++)total=plus8(total,adaptive8(outer,edge[j],edge[j+1],tolerance));
 }
 NumericMatrix out(2,4);for(int t=0;t<4;t++){out(0,t)=total[t];out(1,t)=total[t+4];}return out;
}
// [[Rcpp::export]]
NumericVector comparison_member_training_cpp(NumericMatrix talent,double R1,double R2,double H,
 double beta1,double beta2,double phi1,double phi2,int nx=80){
 PriceTraining z{beta1,beta2,phi1,phi2};NumericVector x(talent.nrow());
 for(int i=0;i<talent.nrow();i++)x[i]=z.choose(talent(i,0),talent(i,1),talent(i,0)>talent(i,1)?R1:R2,H,nx);
 return x;
}


