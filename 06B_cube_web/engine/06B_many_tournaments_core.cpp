// Allocation engine for the separable T-tournament model. Rcpp is the only
// dependency. Concave cases use the common marginal return. Other cases use
// budget-grid dynamic programming with continuous refinement and GLOBAL
// interpolation / Lagrange-multiplier utility bounds.
#include <Rcpp.h>
#include <algorithm>
#include <cmath>
#include <vector>
#include <functional>
using namespace Rcpp;
// [[Rcpp::plugins(cpp11)]]
struct Activity {
 double a,b,p,c;
 double w(double x) const {return c*std::expm1(p*std::log1p(a*std::expm1(b*std::log1p(x))));}
 double d(double x) const {return c*p*a*b*std::pow(1+x,b-1)*std::pow(1+a*std::expm1(b*std::log1p(x)),p-1);}
 double logcurv(double x) const {return (b-1)/(1+x)+(p-1)*a*b*std::pow(1+x,b-1)/(1+a*std::expm1(b*std::log1p(x)));}
 double k(double x) const {return b-1+a*((b*p-1)*std::pow(1+x,b)-(b-1));}
 double dmax(double B) const {
  double ans=std::max(d(0),d(B));
  if(a>0 && std::abs(b*p-1)>1e-14) {
   double z=(b-1)*(a-1)/(a*(b*p-1));
   if(z>0) {double x=std::pow(z,1/b)-1;if(x>0 && x<B)ans=std::max(ans,d(x));}
  }
  return ans;
 }
 double invert(double lambda,double B,double guess) const {
  if(a==0 || d(0)<=lambda)return 0;
  if(d(B)>=lambda)return B;
  double lo=0,hi=B,x=std::max(1e-12,std::min(B-1e-12,guess));
  for(int j=0;j<48;++j){
   double e=std::log(d(x)/lambda);
   if(std::abs(e)<2e-12)return x;
   if(e>0)lo=x;else hi=x;
   double next=x-e/logcurv(x);
   if(!std::isfinite(next)||next<=lo||next>=hi)next=(lo+hi)/2;
   x=next;
  }
  return x;
 }
};
struct Answer {std::vector<double>x;double value,gap;bool concave;int grid;};

// Improve the grid solution continuously, without weakening its global bound.
// Each pairwise move preserves the budget and is accepted only if utility rises.
void polish_actor(std::vector<Activity> &f, std::vector<double> &x, double B) {
 int T=f.size();
 for(int sweep=0;sweep<160;++sweep) {
  double change=0;
  for(int i=0;i<T;++i)for(int j=i+1;j<T;++j) {
   double di=f[i].d(x[i]),dj=f[j].d(x[j]),grad=di-dj;
   if((x[j]<=0 && grad>=0)||(x[i]<=0 && grad<=0))continue;
   double curvature=di*f[i].logcurv(x[i])+dj*f[j].logcurv(x[j]);
   double step=curvature<0 ? -grad/curvature : (grad>0 ? x[j] : -x[i]);
   step=std::max(-x[i],std::min(x[j],step));
   if(std::abs(step)<2e-13*(1+B))continue;
   double old=f[i].w(x[i])+f[j].w(x[j]);
   for(int back=0;back<30;++back,step*=.5) {
    double value=f[i].w(x[i]+step)+f[j].w(x[j]-step);
    // At floating-point precision allow a Newton correction only in a
    // direction of positive first-order gain and with negative curvature.
    if(value>old || (value>=old-4e-15*(1+std::abs(old)) &&
                    curvature<0 && grad*step>0)) {
     x[i]+=step;x[j]-=step;change=std::max(change,std::abs(step));break;
    }
   }
  }
  if(change<2e-11*(1+B))break;
 }
}
// A second global bound follows from a Lagrange multiplier. Each component's
// second derivative changes sign at most once, so its tilted utility can be
// maximized by checking endpoints and the root on its concave interval.
double tilted_max(const Activity &f,double lambda,double B) {
 double best=std::max(0.0,f.w(B)-lambda*B);
 if(f.a==0)return best;
 double lo=0,hi=B,k0=f.k(0),kB=f.k(B);
 if(k0>=0 && kB>=0)return best;
 if(k0>0 || kB>0) {
  double z=(f.b-1)*(f.a-1)/(f.a*(f.b*f.p-1));
  double turning=std::max(0.0,std::min(B,std::pow(z,1/f.b)-1));
  if(k0>0)lo=turning;else hi=turning;
 }
 if(f.d(lo)>lambda && f.d(hi)<lambda) {
  double x=(lo+hi)/2;
  for(int k=0;k<55;++k) {
   double err=std::log(f.d(x)/lambda);
   if(std::abs(err)<2e-13)break;
   if(err>0)lo=x;else hi=x;
   double next=x-err/f.logcurv(x);
   if(!std::isfinite(next)||next<=lo||next>=hi)next=(lo+hi)/2;
   x=next;
  }
  best=std::max(best,f.w(x)-lambda*x);
 }
 return best;
}

Answer solve_actor(std::vector<Activity> &f,double B,int initial_grid,int max_grid,double tol){
 int T=f.size();Answer out;out.x.assign(T,0);out.value=0;out.gap=0;out.concave=true;out.grid=0;
 if(B==0)return out;
 bool convex=true;int positive=0,last=0;
 for(int t=0;t<T;++t){
  if(f[t].a>0){positive++;last=t;}
  if(std::max(f[t].k(0),f[t].k(B))>1e-12 && f[t].a>0)out.concave=false;
  if(std::min(f[t].k(0),f[t].k(B))< -1e-12 && f[t].a>0)convex=false;
 }
 if(positive==0){out.x.assign(T,B/T);return out;}
 if(positive==1){out.x[last]=B;out.value=f[last].w(B);return out;}
 if(out.concave){
  double low=1e300,high=0;
  for(int t=0;t<T;++t)if(f[t].a>0){low=std::min(low,f[t].d(B));high=std::max(high,f[t].d(0));}
  low=std::log(low);high=std::log(high);out.x.assign(T,B/T);
  for(int j=0;j<52;++j){
   double lm=(low+high)/2,lambda=std::exp(lm),sum=0;
   for(int t=0;t<T;++t){out.x[t]=f[t].invert(lambda,B,out.x[t]);sum+=out.x[t];}
   if(std::abs(sum-B)<1e-10*(1+B))break;
   if(sum>B)low=lm;else high=lm;
  }
  double sum=0;int largest=0;
  for(int t=0;t<T;++t){sum+=out.x[t];if(out.x[t]>out.x[largest])largest=t;}
  // A linear component has a set-valued inverse at its constant marginal
  // return. Let such components absorb the remaining budget at that cutoff.
  double correction=B-sum,cutoff=std::exp((low+high)/2);
  for(int t=0;t<T && std::abs(correction)>1e-11*(1+B);++t) {
   if(f[t].a>0 && std::max(std::abs(f[t].k(0)),std::abs(f[t].k(B)))<1e-11 &&
      std::abs(std::log(f[t].d(0)/cutoff))<1e-8) {
    double amount=std::max(-out.x[t],std::min(B-out.x[t],correction));
    out.x[t]+=amount;correction-=amount;
   }
  }
  largest=std::max_element(out.x.begin(),out.x.end())-out.x.begin();
  out.x[largest]+=correction;
  double lambda=f[largest].d(out.x[largest]);
  for(int t=0;t<T;++t){
   out.value+=f[t].w(out.x[t]);
   // Supporting tangents give a global upper bound for concave utilities.
   out.gap+=std::max(0.0,f[t].d(out.x[t])-lambda)*(B-out.x[t]);
   out.gap+=std::max(0.0,lambda-f[t].d(out.x[t]))*out.x[t];
  }
  return out;
 }
 if(convex){
  int best=0;for(int t=1;t<T;++t)if(f[t].w(B)>f[best].w(B))best=t;
  out.x[best]=B;out.value=f[best].w(B);return out;
 }
 // On a grid with B=N*h, the maximum of the piecewise-linear interpolants
 // is attained at a grid allocation. Each interpolation error is at most
 // h^2*sup|f''|/8. This yields a global upper bound, not a local-optimum test.
 double curvature=0;
 for(int t=0;t<T;++t)curvature+=f[t].dmax(B)*(1-f[t].b+std::abs(f[t].p-1)*f[t].a*f[t].b);
 // Process nonconcave components first. Convolution with a concave
 // component has monotone maximizing indices and admits exact divide-and-
 // conquer on the same grid; no approximation is introduced by this shortcut.
 std::vector<int> order(T);std::vector<bool> component_concave(T);
 for(int t=0;t<T;++t){order[t]=t;component_concave[t]=f[t].a==0 || std::max(f[t].k(0),f[t].k(B))<=0;}
 std::stable_sort(order.begin(),order.end(),[&](int i,int j){return component_concave[i]<component_concave[j];});
 std::vector<double> best_x;double best_value=-1e300;
 for(int N=initial_grid;;N=std::min(2*N,max_grid)){
  double h=B/N;
  std::vector<double>previous(N+1,-1e300),next(N+1),pay(N+1);
  std::vector<std::vector<int> > choice(T,std::vector<int>(N+1,0));
  previous[0]=0;
  for(int stage=0;stage<T;++stage){
   int t=order[stage];
   for(int k=0;k<=N;++k)pay[k]=f[t].w(h*k);
   if(stage==0) {
    for(int n=0;n<=N;++n){next[n]=pay[n];choice[t][n]=n;}
   } else if(component_concave[t]) {
    std::function<void(int,int,int,int)> fill=[&](int left,int right,int opt_left,int opt_right){
     if(left>right)return;
     int mid=(left+right)/2,limit=std::min(mid,opt_right),arg=opt_left;
     double best=-1e300;
     for(int j=opt_left;j<=limit;++j){double v=previous[j]+pay[mid-j];if(v>best){best=v;arg=j;}}
     next[mid]=best;choice[t][mid]=mid-arg;
     fill(left,mid-1,opt_left,arg);fill(mid+1,right,arg,opt_right);
    };
    fill(0,N,0,N);
   } else {
    for(int n=0;n<=N;++n){
     double best=-1e300;int arg=0;
     for(int k=0;k<=n;++k){double v=previous[n-k]+pay[k];if(v>best){best=v;arg=k;}}
     next[n]=best;choice[t][n]=arg;
    }
   }
   previous.swap(next);
  }
  double grid_value=previous[N];out.grid=N;
  int n=N;for(int stage=T-1;stage>=0;--stage){int t=order[stage],k=choice[t][n];out.x[t]=h*k;n-=k;}
  polish_actor(f,out.x,B);
  // Near a switch, the best grid point may belong to a different basin than
  // the best continuous allocation. Compare continuously optimized specialist
  // seeds too, so grid rounding does not move the switching boundary.
  if(N==initial_grid) {
   double best=0;for(int t=0;t<T;++t)best+=f[t].w(out.x[t]);
   for(int seed=0;seed<T;++seed) {
    std::vector<double> trial(T,0);trial[seed]=B;
    polish_actor(f,trial,B);
    double value=0;for(int t=0;t<T;++t)value+=f[t].w(trial[t]);
    if(value>best){out.x=trial;best=value;}
   }
  }
  out.value=0;int largest=0;
  for(int t=0;t<T;++t){out.value+=f[t].w(out.x[t]);if(out.x[t]>out.x[largest])largest=t;}
  if(out.value>best_value){best_value=out.value;best_x=out.x;}
  else {out.x=best_x;out.value=best_value;largest=std::max_element(out.x.begin(),out.x.end())-out.x.begin();}
  out.gap=std::max(0.0,grid_value+curvature*h*h/8-out.value);
  double lambda=f[largest].d(out.x[largest]),dual=lambda*B-out.value;
  for(int t=0;t<T;++t)dual+=tilted_max(f[t],lambda,B);
  out.gap=std::min(out.gap,std::max(0.0,dual)+1e-13*(1+std::abs(out.value)));
  if(out.gap<=tol*(1+std::abs(out.value)) || N==max_grid)break;
 }
 return out;
}
// [[Rcpp::export]]
List many_allocations_cpp(NumericMatrix talents,NumericVector budgets,
 NumericVector beta,NumericVector phi,NumericVector scaler,
 int initial_grid=128,int max_grid=2048,double utility_tolerance=1e-4){
 int n=talents.nrow(),T=talents.ncol();NumericMatrix X(n,T);
 NumericVector U(n),gap(n);LogicalVector concave(n);IntegerVector grid(n);
 std::vector<Activity> f(T);
 for(int i=0;i<n;++i){
  if(i%256==0)checkUserInterrupt();
  for(int t=0;t<T;++t)f[t]={talents(i,t),beta[t],phi[t],scaler[t]};
  Answer a=solve_actor(f,budgets[i],initial_grid,max_grid,utility_tolerance);
  for(int t=0;t<T;++t)X(i,t)=a.x[t];U[i]=a.value;gap[i]=a.gap;concave[i]=a.concave;grid[i]=a.grid;
 }
 return List::create(_["allocation"]=X,_["utility"]=U,_["utility_gap_bound"]=gap,
  _["concave"]=concave,_["budget_grid"]=grid);
}
