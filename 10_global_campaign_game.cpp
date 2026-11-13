#include <Rcpp.h>
#include <algorithm>
#include <cmath>
#include <vector>
#include <functional>
using namespace Rcpp;
using Vec=std::vector<double>;
struct Broadcast {
 int n,nx; double B,eps,b1,b2,p1,p2;
 Vec a1,a2;
 double w(double x,double a,double b,double p)const{
  return std::expm1(p*std::log1p(a*std::expm1(b*std::log1p(std::max(0.,x)))));
 }
 template<class F> std::pair<double,double> maximum(F value,double R,int grid,int refine=25)const{
  Vec xx,yy;
  for(int j=0;j<=grid;j++)xx.push_back(R*j/grid);
  if(R>0)for(double z:{1e-7,1e-5,.001,.01}){xx.push_back(R*z);xx.push_back(R*(1-z));}
  std::sort(xx.begin(),xx.end());xx.erase(std::unique(xx.begin(),xx.end()),xx.end());
  int k=0;for(int j=0;j<(int)xx.size();j++){yy.push_back(value(xx[j]));if(yy[j]>yy[k])k=j;}
  double bx=xx[k],bv=yy[k];
  for(int j=1;j<(int)xx.size()-1;j++)if(yy[j]>=yy[j-1]&&yy[j]>=yy[j+1]){
   double lo=xx[j-1],hi=xx[j+1],r=.618033988749895;
   double x=hi-r*(hi-lo),y=lo+r*(hi-lo),fx=value(x),fy=value(y);
   for(int it=0;it<refine;it++)if(fx>fy){hi=y;y=x;fy=fx;x=hi-r*(hi-lo);fx=value(x);}else{
    lo=x;x=y;fx=fy;y=lo+r*(hi-lo);fy=value(y);}
   if(fx>bv){bx=x;bv=fx;}if(fy>bv){bx=y;bv=fy;}
  }
  return {bx,bv};
 }
 struct State {Vec x,v,h,u,material; double regret;int iterations;bool converged;};
 State training(const Vec &v,const Vec &initial,int maxit=350,double tol=2e-5)const{
  Vec R(n),x(n),h(n),c1(n),c2(n),w1(n),w2(n);double sv=0,C1=0,C2=0;
  for(double z:v){sv+=z;C1+=std::max(z,0.);C2+=std::max(-z,0.);}
  for(int i=0;i<n;i++){
   R[i]=B-std::abs(v[i]);x[i]=std::max(0.,std::min(R[i],initial.size()==(size_t)n?initial[i]:R[i]/2));
   h[i]=eps*(sv-v[i]);c1[i]=std::exp(std::min(h[i],0.));c2[i]=std::exp(std::min(-h[i],0.));
  }
  bool symmetry=b1==b2&&p1==p2&&n%2==0;
  if(symmetry)for(int i=0;i<n/2;i++)if(std::abs(v[2*i]+v[2*i+1])>1e-10) symmetry=false;
  if(symmetry)for(int i=0;i<n/2;i++)if(a1[2*i]!=a2[2*i+1]||a2[2*i]!=a1[2*i+1])symmetry=false;
  double regret=1e100;int iteration=0;Vec u(n),mat(n),best(n);
  for(int it=0;it<maxit;it++){
   double E1=C1,E2=C2,M1=0,M2=0;
   for(int i=0;i<n;i++){w1[i]=w(x[i],a1[i],b1,p1);w2[i]=w(R[i]-x[i],a2[i],b2,p2);
    E1+=x[i];E2+=R[i]-x[i];M1+=w1[i];M2+=w2[i];}
   regret=0;
   for(int i=0;i<n;i++){
    double e1=E1-x[i],e2=E2-(R[i]-x[i]),m1=std::max(0.,M1-w1[i]),m2=std::max(0.,M2-w2[i]);
    auto val=[&](double q){double z1=w(q,a1[i],b1,p1),z2=w(R[i]-q,a2[i],b2,p2);
     return c1[i]*(m1+z1>0?(e1+q)*z1/(m1+z1):0)+
      c2[i]*(m2+z2>0?(e2+R[i]-q)*z2/(m2+z2):0);};
    auto opt=maximum(val,R[i],nx);best[i]=opt.first;u[i]=val(x[i]);
    mat[i]=(M1>0?E1*w1[i]/M1:0)+(M2>0?E2*w2[i]/M2:0);
    double gain=opt.second-u[i];regret=std::max(regret,gain);
    if(gain>tol){
     double rate=it<8?.7:1.;double next=x[i]+rate*(best[i]-x[i]);
     double nw1=w(next,a1[i],b1,p1),nw2=w(R[i]-next,a2[i],b2,p2);
     E1+=next-x[i];E2-=next-x[i];M1+=nw1-w1[i];M2+=nw2-w2[i];
     x[i]=next;w1[i]=nw1;w2[i]=nw2;
    }
   }
   iteration=it+1;
   if(regret<tol)return {x,v,h,u,mat,regret,iteration,true};
   if(symmetry)for(int i=0;i<n/2;i++){double q=(x[2*i]+R[2*i]-x[2*i+1])/2;
    x[2*i]=q;x[2*i+1]=R[2*i]-q;}
  }
  return {x,v,h,u,mat,regret,iteration,false};
 }
 struct Choice {double v,value,gain;State state;int failed;double failed_regret;};
 Choice campaign(int i,const State &base,int grid=16,int refine=12,int direction=0)const{
  int failed=0;double failed_regret=0;State last=base;double h=base.h[i];
  auto action=[&](double z){return direction==0?z-B:direction*z;};
  auto value=[&](double z){Vec v=base.v;v[i]=action(z);
   auto q=training(v,base.x);if(!q.converged){failed++;failed_regret=std::max(failed_regret,q.regret);return -1e100;}
   last=q;return q.u[i];};
  auto best=maximum(value,(direction==0?2*B:B)-2e-6,grid,refine);
  // Explicitly retain the current policy and no campaign.
  double zero=value(direction==0?B:0),current=base.u[i];double chosen=action(best.first),bv=best.second;
  if(zero>bv){chosen=0;bv=zero;}
  if(current>=bv){chosen=base.v[i];bv=current;}
  Vec v=base.v;v[i]=chosen;last=training(v,base.x);
  return {chosen,bv,bv-current,last,failed,failed_regret};
 }
};
Broadcast make_broadcast(NumericMatrix a,NumericVector beta,NumericVector phi,double B,double eps,int nx){
 Broadcast f;f.n=a.nrow();f.B=B;f.eps=eps;f.b1=beta[0];f.b2=beta[1];f.p1=phi[0];f.p2=phi[1];f.nx=nx;
 for(int i=0;i<f.n;i++){f.a1.push_back(a(i,0));f.a2.push_back(a(i,1));}return f;
}
List export_broadcast(const Broadcast::State &s){return List::create(_["x1"]=s.x,_["campaign"]=s.v,
 _["received_pressure"]=s.h,_["perceived_scaled"]=s.u,_["material"]=s.material,
 _["training_regret"]=s.regret,_["training_converged"]=s.converged,_["iterations"]=s.iterations);}
// [[Rcpp::export]]
List broadcast_training_cpp(NumericMatrix talent,NumericVector beta,NumericVector phi,
 NumericVector campaign,double B=50,double efficiency=.1,int nx=32,NumericVector initial=NumericVector()){
 auto f=make_broadcast(talent,beta,phi,B,efficiency,nx);
 return export_broadcast(f.training(as<Vec>(campaign),as<Vec>(initial)));
}
// [[Rcpp::export]]
NumericVector broadcast_gradient_cpp(NumericMatrix talent,NumericVector beta,NumericVector phi,
 NumericVector campaign,double B=50,double efficiency=.1,int nx=24,double step=.03){
 auto f=make_broadcast(talent,beta,phi,B,efficiency,nx);Vec v=as<Vec>(campaign);
 auto base=f.training(v,Vec());NumericVector g(f.n);
 if(!base.converged){std::fill(g.begin(),g.end(),NA_REAL);return g;}
 for(int i=0;i<f.n;i++){
  Vec vp=v,vm=v;vp[i]+=step;vm[i]-=step;
  if(vp[i]>=B||vm[i]<=-B){g[i]=NA_REAL;continue;}
  auto plus=f.training(vp,base.x),minus=f.training(vm,base.x);
  if(!plus.converged||!minus.converged){g[i]=NA_REAL;continue;}
  if(std::abs(v[i])>step)g[i]=(plus.u[i]-minus.u[i])/(2*step);
  else{
   double gp=(plus.u[i]-base.u[i])/step,gm=(base.u[i]-minus.u[i])/step;
   g[i]=gp>0?gp:(gm<0?gm:0.);
  }
 }
 return g;
}
// [[Rcpp::export]]
List broadcast_campaign_cpp(NumericMatrix talent,NumericVector beta,NumericVector phi,
 NumericVector campaign,double B=50,double efficiency=.1,int nx=24,int grid=12,int maxit=12,
 double tolerance=.005,bool sequential=true,bool enforce_symmetry=false,int direction=0,double relaxation=.5){
 auto f=make_broadcast(talent,beta,phi,B,efficiency,nx);Vec v=as<Vec>(campaign),gains(f.n),target(f.n);
 auto state=f.training(v,Vec());int failures=0;bool converged=false;
 bool symmetry=enforce_symmetry&&beta[0]==beta[1]&&phi[0]==phi[1]&&f.n%2==0;
 if(symmetry)for(int i=0;i<f.n/2;i++)if(f.a1[2*i]!=f.a2[2*i+1]||f.a2[2*i]!=f.a1[2*i+1])symmetry=false;
 int final_failures=0;double failed_regret=0;
 for(int it=0;it<maxit;it++){
  if(!state.converged)break;
  double maximum_gain=0;final_failures=0;
  for(int i=0;i<f.n;i++){
   if(symmetry&&i%2){target[i]=-target[i-1];gains[i]=gains[i-1];continue;}
   auto br=f.campaign(i,state,grid,12,direction);target[i]=br.v;gains[i]=br.gain;
   failures+=br.failed;final_failures+=br.failed;maximum_gain=std::max(maximum_gain,br.gain);
   failed_regret=std::max(failed_regret,br.failed_regret);
   if(sequential&&!symmetry&&br.gain>tolerance){
    v[i]+=relaxation*(br.v-v[i]);state=f.training(v,state.x);
    if(!state.converged){final_failures++;break;}
   }
  }
  Rcout << "Campaign iteration " << it+1 << "; max gain " << maximum_gain << "; failed training profiles " << failures << std::endl;
  if(maximum_gain<tolerance){converged=final_failures==0;break;}
  if(!sequential||symmetry){double rate=it<4?.3:.12;
   for(int i=0;i<f.n;i++)v[i]+=rate*(target[i]-v[i]);state=f.training(v,state.x);}
  Rcpp::checkUserInterrupt();
 }
 auto out=export_broadcast(state);out["campaign_gain"]=gains;out["campaign_converged"]=converged;
 out["failed_training_profiles"]=final_failures;out["total_failed_training_profiles"]=failures;
 out["max_failed_training_regret"]=failed_regret;return out;
}
