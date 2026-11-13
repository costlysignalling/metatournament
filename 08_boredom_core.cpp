#include <Rcpp.h>
#include <algorithm>
#include <cmath>
using namespace Rcpp;
struct Local {
 double B,c,b1,b2,p1,p2; int ng;
 std::vector<double> grid;
 Local(double B_,double c_,double b1_,double b2_,double p1_,double p2_,int n):B(B_),c(c_),b1(b1_),b2(b2_),p1(p1_),p2(p2_) {
  for(int j=0;j<=n;j++)grid.push_back(c*B*j/n);
  for(int j=1;j<=n;j++)grid.push_back(c*B+(1-c)*B*j/n);
  ng=grid.size();
 }
 double w(double x,double a,double b,double p){return std::expm1(p*std::log1p(a*std::expm1(b*std::log1p(x))));}
 double w1(double q,double a){return w(q,a,b1,p1);}
 double w2(double q,double a){return w(std::min((1-c)*B,B-q),a,b2,p2);}
 double wp(double x,double a,double b,double p){return p*a*b*std::pow(1+x,b-1)*std::pow(1+a*std::expm1(b*std::log1p(x)),p-1);}
 double val(double q,double a1,double a2,double E1,double E2,double M1,double M2) {
  double u=w1(q,a1),v=w2(q,a2);
  return (M1+u>0?(E1+q)*u/(M1+u):0)+(M2+v>0?(E2+B-q)*v/(M2+v):0);
 }
 std::pair<double,double> br(double a1,double a2,double E1,double E2,double M1,double M2,const std::vector<double>*ww1=nullptr,const std::vector<double>*ww2=nullptr) {
  std::vector<double> vv(ng);int best=0;
  for(int j=0;j<ng;j++) {
   double u=ww1?(*ww1)[j]:w1(grid[j],a1),v=ww2?(*ww2)[j]:w2(grid[j],a2);
   vv[j]=(M1+u>0?(E1+grid[j])*u/(M1+u):0)+(M2+v>0?(E2+B-grid[j])*v/(M2+v):0);
   if(vv[j]>vv[best])best=j;
  }
  double bx=grid[best],bv=vv[best];
  std::vector<std::pair<double,double>> intervals;
  for(int j=1;j<ng-1;j++)if(vv[j]>=vv[j-1]&&vv[j]>=vv[j+1])intervals.push_back({grid[j-1],grid[j+1]});
  auto deriv=[&](double q,bool right){double u=w1(q,a1),v=w2(q,a2),s1=u/(M1+u),s2=v/(M2+v);return s1+(E1+q)/(M1+u)*(1-s1)*wp(q,a1,b1,p1)-s2-(right?(E2+B-q)/(M2+v)*(1-s2)*wp(std::min((1-c)*B,B-q),a2,b2,p2):0.);};
  // A narrow maximum can sit between a grid endpoint and its neighbour.
  if(deriv(0,false)>0&&deriv(grid[1],false)<0)intervals.push_back({grid[0],grid[1]});
  if(deriv(B,true)<0&&deriv(grid[ng-2],true)>0)intervals.push_back({grid[ng-2],grid[ng-1]});
  int mid=(ng-1)/2;
  if(deriv(c*B,false)<0&&deriv(grid[mid-1],false)>0)intervals.push_back({grid[mid-1],grid[mid]});
  if(deriv(c*B,true)>0&&deriv(grid[mid+1],true)<0)intervals.push_back({grid[mid],grid[mid+1]});
  for(auto interval:intervals) {
   double lo=interval.first,hi=interval.second,r=.618033988749895;
   double x=hi-r*(hi-lo),y=lo+r*(hi-lo),vx=val(x,a1,a2,E1,E2,M1,M2),vy=val(y,a1,a2,E1,E2,M1,M2);
   for(int k=0;k<38;k++)if(vx>vy){hi=y;y=x;vy=vx;x=hi-r*(hi-lo);vx=val(x,a1,a2,E1,E2,M1,M2);}else{lo=x;x=y;vx=vy;y=lo+r*(hi-lo);vy=val(y,a1,a2,E1,E2,M1,M2);}
   if(vx>bv){bx=x;bv=vx;}if(vy>bv){bx=y;bv=vy;}
  }
  return {bx,bv};
 }
};
// [[Rcpp::export]]
NumericVector boredom_map_cpp(NumericMatrix actor,NumericVector background,double beta1,double beta2,double phi1,double phi2,double B=50,double compulsory=.5,int grid_half=160) {
 Local f(B,compulsory,beta1,beta2,phi1,phi2,grid_half);NumericVector out(actor.nrow());
 for(int i=0;i<actor.nrow();i++)out[i]=f.br(actor(i,0),actor(i,1),background[0],background[1],background[2],background[3]).first;
 return out;
}
// [[Rcpp::export]]
List boredom_density_cpp(NumericMatrix talent,NumericVector weights,NumericVector moment,double beta1,double beta2,double phi1,double phi2,double B=50,double compulsory=.5,double N=30,int grid_half=80,int maxit=300,double tolerance=1e-7,bool responsive=true) {
 Local f(B,compulsory,beta1,beta2,phi1,phi2,grid_half);int n=talent.nrow();
 std::vector<std::vector<double>> tab1(n),tab2(n);
 for(int i=0;i<n;i++){tab1[i].resize(f.ng);tab2[i].resize(f.ng);for(int j=0;j<f.ng;j++){tab1[i][j]=f.w1(f.grid[j],talent(i,0));tab2[i][j]=f.w2(f.grid[j],talent(i,1));}}
 NumericVector state=clone(moment),best=clone(moment),allocation(n);double error=1e100;int iters=0;bool solved=!responsive;
 auto evaluate=[&](NumericVector m,NumericVector &qq){NumericVector ans(4);for(int i=0;i<n;i++){
  auto out=f.br(talent(i,0),talent(i,1),(N-1)*m[0],(N-1)*(B-m[0]),(N-1)*m[1],(N-1)*m[2],&tab1[i],&tab2[i]);double q=out.first;qq[i]=q;
  ans[0]+=weights[i]*q;ans[1]+=weights[i]*f.w1(q,talent(i,0));ans[2]+=weights[i]*f.w2(q,talent(i,1));ans[3]+=weights[i]*std::max(compulsory*B-q,0.);
 }return ans;};
 if(responsive) for(double damping:{.7,.35,.15,.07}) {
  state=clone(best);
  for(int it=0;it<maxit;it++) {
   NumericVector mapped=evaluate(state,allocation);double r=std::abs(mapped[0]-state[0])/B;
   for(int j=1;j<3;j++)r=std::max(r,std::abs(std::log(std::max(mapped[j],1e-20)/std::max(state[j],1e-20))));
   if(r<error){error=r;best=clone(state);iters=it+1;}
   if(r<tolerance){solved=true;break;}
   state[0]+=damping*(mapped[0]-state[0]);
   for(int j=1;j<3;j++)state[j]=std::exp((1-damping)*std::log(std::max(state[j],1e-20))+damping*std::log(std::max(mapped[j],1e-20)));
  }
  if(solved)break;
 }
 if(responsive&&!solved) {
  // Newton correction in logit(E1/B), log(M1), log(M2), with line search.
  NumericVector th(3);th[0]=std::log(best[0]/(B-best[0]));th[1]=std::log(best[1]);th[2]=std::log(best[2]);
  auto unpack=[&](NumericVector t){return NumericVector::create(B/(1+std::exp(-t[0])),std::exp(t[1]),std::exp(t[2]));};
  auto residual=[&](NumericVector t){auto m=unpack(t);NumericVector qq(n),g=evaluate(m,qq),r(3);r[0]=(g[0]-m[0])/B;r[1]=std::log(g[1]/m[1]);r[2]=std::log(g[2]/m[2]);return r;};
  auto norm=[](NumericVector r){double a=0;for(double x:r)a=std::max(a,std::abs(x));return a;};
  for(int it=0;it<55;it++){
   auto r=residual(th);double er=norm(r);
   if(er<error){error=er;best=unpack(th);}
   if(er<tolerance){solved=true;break;}
   double mat[3][4];
   for(int j=0;j<3;j++){auto t=clone(th);t[j]+=1e-4;auto rr=residual(t);for(int i=0;i<3;i++)mat[i][j]=(rr[i]-r[i])/1e-4;}
   for(int i=0;i<3;i++)mat[i][3]=-r[i];bool valid=true;
   for(int k=0;k<3;k++){
    int pivot=k;for(int i=k+1;i<3;i++)if(std::abs(mat[i][k])>std::abs(mat[pivot][k]))pivot=i;
    if(std::abs(mat[pivot][k])<1e-12){valid=false;break;}
    for(int j=k;j<4;j++)std::swap(mat[k][j],mat[pivot][j]);
    double v=mat[k][k];for(int j=k;j<4;j++)mat[k][j]/=v;
    for(int i=0;i<3;i++)if(i!=k){double a=mat[i][k];for(int j=k;j<4;j++)mat[i][j]-=a*mat[k][j];}
   }
   if(!valid)break;bool accepted=false;
   for(double rate=1.;rate>1e-6;rate*=.5){auto cand=clone(th);for(int j=0;j<3;j++)cand[j]+=rate*std::max(-2.,std::min(2.,mat[j][3]));if(norm(residual(cand))<er){th=cand;accepted=true;break;}}
   if(!accepted)break;
  }
 }
 if(!responsive)best=clone(moment);
 NumericVector final=evaluate(best,allocation);
 return List::create(_["moment"]=best,_["mapped"]=final,_["q"]=allocation,_["converged"]=solved,_["residual"]=responsive?error:0.,_["iterations"]=iters);
}
