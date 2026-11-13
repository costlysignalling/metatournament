#include <Rcpp.h>
#include <vector>
#include <cmath>
#include <algorithm>
#ifdef _OPENMP
#include <omp.h>
#endif
using namespace Rcpp;
// [[Rcpp::plugins(cpp11)]]
// [[Rcpp::plugins(openmp)]]
struct Cache {
 int np,nx; double B,b1,b2,p1,p2;
 std::vector<double> x,w1,w2,c1,c2;std::vector<int> cap;
 Cache(int n,int g,double budget,double bb1,double bb2,double pp1,double pp2,
   const NumericMatrix& a):np(n),nx(g),B(budget),b1(bb1),b2(bb2),p1(pp1),p2(pp2),
   x(g),w1(n*g),w2(n*g),c1(n*g),c2(n*g),cap(g) {
  for(int k=0;k<nx;++k) x[k]=B*std::pow(double(k)/(nx-1),2);
  for(int k=0;k<nx;++k) cap[k]=std::upper_bound(x.begin(),x.end(),B-x[k]+1e-11)-x.begin()-1;
  for(int q=0;q<np;++q) for(int k=0;k<nx;++k) {
   w1[q*nx+k]=w(x[k],a(q,0),b1,p1);w2[q*nx+k]=w(x[k],a(q,1),b2,p2);
   c1[q*nx+k]=w(B-x[k],a(q,0),b1,p1);c2[q*nx+k]=w(B-x[k],a(q,1),b2,p2);
  }
 }
 static double w(double x,double a,double b,double p) {
  return std::expm1(p*std::log1p(a*std::expm1(b*std::log1p(x))));
 }
};
struct State {double e0,e1,e2,m1,m2;};
double vertex(double x0,double x1,double x2,double y0,double y1,double y2) {
 double s0=(y1-y0)/(x1-x0),s1=(y2-y1)/(x2-x1);
 double aa=(s1-s0)/(x2-x0);
 if(aa>=-1e-15) return x1;
 double v=(x0+x1-s0/aa)/2;
 return std::max(x0,std::min(x2,v));
}
State state(const Cache& c,const NumericMatrix& a,double k1,double k2,double eta) {
 State st{0,0,0,0,0};int nx=c.nx;
 std::vector<double> prefix(nx);std::vector<int> loc(nx);
 for(int q=0;q<c.np;++q) {
  int off=q*nx;
  for(int k=0;k<nx;++k) {
   double v=k1*c.w1[off+k]-eta*c.x[k];
   if(k==0 || v>prefix[k-1]) {prefix[k]=v;loc[k]=k;}
   else {prefix[k]=prefix[k-1];loc[k]=loc[k-1];}
  }
  double best=0,bx=0,by=0,bw1=0,bw2=0;
  auto consider=[&](double x,double y,double w1,double w2) {
   if(x+y>c.B+1e-9)return;
   double val=k1*w1+k2*w2-eta*(x+y);
   if(val>best) {best=val;bx=x;by=y;bw1=w1;bw2=w2;}
  };
  double pv=-1e300,fv1=-1e300,fv2=-1e300;int pi=0,pj=0,fi1=0,fi2=0;
  for(int j=0;j<nx;++j) {
   int i=loc[c.cap[j]];
   double v=prefix[c.cap[j]]+k2*c.w2[off+j]-eta*c.x[j];
   if(v>pv){pv=v;pi=i;pj=j;}
   double v1=k1*c.w1[off+j]+k2*c.c2[off+j]-eta*c.B;
   double v2=k1*c.c1[off+j]+k2*c.w2[off+j]-eta*c.B;
   if(v1>fv1){fv1=v1;fi1=j;}if(v2>fv2){fv2=v2;fi2=j;}
  }
  consider(c.x[pi],c.x[pj],c.w1[off+pi],c.w2[off+pj]);
  auto free_vertex=[&](int i,double scaler,const std::vector<double>& w) {
   if(i==0 || i==nx-1)return c.x[i];
   return vertex(c.x[i-1],c.x[i],c.x[i+1],scaler*w[off+i-1]-eta*c.x[i-1],
    scaler*w[off+i]-eta*c.x[i],scaler*w[off+i+1]-eta*c.x[i+1]);
  };
  double x=free_vertex(pi,k1,c.w1),y=free_vertex(pj,k2,c.w2);
  if(x+y<=c.B) consider(x,y,Cache::w(x,a(q,0),c.b1,c.p1),Cache::w(y,a(q,1),c.b2,c.p2));
  auto full=[&](int i,bool first) {
   double x=c.x[i];
   auto value=[&](int j) {return first ? k1*c.w1[off+j]+k2*c.c2[off+j] : k1*c.c1[off+j]+k2*c.w2[off+j];};
   if(first) consider(x,c.B-x,c.w1[off+i],c.c2[off+i]);
   else consider(c.B-x,x,c.c1[off+i],c.w2[off+i]);
   if(i>0 && i<nx-1) {
    x=vertex(c.x[i-1],c.x[i],c.x[i+1],value(i-1),value(i),value(i+1));
    double xx=first?x:c.B-x,yy=first?c.B-x:x;
    consider(xx,yy,Cache::w(xx,a(q,0),c.b1,c.p1),Cache::w(yy,a(q,1),c.b2,c.p2));
   }
  };
  full(fi1,true);full(fi2,false);
  st.e0+=std::max(0.0,c.B-bx-by);st.e1+=bx;st.e2+=by;st.m1+=bw1;st.m2+=bw2;
 }
 st.e0/=c.np;st.e1/=c.np;st.e2/=c.np;st.m1/=c.np;st.m2/=c.np;return st;
}

// [[Rcpp::export]]
List optout_energy_grid_cpp(NumericVector betas,double phi1,double phi2,
 NumericMatrix population,NumericVector ref1,NumericVector ref2,double eta,double B=50,
 int nx=321,int maxit=180,double tol=0.001,int threads=6,
 Nullable<NumericMatrix> initial1=R_NilValue,Nullable<NumericMatrix> initial2=R_NilValue,
 Nullable<NumericVector> beta2_values=R_NilValue) {
 int n=betas.size();
 NumericVector betas2=beta2_values.isNotNull()?NumericVector(beta2_values):betas;
 if(betas2.size()!=n)stop("beta axes must have the same length");
 bool sym=std::fabs(phi1-phi2)<1e-14 && beta2_values.isNull();
 NumericMatrix f0(n,n),f1(n,n),f2(n,n),r0(n,n),r1(n,n),r2(n,n),err(n,n),s1(n,n),s2(n,n);
 IntegerMatrix status(n,n);
 NumericMatrix in1=initial1.isNotNull()?NumericMatrix(initial1):NumericMatrix(n,n);
 NumericMatrix in2=initial2.isNotNull()?NumericMatrix(initial2):NumericMatrix(n,n);
#ifdef _OPENMP
 omp_set_num_threads(threads);
#pragma omp parallel for schedule(dynamic)
#endif
 for(int ij=0;ij<n*n;++ij) {
  int i=ij/n,j=ij%n;if(sym&&i>j)continue;
  Cache c(population.nrow(),nx,B,betas[i],betas2[j],phi1,phi2,population);
  State first=state(c,population,ref1[i],ref2[j],eta);
  f0(i,j)=first.e0/B;f1(i,j)=first.e1/B;f2(i,j)=first.e2/B;
  double start1=in1(i,j)>0?in1(i,j):ref1[i],start2=in2(i,j)>0?in2(i,j):ref2[j];
  double k1=start1,k2=start2,besterr=1e100,best1=k1,best2=k2;
  State bestst=first;bool ok=false;
  if(eta==1){bestst=State{B,0,0,0,0};besterr=0;best1=0;best2=0;ok=true;}
  const double damping[4]={.4,.2,.08,.03};
  for(int attempt=0;attempt<4 && !ok;++attempt) {
   k1=attempt==0?start1:best1;k2=attempt==0?start2:best2;
   for(int it=0;it<maxit;++it) {
    State st=state(c,population,k1,k2,eta);
    if(st.e1<=0 || st.e2<=0 || st.m1<=0 || st.m2<=0)break;
    double d1=std::log(st.e1/st.m1/k1),d2=std::log(st.e2/st.m2/k2);
    double ee=std::max(std::fabs(d1),std::fabs(d2));
    if(ee<besterr){besterr=ee;bestst=st;best1=k1;best2=k2;}
    if(ee<tol){ok=true;break;}
    k1*=std::exp(damping[attempt]*d1);k2*=std::exp(damping[attempt]*d2);
   }
  }
  r0(i,j)=bestst.e0/B;r1(i,j)=bestst.e1/B;r2(i,j)=bestst.e2/B;
  s1(i,j)=best1;s2(i,j)=best2;err(i,j)=besterr;status(i,j)=ok?1:0;
  if(sym && i<j) {
   f0(j,i)=f0(i,j);f1(j,i)=f2(i,j);f2(j,i)=f1(i,j);
   r0(j,i)=r0(i,j);r1(j,i)=r2(i,j);r2(j,i)=r1(i,j);
   s1(j,i)=best2;s2(j,i)=best1;err(j,i)=besterr;status(j,i)=status(i,j);
  }
 }
 return List::create(_["first"]=List::create(_["s0"]=f0,_["s1"]=f1,_["s2"]=f2),
  _["responsive"]=List::create(_["s0"]=r0,_["s1"]=r1,_["s2"]=r2),
  _["scaler1"]=s1,_["scaler2"]=s2,_["residual"]=err,_["status"]=status);
}
