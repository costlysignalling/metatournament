#include <Rcpp.h>
#include <cmath>
#include <algorithm>
using namespace Rcpp;
double hw(double z,double a,double b,double p){return std::expm1(p*std::log1p(a*std::expm1(b*std::log1p(z))));}
double hd(double z,double a,double b,double p){return p*a*b*std::pow(1+z,b-1)*std::pow(1+a*std::expm1(b*std::log1p(z)),p-1);}
double normcdf(double z){return .5*std::erfc(-z/std::sqrt(2.));}
double normtail(double z){return .5*std::erfc(z/std::sqrt(2.));}
struct HEnvelope{
 double a,b,p,B,lo,hi,klo,khi;bool chord,linear;
 HEnvelope(double a_,double b_,double p_,double B_):a(a_),b(b_),p(p_),B(B_),lo(0),hi(B_),chord(false),linear(false){
  double A=(b-1)*(1-a),C=a*(b*p-1),c0=A+C,cB=A+C*std::pow(1+B,b),tol=1e-12*(1+std::abs(A)+std::abs(C)*std::pow(1+B,b));
  linear=std::max(std::abs(c0),std::abs(cB))<tol;
  if(linear || std::min(c0,cB)>=-tol){chord=true;lo=hi=0;klo=khi=hw(B,a,b,p)/B;return;}
  if(std::max(c0,cB)<=tol){khi=hd(0,a,b,p);klo=hd(B,a,b,p);return;}
  double inf=std::pow(-A/C,1/b)-1;
  if(c0>0){
   auto f=[&](double z){return z*hd(z,a,b,p)-hw(z,a,b,p);};
   if(f(B)>=0){chord=true;lo=hi=0;klo=khi=hw(B,a,b,p)/B;return;}
   double left=std::max(0.,inf),right=B;
   for(int k=0;k<48;k++){double m=(left+right)/2;if(f(m)>0)left=m;else right=m;}lo=(left+right)/2;
  }else{
   auto f=[&](double z){return (B-z)*hd(z,a,b,p)-(hw(B,a,b,p)-hw(z,a,b,p));};
   if(f(0)<=0){chord=true;lo=hi=0;klo=khi=hw(B,a,b,p)/B;return;}
   double left=0,right=std::min(B,inf);
   for(int k=0;k<48;k++){double m=(left+right)/2;if(f(m)>0)left=m;else right=m;}hi=(left+right)/2;
  }
  khi=lo>1e-10?hw(lo,a,b,p)/lo:hd(0,a,b,p);
  klo=hi<B-1e-10?(hw(B,a,b,p)-hw(hi,a,b,p))/(B-hi):hd(B,a,b,p);
 }
 double alloc(double k)const{
  if(k>=khi)return 0;if(k<=klo)return B;
  double left=lo,right=hi;
  for(int j=0;j<36;j++){double m=(left+right)/2;if(hd(m,a,b,p)>k)left=m;else right=m;}
  return (left+right)/2;
 }
 double slope(double z)const{return chord?khi:(z<lo?khi:(z>hi?klo:hd(z,a,b,p)));}
 bool contact(double z)const{return linear||(!chord&&z>=lo-1e-9&&z<=hi+1e-9);}
};
struct HEngine{
 double b,p1,p2,B,s,EA;bool both;
 NumericVector ax,aw,gx,gw;std::vector<HEnvelope> env;
 std::vector<double> cum,tail,der;int nt;double step;
 HEngine(double b_,double p1_,double p2_,double B_,NumericVector ax_,NumericVector aw_,NumericVector gx_,NumericVector gw_,bool both_):b(b_),p1(p1_),p2(p2_),B(B_),both(both_),ax(ax_),aw(aw_),gx(gx_),gw(gw_){
  s=std::sqrt(std::log((1+std::sqrt(2.))/2.));nt=4097;step=40./(nt-1);cum.resize(nt);tail.resize(nt);der.resize(nt);std::vector<double> increments(nt-1);
  auto f=[&](double z){return std::expm1(p1*std::log1p(std::exp(s*z)))*std::exp(-z*z/2)/std::sqrt(2*M_PI);};
  for(int i=0;i<nt;i++)der[i]=f(-20+i*step);
  cum[0]=0;for(int i=0;i<nt-1;i++){increments[i]=step/6*(der[i]+4*f(-20+(i+.5)*step)+der[i+1]);cum[i+1]=cum[i]+increments[i];}
  EA=cum[nt-1];tail[nt-1]=0;for(int i=nt-2;i>=0;i--)tail[i]=tail[i+1]+increments[i];
  for(double a:ax)env.emplace_back(a,b,p2,B);
 }
 double inv(double v)const{return std::expm1(std::log1p(std::max(0.,v))/p1);}
 double W1(double a)const{return std::expm1(p1*std::log1p(a));}
 double tm(double a)const{
  if(a<=0)return EA;double z=std::log(a)/s;
  if(z<=-20)return EA;if(z>=20)return 0;
  int i=std::max(0,std::min(nt-2,int((z+20)/step)));double u=(z+20-i*step)/step;
  // Cubic Hermite interpolation of the directly accumulated upper moment.
  return std::max(0.,(2*u*u*u-3*u*u+1)*tail[i]+(u*u*u-2*u*u+u)*step*(-der[i])+(-2*u*u*u+3*u*u)*tail[i+1]+(u*u*u-u*u)*step*(-der[i+1]));
 }
 NumericVector moments(double ratio)const{
  NumericVector ans(4);
  for(int j=0;j<ax.size();j++){
   double A2=std::expm1(p2*std::log1p(ax[j])),cut=inv(A2/ratio);
   if(both){double z=std::log(std::max(cut,1e-300))/s;ans[0]+=aw[j]*B*normtail(z);ans[1]+=aw[j]*B*normcdf(z);ans[2]+=aw[j]*B*tm(cut);ans[3]+=aw[j]*B*A2*normcdf(z);continue;}
   const auto &e=env[j];double al=inv(e.klo/ratio),ah=inv(e.khi/ratio);
   double zl=std::log(std::max(al,1e-300))/s,zh=std::log(std::max(ah,1e-300))/s;
   double pl=normcdf(zl),ph=normtail(zh),mass=std::max(0.,1-pl-ph);
   double E1=B*ph,E2=B*pl,M1=B*tm(ah),M2=hw(B,ax[j],b,p2)*pl;
   if(mass>1e-15 && zh-zl>1e-10){
    double l=std::max(-12.,zl),h=std::min(12.,zh),normalization=0;
    int chunks=std::max(1,int(std::ceil((h-l)/2.)));int nk=chunks*gx.size();
    std::vector<double> ws(nk),qs(nk),As(nk),w2(nk);
    for(int k=0;k<nk;k++){
     int block=k/gx.size(),node=k%gx.size();double width=(h-l)/chunks;
     double z=l+width*(block+.5)+width*gx[node]/2,a=std::exp(s*z);ws[k]=gw[node]*std::exp(-z*z/2);normalization+=ws[k];
     As[k]=W1(a);qs[k]=e.alloc(ratio*As[k]);w2[k]=hw(qs[k],ax[j],b,p2);
    }
    if(normalization>0)for(int k=0;k<nk;k++){double wt=mass*ws[k]/normalization;E1+=wt*(B-qs[k]);E2+=wt*qs[k];M1+=wt*(B-qs[k])*As[k];M2+=wt*w2[k];}
   }
   ans[0]+=aw[j]*E1;ans[1]+=aw[j]*E2;ans[2]+=aw[j]*M1;ans[3]+=aw[j]*M2;
  }return ans;
 }
};
// [[Rcpp::export]]
SEXP harvesting_engine_cpp(double beta2,double phi1,double phi2,double B,NumericVector alpha2,NumericVector alpha_weights,NumericVector gl_nodes,NumericVector gl_weights,bool both_fixed=false){
 XPtr<HEngine> ptr(new HEngine(beta2,phi1,phi2,B,alpha2,alpha_weights,gl_nodes,gl_weights,both_fixed),true);return ptr;
}
// [[Rcpp::export]]
NumericVector harvesting_moments_cpp(SEXP ptr,double ratio){XPtr<HEngine> e(ptr);return e->moments(ratio);}
// [[Rcpp::export]]
NumericVector harvesting_reference_cpp(SEXP ptr){
 XPtr<HEngine> e(ptr);double m2=0;
 if(e->both)for(int j=0;j<e->ax.size();j++)m2+=e->aw[j]*e->B/2*std::expm1(e->p2*std::log1p(e->ax[j]));
 else for(int j=0;j<e->ax.size();j++)for(int k=0;k<e->gx.size();k++)m2+=e->aw[j]*e->gw[k]/2*hw(e->B*(e->gx[k]+1)/2,e->ax[j],e->b,e->p2);
 return NumericVector::create(e->B/2,e->B/2,e->B/2*e->EA,m2);
}
// [[Rcpp::export]]
NumericVector harvesting_map_cpp(SEXP ptr,double ratio,NumericMatrix actor){
 XPtr<HEngine> e(ptr);NumericVector out(actor.nrow());
 for(int i=0;i<actor.nrow();i++){
  double k=ratio*e->W1(actor(i,0));
  if(e->both){double a=std::expm1(e->p2*std::log1p(actor(i,1)));out[i]=k>a?e->B:(k<a?0:e->B/2);}
  else {HEnvelope env(actor(i,1),e->b,e->p2,e->B);out[i]=e->B-env.alloc(k);}
 }return out;
}
// [[Rcpp::export]]
NumericMatrix harvesting_boundaries_cpp(SEXP ptr,double ratio,NumericVector alpha2){
 XPtr<HEngine> e(ptr);NumericMatrix out(alpha2.size(),4);
 for(int i=0;i<alpha2.size();i++){
  if(e->both){double c=e->inv(std::expm1(e->p2*std::log1p(alpha2[i]))/ratio);for(int k=0;k<3;k++)out(i,k)=c;out(i,3)=1;}
  else{HEnvelope env(alpha2[i],e->b,e->p2,e->B);out(i,0)=e->inv(env.klo/ratio);out(i,1)=e->inv(env.khi/ratio);out(i,2)=e->inv(hd(e->B/2,alpha2[i],e->b,e->p2)/ratio);out(i,3)=env.contact(e->B/2);}
 }return out;
}
// [[Rcpp::export]]
NumericVector harvesting_presence_cpp(SEXP ptr,double ratio,double radius,int n=193){
 XPtr<HEngine> e(ptr);NumericVector score(4,-100.);
 for(int i=0;i<n;i++){
  double th=-M_PI/2+M_PI*i/(n-1),y=std::exp(-e->s*e->s+e->s*radius*std::sin(th));
  double low=-e->s*e->s-e->s*radius*std::cos(th),high=-e->s*e->s+e->s*radius*std::cos(th);
  double cL,cF,cH,ready1=-1,ready2=-1;
  if(e->both){cL=cF=cH=e->inv(std::expm1(e->p2*std::log1p(y))/ratio);}
  else{HEnvelope en(y,e->b,e->p2,e->B);cL=e->inv(en.klo/ratio);cF=e->inv(en.khi/ratio);cH=e->inv(en.slope(e->B/2)/ratio);
   // Signed support margins continue below zero when the entire leaning
   // class disappears. Flat zero plateaus would generate spurious contours.
   double B=e->B,H=B/2,wB=hw(B,y,e->b,e->p2),wH=hw(H,y,e->b,e->p2),d0=hd(0,y,e->b,e->p2),dH=hd(H,y,e->b,e->p2),dB=hd(B,y,e->b,e->p2);
   double scale=1+wB+B*d0+B*dB;
   double A=(e->b-1)*(1-y),C=y*(e->b*e->p2-1);
   double c0=A+C,cB=A+C*std::pow(1+B,e->b);
   // Convex-then-concave: entry of an initial chord's contact into z<H.
   // Concave-then-convex: exit of a terminal chord's contact from z>H.
   ready1=(c0>0&&cB<0?wH-H*dH:B*d0-wB)/scale;
   ready2=(c0<0&&cB>0?(B-H)*dH-(wB-wH):wB-B*dB)/scale;
  }
  double l=std::log(std::max(cL,1e-300)),f=std::log(std::max(cF,1e-300)),h=std::log(std::max(cH,1e-300));
  score[0]=std::max(score[0],std::min(high,l)-low); // 2max
  score[1]=std::max(score[1],std::min(std::min(high,h)-std::max(low,l),ready2)); // 2-lean
  score[2]=std::max(score[2],std::min(std::min(high,f)-std::max(low,h),ready1)); // 1-lean
  score[3]=std::max(score[3],high-std::max(low,f)); // 1max
 }return score;
}
// [[Rcpp::export]]
NumericMatrix harvesting_crossings_cpp(NumericVector phi1,NumericVector phi2,NumericVector ratio,double beta2,double B,double median_alpha){
 int n=phi1.size();NumericMatrix out(n,7);double g=std::expm1(beta2*std::log1p(B/2));
 auto l1pe=[](double x){return x>30?x:std::log1p(std::exp(x));};
 for(int i=0;i<n;i++){
  double p=phi1[i],q=phi2[i],r=std::log(ratio[i]);
  auto f=[&](double t){double u=p*l1pe(t);double lw=t< -25?std::log(p)+t:(u>35?u:std::log(std::expm1(u)));return r+lw-std::log(q*beta2)-t-(beta2-1)*std::log1p(B/2)-(q-1)*l1pe(t+std::log(g));};
  std::vector<double> roots;double last_t=-40,last_y=f(last_t);
  for(int k=1;k<=160;k++){double t=-40+k*.5,y=f(t);if(y*last_y<0){double l=last_t,h=t,fl=last_y;for(int it=0;it<40;it++){double m=(l+h)/2,fm=f(m);if(fl*fm>0){l=m;fl=fm;}else h=m;}roots.push_back((l+h)/2);}last_t=t;last_y=y;}
  double asym=r-std::log(q*beta2)-(beta2-1)*std::log1p(B/2)-(q-1)*std::log(g);
  if(std::abs(p-q)>1e-12){double t=-asym/(p-q);if(t>40)roots.push_back(t);}
  int sign=std::abs(p-q)>1e-12?(p>q?1:-1):(asym>=0?1:-1);
  if(roots.empty())sign=f(0)>=0?1:-1;
  out(i,0)=roots.size();out(i,1)=roots.empty()?NA_REAL:roots.back();out(i,2)=sign;out(i,3)=0;out(i,4)=1;
  double a=median_alpha;
  double derivative=p*a*std::pow(1+a,p-1)/std::expm1(p*std::log1p(a))-1-(q-1)*a*g/(1+a*g);
  out(i,5)=f(std::log(a));out(i,6)=(!roots.empty() && derivative*sign>0)?1:0;
  if(!roots.empty()){
   a=std::exp(std::min(roots.back(),100.));HEnvelope env(a,beta2,q,B);
   out(i,3)=env.contact(B/2)?0:1;out(i,4)=B*hd(B/2,a,beta2,q)>=hw(B,a,beta2,q)?1:-1;
  }
 }return out;
}
