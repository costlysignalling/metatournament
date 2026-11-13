// One persuader, one recipient-response round. Deterministic expected rivals.
// Flexible time, no compulsory attendance and no boredom.
#include <Rcpp.h>
#include <algorithm>
#include <cmath>
#include <map>
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
 double w2(double q,double a){return w((B-q),a,b2,p2);}
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
  auto deriv=[&](double q,bool right){double u=w1(q,a1),v=w2(q,a2),s1=u/(M1+u),s2=v/(M2+v);return s1+(E1+q)/(M1+u)*(1-s1)*wp(q,a1,b1,p1)-s2-((E2+B-q)/(M2+v)*(1-s2)*wp((B-q),a2,b2,p2));};
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

struct Persuasion {
 Local f; double N;
 Persuasion(double B,double b1,double b2,double p1,double p2,double n,int nx):f(B,.5,b1,b2,p1,p2,nx),N(n){}
 double part(double x,double a,int t,const NumericVector &bg) {
  double w=t==0?f.w1(x,a):f.w2(f.B-x,a);
  if(N>1e30)return bg[t]/bg[t+2]*w;
  return (bg[t]+x)*w/(bg[t+2]+w);
 }
 std::pair<double,double> receiver(double a1,double a2,double pressure,const NumericVector &bg,
    const std::vector<double>* u=nullptr,const std::vector<double>* v=nullptr) {
  double tilt=std::exp(pressure/2),other=std::exp(-pressure/2),bx=0,bv=-1; std::vector<double> vv(f.ng);
  auto value=[&](double x){return tilt*part(x,a1,0,bg)+other*part(f.B-x,a2,1,bg);};
  for(int j=0;j<f.ng;j++) {
   vv[j]=u?tilt*(*u)[j]+other*(*v)[j]:value(f.grid[j]);
   if(vv[j]>bv){bv=vv[j];bx=f.grid[j];}
  }
  std::vector<std::pair<double,double>> intervals;
  for(int j=1;j<f.ng-1;j++)if(vv[j]>=vv[j-1]&&vv[j]>=vv[j+1])intervals.push_back({f.grid[j-1],f.grid[j+1]});
  auto prime=[&](double x,double a,int t){
   double b=t==0?f.b1:f.b2,p=t==0?f.p1:f.p2,w=f.w(x,a,b,p),M=bg[t+2],den=M+w;
   if(N>1e30)return bg[t]/M*f.wp(x,a,b,p);
   return w/den+(bg[t]+x)*M*f.wp(x,a,b,p)/(den*den);
  };
  auto derivative=[&](double x){return tilt*prime(x,a1,0)-other*prime(f.B-x,a2,1);};
  if(derivative(0)>0&&derivative(f.grid[1])<0)intervals.push_back({0,f.grid[1]});
  if(derivative(f.B)<0&&derivative(f.grid[f.ng-2])>0)intervals.push_back({f.grid[f.ng-2],f.B});
  for(auto interval:intervals){
   double l=interval.first,h=interval.second,r=.618033988749895;
   double x=h-r*(h-l),y=l+r*(h-l),vx=value(x),vy=value(y);
   for(int k=0;k<32;k++)if(vx>vy){h=y;y=x;vy=vx;x=h-r*(h-l);vx=value(x);}else{l=x;x=y;vx=vy;y=l+r*(h-l);vy=value(y);}
   if(vx>bv){bv=vx;bx=x;}if(vy>bv){bv=vy;bx=y;}
  }
  return {bx,bv};
 }
};
// [[Rcpp::export]]
NumericMatrix convincing_receivers_cpp(NumericMatrix actor,NumericVector moment,NumericVector pressure,
 double beta1,double beta2,double phi1,double phi2,double B=50,double N=30,int nx=80){
 Persuasion z(B,beta1,beta2,phi1,phi2,N,nx);
 NumericVector bg=NumericVector::create((N-1)*moment[0],(N-1)*(B-moment[0]),(N-1)*moment[1],(N-1)*moment[2]);
 NumericMatrix out(actor.nrow(),pressure.size());
 for(int i=0;i<actor.nrow();i++){
  std::vector<double> u(z.f.ng),v(z.f.ng);
  for(int j=0;j<z.f.ng;j++){u[j]=z.part(z.f.grid[j],actor(i,0),0,bg);v[j]=z.part(B-z.f.grid[j],actor(i,1),1,bg);}
  for(int k=0;k<pressure.size();k++)out(i,k)=z.receiver(actor(i,0),actor(i,1),pressure[k],bg,&u,&v).first;
 }
 return out;
}
// [[Rcpp::export]]
NumericMatrix convincing_focal_cpp(NumericMatrix actor,NumericMatrix moments,NumericVector effort,
 double beta1,double beta2,double phi1,double phi2,double B=50,double N=30,int nx=40,bool refine=true,bool aligned=false,bool cache_claims=true){
 if(moments.ncol()!=4||moments.nrow()!=effort.size())stop("Recreate the preparation with the current tail-aware integration.");
 Persuasion z(B,beta1,beta2,phi1,phi2,N,nx);
 int na=actor.nrow(),nc=effort.size();NumericMatrix out(na,10);
 // Positive, shape-preserving interpolation. Log claims retain relative
 // accuracy in rare tails; logit investments retain both positive tails.
 std::vector<std::vector<double>> value(3,std::vector<double>(nc)),slope=value;
 for(int j=0;j<nc;j++){
  if(moments(j,0)<=0||moments(j,3)<=0||moments(j,1)<=0||moments(j,2)<=0)
   stop("An empty integration tail was found. Recalculate the stratified response table.");
  value[0][j]=std::log(moments(j,0))-std::log(moments(j,3));
  value[1][j]=std::log(moments(j,1));value[2][j]=std::log(moments(j,2));
 }
 for(int t=0;t<3;t++){
  std::vector<double> d(nc-1),h(nc-1);
  for(int j=0;j<nc-1;j++){h[j]=effort[j+1]-effort[j];d[j]=(value[t][j+1]-value[t][j])/h[j];}
  for(int j=1;j<nc-1;j++){
   if(d[j-1]*d[j]<=0)slope[t][j]=0;
   else {double w1=2*h[j]+h[j-1],w2=h[j]+2*h[j-1];slope[t][j]=(w1+w2)/(w1/d[j-1]+w2/d[j]);}
  }
  auto endpoint=[](double a,double b,double h0,double h1){double m=((2*h0+h1)*a-h0*b)/(h0+h1);if(m*a<=0)return 0.;if(a*b<0&&std::abs(m)>3*std::abs(a))return 3*a;return m;};
  slope[t][0]=endpoint(d[0],d[1],h[0],h[1]);
  slope[t][nc-1]=endpoint(d[nc-2],d[nc-3],h[nc-2],h[nc-3]);
 }
 auto background=[&](double c){
  int j=std::upper_bound(effort.begin(),effort.end(),c)-effort.begin()-1;j=std::max(0,std::min(nc-2,j));
  double h=effort[j+1]-effort[j],u=(c-effort[j])/h,u2=u*u,u3=u2*u;
  NumericVector bg(4);double v[3];
  for(int t=0;t<3;t++)v[t]=(2*u3-3*u2+1)*value[t][j]+(u3-2*u2+u)*h*slope[t][j]+(-2*u3+3*u2)*value[t][j+1]+(u3-u2)*h*slope[t][j+1];
  bg[0]=(N-1)*B/(1+std::exp(-v[0]));bg[1]=(N-1)*B/(1+std::exp(v[0]));
  bg[2]=(N-1)*std::exp(v[1]);bg[3]=(N-1)*std::exp(v[2]);return bg;
 };
 std::vector<NumericVector> bgs;for(double c:effort)bgs.push_back(background(c));
 std::vector<double> unique1(actor.column(0).begin(),actor.column(0).end()),unique2(actor.column(1).begin(),actor.column(1).end());
 std::sort(unique1.begin(),unique1.end());unique1.erase(std::unique(unique1.begin(),unique1.end()),unique1.end());
 std::sort(unique2.begin(),unique2.end());unique2.erase(std::unique(unique2.begin(),unique2.end()),unique2.end());
 int zero=std::min_element(effort.begin(),effort.end(),[](double a,double b){return std::abs(a)<std::abs(b);})-effort.begin();
 bool symmetric_grid=nc%2==1&&zero==nc/2;
 for(int j=0;j<nc&&symmetric_grid;j++)if(std::abs(effort[j]+effort[nc-1-j])>1e-9*B)symmetric_grid=false;
 double bytes=8.*(unique1.size()+unique2.size())*(zero+1)*(nx+1);
 bool cache=cache_claims&&symmetric_grid&&bytes<256.*1024*1024&&unique1.size()+unique2.size()<2*na;
 std::vector<std::vector<double>> w1,w2;
 if(cache){
  for(double a:unique1){std::vector<double> w;for(int k=zero;k<nc;k++)for(int j=0;j<=nx;j++)w.push_back(z.f.w((B-std::abs(effort[k]))*j/nx,a,beta1,phi1));w1.push_back(w);}
  for(double a:unique2){std::vector<double> w;for(int k=zero;k<nc;k++)for(int j=0;j<=nx;j++)w.push_back(z.f.w((B-std::abs(effort[k]))*(nx-j)/nx,a,beta2,phi2));w2.push_back(w);}
 }

 for(int i=0;i<na;i++){
  double a1=actor(i,0),a2=actor(i,1);
  auto basebr=z.receiver(a1,a2,0,bgs[zero]);double base=basebr.second,basex=basebr.first,preferred=2*basex-B;
  int i1=0,i2=0;if(cache){i1=std::lower_bound(unique1.begin(),unique1.end(),a1)-unique1.begin();i2=std::lower_bound(unique2.begin(),unique2.end(),a2)-unique2.begin();}
  auto val=[&](double c,double r){
   if(c< -B||c>B||r<0||r>1||(aligned&&c*preferred< -1e-7))return -1e100;
   auto bg=background(c);double remaining=B-std::abs(c),x=remaining*r,y=remaining*(1-r);
   double u=z.f.w(x,a1,beta1,phi1),v=z.f.w(y,a2,beta2,phi2);
   if(N>1e30)return bg[0]/bg[2]*u+bg[1]/bg[3]*v;
   return (bg[0]+x+std::max(c,0.))*u/(bg[2]+u)+(bg[1]+y+std::max(-c,0.))*v/(bg[3]+v);
  };
  double bc[2]={0,0},br[2]={basex/B,basex/B},bv[2]={base,base};
  std::vector<double> costval(nc,-1e100),costshare(nc,.5);
  // Search and refine each direction separately. Share coordinates transform
  // as r -> 1-r under relabelling, unlike an absolute-x pattern search.
  for(int k=0;k<nc;k++)for(int j=0;j<=nx;j++){
   double c=effort[k];if(aligned&&c*preferred< -1e-7)continue;
   double x=(B-std::abs(c))*j/nx,y=(B-std::abs(c))*(nx-j)/nx;
   double u=cache?w1[i1][std::abs(k-zero)*(nx+1)+j]:z.f.w(x,a1,beta1,phi1);
   double v=cache?w2[i2][std::abs(k-zero)*(nx+1)+j]:z.f.w(y,a2,beta2,phi2);
   const auto &bg=bgs[k];double vlu=N>1e30?bg[0]/bg[2]*u+bg[1]/bg[3]*v:(bg[0]+x+std::max(c,0.))*u/(bg[2]+u)+(bg[1]+y+std::max(-c,0.))*v/(bg[3]+v);
   if(vlu>costval[k]){costval[k]=vlu;costshare[k]=(double)j/nx;}
   for(int d=0;d<2;d++)if((d==0?c<=0:c>=0)&&vlu>bv[d]){bv[d]=vlu;bc[d]=c;br[d]=(double)j/nx;}
  }
  if(refine)for(int d=0;d<2;d++){
   // Keep the symmetric pattern-search candidate as an independent safeguard.
   double stepc=2*B/(nc-1),stepr=1./nx;
   for(int it=0;it<240&&(stepc>1e-7*B||stepr>1e-7);it++){
    double cc=bc[d],rr=br[d],vv=bv[d];
    for(int dc=-1;dc<=1;dc++)for(int dr=-1;dr<=1;dr++){
     double c=std::max(d==0?-B:0.,std::min(d==0?0.:B,bc[d]+dc*stepc));
     double r=std::max(0.,std::min(1.,br[d]+dr*stepr)),v=val(c,r);
     if(v>vv+1e-13){vv=v;cc=c;rr=r;}
    }
    if(vv>bv[d]+1e-13){bv[d]=vv;bc[d]=cc;br[d]=rr;}else{stepc*=.5;stepr*=.5;}
   }
   std::vector<int> peaks;int first=d==0?0:zero,last=d==0?zero:nc-1;
   for(int k=first;k<=last;k++)if((k==first||costval[k]>=costval[k-1])&&(k==last||costval[k]>=costval[k+1]))peaks.push_back(k);
   std::sort(peaks.begin(),peaks.end(),[&](int k,int l){return costval[k]>costval[l];});
   if(peaks.size()>4)peaks.resize(4);if(std::find(peaks.begin(),peaks.end(),zero)==peaks.end())peaks.push_back(zero);
   for(int seed:peaks){
    struct Point{double c,r,v;};
    double lo=d==0?-1.:0.,hi=d==0?0.:1.;
    auto point=[&](double c,double r){c=std::max(lo,std::min(hi,c));r=std::max(0.,std::min(1.,r));return Point{c,r,val(B*c,r)};};
    double cc=effort[seed]/B,rr=seed==zero?basex/B:costshare[seed];
    double dc=2./(nc-1)*(d==0?1.:-1.);if(cc+dc<lo||cc+dc>hi)dc=-dc;
    double dr=1./nx*(rr<=.5?1.:-1.);
    std::vector<Point> simplex={point(cc,rr),point(cc+dc,rr),point(cc,rr+dr)};
    for(int it=0;it<420;it++){
     std::sort(simplex.begin(),simplex.end(),[](Point a,Point b){return a.v>b.v;});
     double spread=0;for(int k=1;k<3;k++)spread=std::max(spread,std::max(std::abs(simplex[k].c-simplex[0].c),std::abs(simplex[k].r-simplex[0].r)));
     if(spread<1e-8&&simplex[0].v-simplex[2].v<1e-11*(1+std::abs(simplex[0].v)))break;
     double mc=(simplex[0].c+simplex[1].c)/2,mr=(simplex[0].r+simplex[1].r)/2;
     Point reflect=point(2*mc-simplex[2].c,2*mr-simplex[2].r);
     if(reflect.v>simplex[0].v){Point expand=point(mc+2*(reflect.c-mc),mr+2*(reflect.r-mr));simplex[2]=expand.v>reflect.v?expand:reflect;}
     else if(reflect.v>simplex[1].v)simplex[2]=reflect;
     else{
      bool outside=reflect.v>simplex[2].v;
      Point contract=outside?point(mc+.5*(reflect.c-mc),mr+.5*(reflect.r-mr)):point(mc+.5*(simplex[2].c-mc),mr+.5*(simplex[2].r-mr));
      if(contract.v>(outside?reflect.v:simplex[2].v))simplex[2]=contract;
      else for(int k=1;k<3;k++)simplex[k]=point((simplex[k].c+simplex[0].c)/2,(simplex[k].r+simplex[0].r)/2);
     }
    }
    for(auto candidate:simplex)if(candidate.v>bv[d]){bv[d]=candidate.v;bc[d]=B*candidate.c;br[d]=candidate.r;}
   }
  }
  double tol=1e-10*(1+std::max(bv[0],bv[1]));int best=bv[1]>bv[0]?1:0;
  bool tied=std::abs(bv[0]-bv[1])<=tol&&std::abs(bc[0]-bc[1])>1e-4*B;
  if(tied){
   double x0=(B-std::abs(bc[0]))*br[0],y0=(B-std::abs(bc[0]))*(1-br[0]);
   double x1=(B-std::abs(bc[1]))*br[1],y1=(B-std::abs(bc[1]))*(1-br[1]);
   double side=a1>=a2?1.:-1.;best=side*((x1-y1)-(x0-y0))>0?1:0;
  }
  double c=bc[best],r=br[best],v=bv[best],x=(B-std::abs(c))*r,y=(B-std::abs(c))*(1-r);
  if(v<=base+1e-8*(1+base)){c=0;x=basex;y=B-basex;v=base;tied=false;}
  out(i,0)=x;out(i,1)=y;out(i,2)=c;out(i,3)=v;out(i,4)=basex;out(i,5)=base;out(i,6)=v-base;
  out(i,7)=bv[1-best];out(i,8)=std::abs(bv[0]-bv[1]);out(i,9)=tied?1:0;
 }
 return out;
}

// Adaptive integration of recipient decisions, including strategy switches.
// Uses the recipient optimizer defined above.
#include <array>
#include <functional>
namespace convincing_quad {
using Vec=std::array<double,8>;
Vec add(Vec a,Vec b){for(int t=0;t<8;t++)a[t]+=b[t];return a;}
Vec scale(Vec a,double v){for(int t=0;t<8;t++)a[t]*=v;return a;}
struct Integrator {
 double tol; Vec unit;Vec error={0,0,0,0};long long evaluations=0;int unresolved=0;
 Integrator(double t,Vec s):tol(t),unit(s){}
 Vec apply(const std::function<Vec(double)>&f,double lo,double hi,int depth=0){
  static const double node[]={.9914553711208126,.9491079123427585,.8648644233597691,
   .7415311855993945,.5860872354676911,.4058451513773972,.2077849550078985,0};
  static const double kw[]={.02293532201052922,.06309209262997855,.1047900103222502,
   .1406532597155259,.1690047266392679,.1903505780647854,.2044329400752989,.2094821410847278};
  static const double gw[]={.1294849661688697,.2797053914892767,.3818300505051189,.4179591836734694};
  double mid=(lo+hi)/2,h=(hi-lo)/2;Vec center=f(mid),k=scale(center,kw[7]),g=scale(center,gw[3]);
  evaluations++;
  for(int j=0;j<7;j++){
   Vec v=add(f(mid-h*node[j]),f(mid+h*node[j]));evaluations+=2;
   k=add(k,scale(v,kw[j]));if(j%2==1)g=add(g,scale(v,gw[(j-1)/2]));
  }
  k=scale(k,h);g=scale(g,h);bool accurate=true;
  for(int t=0;t<8;t++)if(std::abs(k[t]-g[t])>tol*(std::abs(k[t])+unit[t]*1e-3))accurate=false;
  if(accurate||depth>=28){
   if(!accurate)unresolved++;
   for(int t=0;t<8;t++)error[t]+=std::abs(k[t]-g[t]);
   return k;
  }
  return add(apply(f,lo,mid,depth+1),apply(f,mid,hi,depth+1));
 }
 Vec whole(const std::function<Vec(double)>&f){
  // Explicit Gaussian tails prevent a central rule from omitting rare entrants.
  const double edge[]={-24,-20,-16,-12,-8,-6,-4,-2,0,2,4,6,8,12,16,20,24};Vec ans={0,0,0,0};
  for(int j=0;j<16;j++)ans=add(ans,apply(f,edge[j],edge[j+1]));return ans;
 }
};
}
// [[Rcpp::export]]
Rcpp::List convincing_adaptive_moments_cpp(Rcpp::NumericVector moment,
 Rcpp::NumericVector pressure,double beta1,double beta2,double phi1,double phi2,
 double B=50,double N=30,double slog=.434,int nx=80,double tolerance=1e-5,bool cache_claims=true){
 using namespace convincing_quad;
 Persuasion z(B,beta1,beta2,phi1,phi2,N,nx);
 Rcpp::NumericVector bg=Rcpp::NumericVector::create((N-1)*moment[0],
  (N-1)*(B-moment[0]),(N-1)*moment[1],(N-1)*moment[2]);
 Rcpp::NumericMatrix moments(pressure.size(),4),errors(pressure.size(),4);Rcpp::NumericVector evaluations(pressure.size());
 Rcpp::IntegerVector unresolved(pressure.size());
 Vec unit={B,std::max(1.,moment[1]),std::max(1.,moment[2]),B,1e6,1e6,1e6,1e6};
 // A grid claim depends on just one talent, not the other talent or pressure.
 // Reuse identical quadrature coordinates; eviction only recomputes values.
 std::map<double,std::vector<double>> claim_cache[2];
 auto grid_claims=[&](double normal,int t)->const std::vector<double>&{
  auto &cache=claim_cache[t];
  if(!cache_claims||cache.size()>100000)cache.clear();
  auto found=cache.find(normal);if(found!=cache.end())return found->second;
  double a=std::exp(slog*normal);std::vector<double> values(z.f.ng);
  for(int j=0;j<z.f.ng;j++)values[j]=z.part(t==0?z.f.grid[j]:B-z.f.grid[j],a,t,bg);
  return cache.emplace(normal,std::move(values)).first->second;
 };

 for(int k=0;k<pressure.size();k++){
  Rcpp::checkUserInterrupt();Integrator outer(tolerance,unit);long long inner_eval=0;int inner_failed=0;
  auto f2=[&](double s2){
   double a2=std::exp(slog*s2),weight2=std::exp(-s2*s2/2)/std::sqrt(2*M_PI);
   Integrator inner(tolerance/5,unit);
   // The tournament-2 grid terms are reused along the inner talent integral.
   const auto &v=grid_claims(s2,1);
   auto f1=[&](double s1){
    double a1=std::exp(slog*s1),weight1=std::exp(-s1*s1/2)/std::sqrt(2*M_PI);
    const auto &u=grid_claims(s1,0);
    double x=z.receiver(a1,a2,pressure[k],bg,&u,&v).first;
    return scale(Vec{x,z.f.w(x,a1,beta1,phi1),z.f.w(B-x,a2,beta2,phi2),B-x},weight1);
   };
   Vec ans=inner.whole(f1);inner_eval+=inner.evaluations;inner_failed+=inner.unresolved;
   for(int t=0;t<4;t++)ans[t+4]=inner.error[t];
   return scale(ans,weight2);
  };
  Vec ans=outer.whole(f2);for(int t=0;t<4;t++){moments(k,t)=ans[t];errors(k,t)=outer.error[t]+ans[t+4];}
  evaluations[k]=inner_eval;unresolved[k]=inner_failed+outer.unresolved;
 }
 return Rcpp::List::create(Rcpp::_["moments"]=moments,Rcpp::_["evaluations"]=evaluations,
  Rcpp::_["unresolved"]=unresolved,Rcpp::_["error_estimate"]=errors,Rcpp::_["tolerance"]=tolerance);
}



