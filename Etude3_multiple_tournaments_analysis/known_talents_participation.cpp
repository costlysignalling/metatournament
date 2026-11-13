//Analytical support conditions and globally bounded branch optimization.
//The executable integrates indicators using deterministic shifted Halton rules.
#include <algorithm>
#include <array>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <queue>
#include <stdexcept>
#include <string>
#include <vector>
using namespace std;
struct Model {
 double beta,phi,B=50;
 double g(double x) const {return expm1(beta*log1p(x));}
 double w(double x,double a) const {return expm1(phi*log1p(a*g(x)));}
 double lq(double x,double a) const {
  return log(a)+(beta-1)*log1p(x)+(phi-1)*log1p(a*g(x));
 }
 double inv(double logc,double a,double lo,double hi) const {
  if(hi<=lo) return lo;
  if(logc>=lq(lo,a)) return lo;
  if(logc<=lq(hi,a)) return hi;
  double x=(lo+hi)/2;
  for(int n=0;n<45;n++){
   double lx=log1p(x), z=exp(beta*lx),h=1+a*(z-1);
   double f=log(a)+(beta-1)*lx+(phi-1)*log(h)-logc;
   if(f>0) lo=x; else hi=x;
   if(hi-lo<1e-9) break;
   double d=(beta-1)/(1+x)+(phi-1)*a*beta*z/((1+x)*h);
   double nx=x-f/d;
   if(!isfinite(nx)||nx<=lo||nx>=hi) nx=(lo+hi)/2;
   if(abs(nx-x)<1e-10) {x=nx;break;}
   x=nx;
  }
  return x;
 }
};
struct Sol {double value=-1;int count=0;vector<double>x;};
Sol concave(const Model&m,const vector<double>&a,const vector<double>&lo,
            const vector<double>&hi,double budget){
 double slo=0,shi=0;for(size_t i=0;i<a.size();i++){slo+=lo[i];shi+=hi[i];}
 if(budget<slo-1e-8||budget>shi+1e-8)return {};
 Sol out;out.x=lo;
 if(budget>slo+1e-10&&budget<shi-1e-10){
  double left=1e100,right=-1e100;
  for(size_t i=0;i<a.size();i++)if(hi[i]>lo[i]){
   left=min(left,m.lq(hi[i],a[i]));right=max(right,m.lq(lo[i],a[i]));
  }
  for(int j=0;j<48;j++){
   double c=(left+right)/2,s=0;
   for(size_t i=0;i<a.size();i++){out.x[i]=m.inv(c,a[i],lo[i],hi[i]);s+=out.x[i];}
   if(abs(s-budget)<1e-8)break;
   if(s>budget)left=c;else right=c;
  }
 }else if(budget>=shi-1e-10)out.x=hi;
 out.value=0;
 for(size_t i=0;i<a.size();i++){out.value+=m.w(out.x[i],a[i]);out.count+=(out.x[i]>1e-7);}
 return out;
}
int concave_count(const Model&m,const vector<double>&a,bool strong_bound=false){
 int T=a.size();
 auto at_least=[&](int k){
  double total=0,c=log(a[k-1]);
  for(int i=0;i<k-1;i++){
   double lo=(strong_bound&&i==0)?m.B/T:0;
   total+=m.inv(c,a[i],lo,m.B);
   if(total>=m.B-1e-9)return false;
  }
  return true;
 };
 if(at_least(T))return T;
 int lo=1,hi=T-1;
 while(lo<hi){int k=(lo+hi+1)/2;if(at_least(k))lo=k;else hi=k-1;}
 return lo;
}
struct Node {double l,r,upper;bool operator<(const Node&o)const{return upper<o.upper;}};
long long branch_cases=0,branch_nodes=0;
//One component is convex on [left,right]; all other components are concave.
//Its chord gives a valid upper bound. Maximize chord+concave profile exactly,
//then bisect until the global value gap is below the prescribed tolerance.
void convex_profile(const Model&m,const vector<double>&a,int convex_index,
                    const vector<double>&lo,const vector<double>&hi,
                    double left,double right,Sol&best){
 if(right<left+1e-12)return;
 branch_cases++;
 auto profile=[&](double t){
  Sol s=concave(m,a,lo,hi,m.B-t);
  if(s.value>=0){s.value+=m.w(t,a[convex_index]);s.x[convex_index]=t;s.count+=(t>1e-7);}
  if(s.value>best.value)best=s;
  return s;
 };
 profile(left);profile(right);
 auto bound=[&](double l,double r){
  double wl=m.w(l,a[convex_index]),wr=m.w(r,a[convex_index]);
  double slope=(wr-wl)/(r-l);
  double logc=log(slope/(m.phi*m.beta)),demand=0;
  for(size_t i=0;i<a.size();i++)demand+=m.inv(logc,a[i],lo[i],hi[i]);
  double t=max(l,min(r,m.B-demand));
  Sol s=profile(t);
  double upper=s.value-m.w(t,a[convex_index])+wl+slope*(t-l);
  return Node{l,r,upper};
 };
 priority_queue<Node>q;q.push(bound(left,right));int iterations=0;
 while(!q.empty()&&q.top().upper>best.value+1e-9*(1+abs(best.value))){
  Node n=q.top();q.pop();double mid=(n.l+n.r)/2;
  if(n.r-n.l<1e-10)continue;
  q.push(bound(n.l,mid));q.push(bound(mid,n.r));
  branch_nodes++;
  if(++iterations>100000)throw runtime_error("Global branch bound did not converge");
 }
}
int count_global(Model m,vector<double>a){
 sort(a.begin(),a.end(),greater<double>());int T=a.size();
 if(m.phi==1){
  double sum=0,p=1/(1-m.beta);int k=1;
  for(int j=0;j<T;j++){double z=pow(a[j]/a[0],p);sum+=z;
   if((m.B+j+1)*z>sum)k=j+1;else break;}
  return k;
 }
 if(m.phi<=1)return concave_count(m,a);
 if(m.beta==0.8&&m.phi==1.5){
  //Here w_xx >= 0 at every x whenever a>=1/2. The strongest
  //talent's convex payoff dominates all others and is superadditive.
  if(a[0]>=0.5)return 1;
  vector<double>lo(T,0),hi(T);
  for(int i=0;i<T;i++)hi[i]=min(m.B,pow(1/a[i]-1,1/m.beta)-1);
  Sol best=concave(m,a,lo,hi,m.B);
  //Any convex allocation must belong to the strongest talent:
  //allocations are talent-ordered and inflections decrease with talent.
  double inf=pow(1/a[0]-1,1/m.beta)-1;
  hi[0]=0;
  double othermax=0;for(double v:hi)othermax+=v;
  convex_profile(m,a,0,lo,hi,max({inf,m.B/T,m.B-othermax}),m.B,best);
  if(best.value<0)throw runtime_error("No feasible global solution");
  return best.count;
 }
 if(m.beta==0.4&&m.phi==1.5){
  int h=0;while(h<T&&a[h]>3)h++;
  if(h==0)return concave_count(m,a);
  //All inflections are below 1.756, while the strongest actor's
  //allocation is at least B/T>=50/7. It is on the concave branch.
  if(h==1)return concave_count(m,a,true);
  vector<double>inf(T,0);
  for(int i=0;i<h;i++)inf[i]=pow(1.5*(1-1/a[i]),1/m.beta)-1;
  Sol best;
  //At most one allocation can lie on a strictly convex branch:
  //otherwise transferring time between two such components increases value.
  for(int j=-1;j<h;j++){
   if(j==0)continue;
   for(int mask=0;mask<(1<<(h-1));mask++){
    if(j>0&&(mask&(1<<(j-1))))continue;
    vector<double>lo(T,0),hi(T,m.B);lo[0]=m.B/T;
    for(int i=1;i<h;i++){
     if(mask&(1<<(i-1)))lo[i]=inf[i];else hi[i]=0;
    }
    if(j<0){Sol s=concave(m,a,lo,hi,m.B);if(s.value>best.value)best=s;}
    else convex_profile(m,a,j,lo,hi,0,inf[j],best);
   }
  }
  if(best.value<0)throw runtime_error("No feasible global solution");
  return best.count;
 }
 throw runtime_error("Unsupported nonconcave parameter combination");
}
double halton(unsigned n,int base){double f=1,v=0;while(n){f/=base;v+=f*(n%base);n/=base;}return v;}
//Acklam inverse-normal rational approximation; error below 1e-8 here.
double normal_quantile(double p){
 const double a[]={-39.6968302866538,220.946098424521,-275.928510446969,138.357751867269,-30.6647980661472,2.50662827745924};
 const double b[]={-54.4760987982241,161.585836858041,-155.698979859887,66.8013118877197,-13.2806815528857};
 const double c[]={-.00778489400243029,-.322396458041136,-2.40075827716184,-2.54973253934373,4.37466414146497,2.93816398269878};
 const double d[]={.00778469570904146,.32246712907004,2.445134137143,3.75440866190742};
 if(p<.02425){double q=sqrt(-2*log(p));return (((((c[0]*q+c[1])*q+c[2])*q+c[3])*q+c[4])*q+c[5])/((((d[0]*q+d[1])*q+d[2])*q+d[3])*q+1);}
 if(p>1-.02425)return -normal_quantile(1-p);
 double q=p-.5,r=q*q;
 return (((((a[0]*r+a[1])*r+a[2])*r+a[3])*r+a[4])*r+a[5])*q/(((((b[0]*r+b[1])*r+b[2])*r+b[3])*r+b[4])*r+1);
}
int run_program(int argc,char**argv){
 try{
  if(argc>1&&string(argv[1])=="point"){
   Model m{stod(argv[2]),stod(argv[3])};vector<double>a;
   for(int i=4;i<argc;i++)a.push_back(stod(argv[i]));
   cout<<count_global(m,a)<<"\n";return 0;
  }
  int N=argc>1?stoi(argv[1]):65536,reps=argc>2?stoi(argv[2]):4;
  string output=argc>3?argv[3]:"participation_many.csv";
  ofstream out(output);out<<"beta,phi,T,rep,N,k,proportion\n";out<<setprecision(12);
  int primes[]={2,3,5,7,11,13,17};double slog=sqrt(log((1+sqrt(2.0))/2));
  for(double beta:{0.4,0.8})for(double phi:{0.66,1.0,1.5})for(int T=2;T<=7;T++){
   for(int rep=0;rep<reps;rep++){
    array<long long,8>counts{},small{};
    for(int n=1;n<=N;n++){
     vector<double>a(T);
     for(int j=0;j<T;j++){
      double shift=fmod(sqrt(double(primes[j]+23))*(rep+1),1.0);
      double u=fmod(halton(n+64,primes[j])+shift,1.0);
      u=max(1e-14,min(1-1e-14,u));a[j]=exp(slog*normal_quantile(u));
     }
     int k=count_global(Model{beta,phi},a);counts[k]++;
     if(n<=N/4)small[k]++;
    }
    for(int k=1;k<=T;k++){
     out<<beta<<','<<phi<<','<<T<<','<<rep<<','<<N<<','<<k<<','<<double(counts[k])/N<<'\n';
     out<<beta<<','<<phi<<','<<T<<','<<rep<<','<<N/4<<','<<k<<','<<double(small[k])/(N/4)<<'\n';
    }
    out.flush();
   }
   cerr<<"Finished beta="<<beta<<" phi="<<phi<<" T="<<T<<"\n";
  }
  cerr<<"Global convex-branch cases="<<branch_cases<<" nodes="<<branch_nodes<<"\n";
  return 0;
 }catch(const exception&e){cerr<<e.what()<<'\n';return 1;}
}
#ifdef RCPP_INTERFACE
#include <Rcpp.h>
// [[Rcpp::plugins(cpp17)]]
// [[Rcpp::export]]
void participation_compute(int N=65536,int reps=4,std::string output="participation_many.csv"){
 vector<string>args={"participation",to_string(N),to_string(reps),output};
 vector<char*>argv;for(auto &s:args)argv.push_back(&s[0]);
 if(run_program(argv.size(),argv.data())!=0)Rcpp::stop("Participation calculation failed");
}
// [[Rcpp::export]]
Rcpp::IntegerVector participation_count(Rcpp::NumericMatrix talents,double beta,double phi){
 Rcpp::IntegerVector result(talents.nrow());
 for(int i=0;i<talents.nrow();i++){
  vector<double>a(talents.ncol());for(int j=0;j<talents.ncol();j++)a[j]=talents(i,j);
  result[i]=count_global(Model{beta,phi},a);
 }
 return result;
}
#else
int main(int argc,char**argv){return run_program(argc,argv);}
#endif
