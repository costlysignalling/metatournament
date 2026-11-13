// Continuous endpoint expectations plus a stratified QMC diversification residual.
// Used for phi1=phi2=1.5. For beta*phi<=1, the residual covers all talents.
// Above a_t=(1-beta_t)/(beta_t*(phi_t-1)), each tournament weight is convex
// for the whole investment interval. If both talents exceed those thresholds,
// the global optimum is an endpoint and the residual is identically zero.
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
inline double weight(double x,double a,double b,double p) {
  return std::expm1(p*std::log1p(a*std::expm1(b*std::log1p(x))));
}
inline double phi_cdf(double a,double sd) {
  if(a<=0) return 0;
  return R::pnorm(std::log(a)/sd,0,1,true,false);
}
struct CVState { double e1,m1,m2; };
struct CVCell {
  double b1,b2,p,B,sd,G1,G2;
  std::vector<double> a1,a2,pw,w1,w2;
  int nx;
  CVCell(double beta1,double beta2,double budget,double s,
         const NumericVector& u,const NumericVector& v,int ngrid):
    b1(beta1),b2(beta2),p(1.5),B(budget),sd(s),nx(ngrid) {
    G1=std::expm1(b1*std::log1p(B));G2=std::expm1(b2*std::log1p(B));
    double q1=b1*p>1 ? phi_cdf((1-b1)/(b1*(p-1)),sd) : 1.0;
    double q2=b2*p>1 ? phi_cdf((1-b2)/(b2*(p-1)),sd) : 1.0;
    int n=u.size();
    if(q1>1e-15) for(int i=0;i<n;++i) {
      a1.push_back(std::exp(sd*R::qnorm(q1*u[i],0,1,true,false)));
      a2.push_back(std::exp(sd*R::qnorm(v[i],0,1,true,false)));
      pw.push_back(q1/n);
    }
    if((1-q1)*q2>1e-15) for(int i=0;i<n;++i) {
      a1.push_back(std::exp(sd*R::qnorm(q1+(1-q1)*u[i],0,1,true,false)));
      a2.push_back(std::exp(sd*R::qnorm(q2*v[i],0,1,true,false)));
      pw.push_back((1-q1)*q2/n);
    }
    w1.resize(a1.size()*nx);w2.resize(a1.size()*nx);
    for(size_t i=0;i<a1.size();++i) for(int k=0;k<nx;++k) {
      double x=B*k/(nx-1.0);
      w1[i*nx+k]=weight(x,a1[i],b1,p);
      w2[i*nx+k]=weight(B-x,a2[i],b2,p);
    }
  }
  double utility(double x,size_t i,double ratio) const {
    return ratio*weight(x,a1[i],b1,p)+weight(B-x,a2[i],b2,p);
  }
  double refine(double lo,double hi,size_t i,double ratio) const {
    const double gr=.6180339887498948482;
    double c=hi-gr*(hi-lo),d=lo+gr*(hi-lo);
    double fc=utility(c,i,ratio),fd=utility(d,i,ratio);
    for(int k=0;k<25;++k) {
      if(fc<fd) {lo=c;c=d;fc=fd;d=lo+gr*(hi-lo);fd=utility(d,i,ratio);}
      else {hi=d;d=c;fd=fc;c=hi-gr*(hi-lo);fc=utility(c,i,ratio);}
    }
    return (lo+hi)/2;
  }
  CVState state(double ratio,const NumericVector& talents,const NumericVector& qw) const {
    CVState out{0,0,0};
    // Three smooth one-dimensional integrals. Swapping the integration
    // variable gives the tournament-1 tail weight without a nested integral.
    for(int j=0;j<talents.size();++j) {
      double a=talents[j];
      double wB1=weight(B,a,b1,p),wB2=weight(B,a,b2,p);
      double threshold1=std::expm1(std::log1p(wB2/ratio)/p)/G1;
      double threshold2=std::expm1(std::log1p(ratio*wB1)/p)/G2;
      double F1=phi_cdf(threshold1,sd),F2=phi_cdf(threshold2,sd);
      out.e1+=qw[j]*B*(1-F1);
      out.m1+=qw[j]*wB1*F2;
      out.m2+=qw[j]*wB2*F1;
    }
    for(size_t i=0;i<a1.size();++i) {
      double endpoint1=ratio*w1[i*nx+nx-1],endpoint2=w2[i*nx];
      bool choose1=endpoint1>endpoint2;
      double xe=choose1?B:0;
      int bestk=choose1?nx-1:0;
      double best=std::max(endpoint1,endpoint2);
      // For beta*phi>1, w'' changes at most once, from negative to
      // positive. Consequently w(x)/x has no interior maximum, and
      // w(x) <= x*max(w'(0),w(B)/B). If the resulting linear utility
      // bound cannot beat the best endpoint, the residual is exactly zero.
      double bound1=ratio*std::max(B*p*a1[i]*b1,w1[i*nx+nx-1]);
      double bound2=std::max(B*p*a2[i]*b2,w2[i*nx]);
      if(b1*p>1 && b2*p>1 && bound1<=best && bound2<=best) continue;
      for(int k=1;k<nx-1;++k) {
        double val=ratio*w1[i*nx+k]+w2[i*nx+k];
        if(val>best) {best=val;bestk=k;}
      }
      double x=xe;
      if(bestk>0 && bestk<nx-1) {
        double candidate=refine(B*(bestk-1)/(nx-1.0),B*(bestk+1)/(nx-1.0),i,ratio);
        double val=utility(candidate,i,ratio);
        if(val>std::max(endpoint1,endpoint2)) {x=candidate;best=val;}
      }
      // Also bracket narrow endpoint-adjacent maxima with endpoint slopes.
      double d0=ratio*p*a1[i]*b1-
        p*a2[i]*b2*std::pow(1+B,b2-1)*std::pow(1+a2[i]*G2,p-1);
      double dB=ratio*p*a1[i]*b1*std::pow(1+B,b1-1)*std::pow(1+a1[i]*G1,p-1)-p*a2[i]*b2;
      if(d0>0) {
        double candidate=refine(0,B/(nx-1.0),i,ratio),val=utility(candidate,i,ratio);
        if(val>best) {x=candidate;best=val;}
      }
      if(dB<0) {
        double candidate=refine(B-B/(nx-1.0),B,i,ratio),val=utility(candidate,i,ratio);
        if(val>best) {x=candidate;best=val;}
      }
      out.e1+=pw[i]*(x-xe);
      out.m1+=pw[i]*(weight(x,a1[i],b1,p)-(choose1?w1[i*nx+nx-1]:0));
      out.m2+=pw[i]*(weight(B-x,a2[i],b2,p)-(choose1?0:w2[i*nx]));
    }
    return out;
  }
};
// [[Rcpp::export]]
List refine_equal_talent_nash_cpp(NumericVector beta1,NumericVector beta2,
    NumericVector initial,NumericVector u,NumericVector v,
    NumericVector talents,NumericVector weights,double B,double sd,
    int nx=81,int threads=4) {
  int n=beta1.size();NumericVector ratios(n),residuals(n),widths(n),lower_residuals(n),upper_residuals(n);
  IntegerVector status(n);
#ifdef _OPENMP
  omp_set_num_threads(threads);
#pragma omp parallel for schedule(dynamic)
#endif
  for(int i=0;i<n;++i) {
    if(beta1[i]<=0 || beta2[i]<=0 || beta1[i]>=1 || beta2[i]>=1) {status[i]=-1;continue;}
    CVCell cell(beta1[i],beta2[i],B,sd,u,v,nx);
    auto f=[&](double z) {
      CVState st=cell.state(std::exp(z),talents,weights);
      return std::log(st.e1)+std::log(st.m2)-std::log(B-st.e1)-std::log(st.m1)-z;
    };
    double centre=std::log(initial[i]),lo=centre-.2,hi=centre+.2;
    double flo=f(lo),fhi=f(hi);
    if(flo*fhi>0) {lo=centre-2;hi=centre+2;flo=f(lo);fhi=f(hi);}
    if(!std::isfinite(flo)||!std::isfinite(fhi)||flo*fhi>0) {status[i]=0;continue;}
    // Bracketed Illinois secant iteration: retain a root bracket, but avoid
    // spending 32 full population integrations on a smooth scalar equation.
    double root=(lo+hi)/2;
    int previous_side=0;
    for(int k=0;k<80;++k) {
      double candidate=(lo*fhi-hi*flo)/(fhi-flo);
      if(!std::isfinite(candidate) || candidate<=lo || candidate>=hi)
        candidate=(lo+hi)/2;
      double fm=f(candidate);
      root=candidate;
      if(!std::isfinite(fm)) break;
      if(std::abs(fm)<1e-11 || hi-lo<1e-11) break;
      if(flo*fm>0) {
        lo=candidate;flo=fm;
        if(previous_side==1) fhi*=.5;
        previous_side=1;
      } else {
        hi=candidate;fhi=fm;
        if(previous_side==-1) flo*=.5;
        previous_side=-1;
      }
    }
    ratios[i]=std::exp(root);residuals[i]=f(root);
    widths[i]=hi-lo;lower_residuals[i]=f(lo);upper_residuals[i]=f(hi);
    // A finite weighted QMC rule may switch a whole quadrature atom between
    // two tied optima. A tiny sign-changing bracket then locates the switch
    // even though its one-sided residuals jump over zero. Keep this distinct
    // from a smooth root, and check integration resolution separately.
    if(std::isfinite(residuals[i]) && std::abs(residuals[i])<1e-8) status[i]=1;
    else if(widths[i]<1e-10 && lower_residuals[i]*upper_residuals[i]<=0) status[i]=2;
    else status[i]=0;
  }
  return List::create(_["ratio"]=ratios,_["residual"]=residuals,_["status"]=status,_["log_bracket_width"]=widths,
    _["lower_residual"]=lower_residuals,_["upper_residual"]=upper_residuals);
}
