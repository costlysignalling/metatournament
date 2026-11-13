// Deterministic HDR coexistence scores; no sampled-population classification.
#include <Rcpp.h>
#include <cmath>
#include <algorithm>
#ifdef _OPENMP
#include <omp.h>
#endif
using namespace Rcpp;
// [[Rcpp::plugins(cpp11)]]
// [[Rcpp::plugins(openmp)]]
inline double gval(double x,double b){return std::expm1(b*std::log1p(x));}
inline double wval(double x,double a,double b,double p){return std::expm1(p*std::log1p(a*gval(x,b)));}
inline double marginal(double x,double a,double b,double p){
 return p*a*b*std::pow(1+x,b-1)*std::pow(1+a*gval(x,b),p-1);
}
inline bool concave(double a,double b,double p,double B){
 double k0=b-1+a*b*(p-1),kb=b-1+a*((b*p-1)*std::pow(1+B,b)-(b-1));
 return std::max(k0,kb)<=0;
}
inline bool convex(double a,double b,double p,double B){
 double k0=b-1+a*b*(p-1),kb=b-1+a*((b*p-1)*std::pow(1+B,b)-(b-1));
 return std::min(k0,kb)>=0;
}
inline double inverse_marginal(double target,double b,double p,double x){
 double G=gval(x,b),v=std::log(target/(p*b))-(b-1)*std::log1p(x),z=v;
 for(int k=0;k<30;++k){double ag=std::exp(z)*G;
  double step=(z+(p-1)*std::log1p(ag)-v)*(1+ag)/(1+p*ag);
  z-=step;if(std::abs(step)<1e-13)break;
 };return z;
}
inline double tie_a1(
    const double x,const double a2,
    const double b1,const double b2,
    const double p1,const double p2,
    const double rho,const double B){

  double lb=p2*std::log1p(a2*gval(B,b2));
  double lx=p2*std::log1p(a2*gval(B-x,b2));
  double d=std::exp(lb)*(-std::expm1(lx-lb));
  double den=gval(x,b1);
  if(den<=0)return R_PosInf;
  return std::expm1(std::log1p(d/rho)/p1)/den;
}

inline double critical_a1(
    const double a2,const double b1,const double b2,
    const double p1,const double p2,const double rho,
    const NumericVector& x_tie,const double B,double* entry=nullptr){

  // Exact x -> 0 endpoint-slope limit.
  double gb=gval(B,b2);
  double best=p2*a2*b2*std::pow(1.0+B,b2-1.0)*
    std::pow(1.0+a2*gb,p2-1.0)/(rho*p1*b1);

  if(entry)*entry=0;
  if(concave(best,b1,p1,B)&&concave(a2,b2,p2,B))return best;
  double endpoint=std::expm1(std::log1p(wval(B,a2,b2,p2)/rho)/p1)/gval(B,b1);
  if(convex(endpoint,b1,p1,B)&&convex(a2,b2,p2,B)){if(entry)*entry=B;return endpoint;}
  int bestk=-1;
  double bestgrid=R_PosInf;

  // Coarse global search.
  for(int k=1;k<x_tie.size();++k){
    double a=tie_a1(x_tie[k],a2,b1,b2,p1,p2,rho,B);

    if(a<bestgrid){
      bestgrid=a;
      bestk=k;
    }
  }

  if(bestgrid<best){best=bestgrid;if(entry)*entry=x_tie[bestk];}

  // Refine the winning interior interval continuously.
  if(bestk>=1 && bestk<x_tie.size()-1){
    double lo=std::max(1e-12,(double)x_tie[bestk-1]);
    double hi=x_tie[bestk+1];
    const double gr=0.6180339887498948482;

    double c=hi-gr*(hi-lo);
    double d=lo+gr*(hi-lo);
    double fc=tie_a1(c,a2,b1,b2,p1,p2,rho,B);
    double fd=tie_a1(d,a2,b1,b2,p1,p2,rho,B);

    for(int it=0;it<20;++it){
      if(fc>fd){
        lo=c;
        c=d;
        fc=fd;
        d=lo+gr*(hi-lo);
        fd=tie_a1(d,a2,b1,b2,p1,p2,rho,B);
      } else {
        hi=d;
        d=c;
        fd=fc;
        c=hi-gr*(hi-lo);
        fc=tie_a1(c,a2,b1,b2,p1,p2,rho,B);
      }
    }

    double ar=tie_a1(
      0.5*(lo+hi),a2,b1,b2,p1,p2,rho,B
    );

    if(ar<best){best=ar;if(entry)*entry=0.5*(lo+hi);}
  }

  return best;
}

// Signed absence margin: compare entry through the left half with the
// globally cheapest entry. Positive-width diversifier slices are unchanged.
inline double left_entry_a1(double a2,double b1,double b2,double p1,double p2,
 double rho,const NumericVector& grid,double B){
 double best=marginal(B,a2,b2,p2)/(rho*p1*b1);
 int kbest=-1,last=0;
 for(int k=1;k<grid.size() && grid[k]<=B/2+1e-12;++k){
  last=k;double v=tie_a1(grid[k],a2,b1,b2,p1,p2,rho,B);
  if(v<best){best=v;kbest=k;}
 }
 if(kbest>0){
  double lo=grid[std::max(0,kbest-1)],hi=grid[std::min(last,kbest+1)];
  lo=std::max(lo,1e-10);const double gr=.6180339887498948482;
  double c=hi-gr*(hi-lo),d=lo+gr*(hi-lo);
  double fc=tie_a1(c,a2,b1,b2,p1,p2,rho,B),fd=tie_a1(d,a2,b1,b2,p1,p2,rho,B);
  for(int k=0;k<28;++k){
   if(fc>fd){lo=c;c=d;fc=fd;d=lo+gr*(hi-lo);fd=tie_a1(d,a2,b1,b2,p1,p2,rho,B);}
   else{hi=d;d=c;fd=fc;c=hi-gr*(hi-lo);fc=tie_a1(c,a2,b1,b2,p1,p2,rho,B);}
  }
  best=std::min(best,std::min(fc,fd));
 }
 return best;
}

inline double indiv_u(
    const double x,const double a1,const double a2,
    const double b1,const double b2,
    const double p1,const double p2,
    const double rho,const double B){

  return rho*wval(x,a1,b1,p1)+wval(B-x,a2,b2,p2);
}

inline double refine_u_peak(
    const double lo0,const double hi0,
    const double a1,const double a2,
    const double b1,const double b2,
    const double p1,const double p2,
    const double rho,const double B){

  double lo=lo0,hi=hi0;
  const double gr=0.6180339887498948482;

  double c=hi-gr*(hi-lo);
  double d=lo+gr*(hi-lo);
  double fc=indiv_u(c,a1,a2,b1,b2,p1,p2,rho,B);
  double fd=indiv_u(d,a1,a2,b1,b2,p1,p2,rho,B);

  for(int it=0;it<22;++it){
    if(fc<fd){
      lo=c;
      c=d;
      fc=fd;
      d=lo+gr*(hi-lo);
      fd=indiv_u(d,a1,a2,b1,b2,p1,p2,rho,B);
    } else {
      hi=d;
      d=c;
      fd=fc;
      c=hi-gr*(hi-lo);
      fc=indiv_u(c,a1,a2,b1,b2,p1,p2,rho,B);
    }
  }

  return indiv_u(
    0.5*(lo+hi),a1,a2,b1,b2,p1,p2,rho,B
  );
}

inline double max_u_interval(
    const NumericVector& grid,
    const double a1,const double a2,
    const double b1,const double b2,
    const double p1,const double p2,
    const double rho,const double B){

  const int n=grid.size();
  if(convex(a1,b1,p1,B)&&convex(a2,b2,p2,B))
   return std::max(indiv_u(grid[0],a1,a2,b1,b2,p1,p2,rho,B),indiv_u(grid[n-1],a1,a2,b1,b2,p1,p2,rho,B));
  if(concave(a1,b1,p1,B)&&concave(a2,b2,p2,B)){
   double lo=grid[0],hi=grid[n-1];
   auto deriv=[&](double x){return rho*marginal(x,a1,b1,p1)-marginal(B-x,a2,b2,p2);};
   if(deriv(lo)<=0)return indiv_u(lo,a1,a2,b1,b2,p1,p2,rho,B);
   if(deriv(hi)>=0)return indiv_u(hi,a1,a2,b1,b2,p1,p2,rho,B);
   for(int k=0;k<35;++k){double mid=(lo+hi)/2;if(deriv(mid)>0)lo=mid;else hi=mid;}
   return indiv_u((lo+hi)/2,a1,a2,b1,b2,p1,p2,rho,B);
  }
  std::vector<double> vv(n);

  // Evaluate the full interval first.  This is important because in the
  // discontinuous-switching regime there may be several separated local
  // maxima of almost equal height.
  for(int k=0;k<n;++k){
    vv[k]=indiv_u(grid[k],a1,a2,b1,b2,p1,p2,rho,B);
  }

  double best=std::max(vv[0],vv[n-1]);

  // Refine EVERY local maximum, not merely the best coarse grid point.
  // The global maximum is then selected from the continuously refined peaks.
  for(int k=1;k<n-1;++k){
    if(vv[k]>=vv[k-1] && vv[k]>=vv[k+1]){
      double z=refine_u_peak(
        grid[k-1],grid[k+1],
        a1,a2,b1,b2,p1,p2,rho,B
      );
      if(z>best) best=z;
    }
  }

  // A maximum can lie between an interval endpoint and its nearest grid
  // point, even when the endpoint is the best coarse-grid value.
  double edge1=refine_u_peak(grid[0],grid[1],a1,a2,b1,b2,p1,p2,rho,B);
  double edge2=refine_u_peak(grid[n-2],grid[n-1],a1,a2,b1,b2,p1,p2,rho,B);
  return std::max(best,std::max(edge1,edge2));
}

inline double halfdiff(
    const double y1,const double a2,const double b1,const double b2,
    const double p1,const double p2,const double rho,
    const NumericVector& xL,const NumericVector& xR,const double B){

  double a1=std::exp(y1);

  double ml=max_u_interval(
    xL,a1,a2,b1,b2,p1,p2,rho,B
  );

  double mr=max_u_interval(
    xR,a1,a2,b1,b2,p1,p2,rho,B
  );

  return ml-mr;
}

inline void hroot_slice(
    const double L,const double U,const double a2,
    const double b1,const double b2,const double p1,const double p2,
    const double rho,const NumericVector& xL,const NumericVector& xR,
    const double B,double& root,int& status){

  double midroot=inverse_marginal(marginal(B/2,a2,b2,p2)/rho,b1,p1,B/2);
  bool continuous=concave(std::exp(U),b1,p1,B)&&concave(a2,b2,p2,B);
  if(!continuous && std::isfinite(midroot) && std::abs(midroot)<100){
   double a1=std::exp(midroot),umid=indiv_u(B/2,a1,a2,b1,b2,p1,p2,rho,B);
   double ml=max_u_interval(xL,a1,a2,b1,b2,p1,p2,rho,B);
   double mr=max_u_interval(xR,a1,a2,b1,b2,p1,p2,rho,B);
   continuous=std::max(ml,mr)<=umid+1e-12*(1+std::abs(umid));
  }
  if(continuous){root=midroot;status=root<=L?0:(root>=U?2:1);return;}
  double fL=halfdiff(L,a2,b1,b2,p1,p2,rho,xL,xR,B);
  double fU=halfdiff(U,a2,b1,b2,p1,p2,rho,xL,xR,B);

  if(fL<=0.0){
    root=L-1.0;
    status=0;
    return;
  }

  if(fU>=0.0){
    root=U+1.0;
    status=2;
    return;
  }

  double lo=L,hi=U;

  for(int it=0;it<32;++it){
    double mid=.5*(lo+hi);
    double fm=halfdiff(mid,a2,b1,b2,p1,p2,rho,xL,xR,B);

    if(fm>0.0) lo=mid;
    else hi=mid;
  }

  root=.5*(lo+hi);
  status=1;
}

inline void slice_scores_theta(
    const double th,
    const double b1,const double b2,const double p1,const double p2,
    const double rho,const NumericVector& x_tie,
    const NumericVector& xL,const NumericVector& xR,
    const double hc,const double hr,const double B,
    double& vm,double& width){

  double y2=hc+hr*std::sin(th);
  double d=hr*std::max(0.0,std::cos(th));
  double L=hc-d,U=hc+d;
  double a2=std::exp(y2);

  double entry=0;
  double bc=critical_a1(a2,b1,b2,p1,p2,rho,x_tie,B,&entry);
  double blog=std::log(std::max(bc,1e-300));

  vm=blog-L;

  if(blog>=U){width=U-blog;return;}
  // Increasing differences imply a nondecreasing global allocation in a1.
  // If entry already jumps to B/2 or farther, there is no 2-leaning interval.
  if(entry>=B/2){
   double left=left_entry_a1(a2,b1,b2,p1,p2,rho,x_tie,B);
   double margin=std::min(0.0,blog-std::log(std::max(left,1e-300)));
   width=std::min(margin,std::min(U,blog)-std::max(L,blog));return;
  }
  double h;
  int hstatus;
  hroot_slice(
    L,U,a2,b1,b2,p1,p2,rho,xL,xR,B,h,hstatus
  );

  if(hstatus==2) width=U-std::max(L,blog);
  else width=std::min(U,h)-std::max(L,blog);
}

inline double theta_score(
    const double th,const int which,
    const double b1,const double b2,const double p1,const double p2,
    const double rho,const NumericVector& x_tie,
    const NumericVector& xL,const NumericVector& xR,
    const double hc,const double hr,const double B){

  double vm,width;

  slice_scores_theta(
    th,b1,b2,p1,p2,rho,x_tie,xL,xR,hc,hr,B,vm,width
  );

  return which==0 ? vm : width;
}

inline double refine_theta_max(
    const double lo0,const double hi0,const int which,
    const double b1,const double b2,const double p1,const double p2,
    const double rho,const NumericVector& x_tie,
    const NumericVector& xL,const NumericVector& xR,
    const double hc,const double hr,const double B){

  double lo=lo0,hi=hi0;
  const double gr=0.6180339887498948482;

  double c=hi-gr*(hi-lo);
  double d=lo+gr*(hi-lo);

  double fc=theta_score(
    c,which,b1,b2,p1,p2,rho,x_tie,xL,xR,hc,hr,B
  );

  double fd=theta_score(
    d,which,b1,b2,p1,p2,rho,x_tie,xL,xR,hc,hr,B
  );

  for(int it=0;it<16;++it){
    if(fc<fd){
      lo=c;
      c=d;
      fc=fd;
      d=lo+gr*(hi-lo);
      fd=theta_score(
        d,which,b1,b2,p1,p2,rho,x_tie,xL,xR,hc,hr,B
      );
    } else {
      hi=d;
      d=c;
      fd=fc;
      c=hi-gr*(hi-lo);
      fc=theta_score(
        c,which,b1,b2,p1,p2,rho,x_tie,xL,xR,hc,hr,B
      );
    }
  }

  return theta_score(
    0.5*(lo+hi),which,b1,b2,p1,p2,rho,x_tie,xL,xR,hc,hr,B
  );
}

inline void scores_one(
    const double b1,const double b2,const double p1,const double p2,
    const double rho,const NumericVector& x_tie,
    const NumericVector& xL,const NumericVector& xR,
    const NumericVector& theta_grid,const double hc,const double hr,
    const double B,double& sm,double& sd){

  const int n=theta_grid.size();
  std::vector<double> vm(n),wd(n);

  for(int iy=0;iy<n;++iy){
    slice_scores_theta(
      theta_grid[iy],b1,b2,p1,p2,rho,
      x_tie,xL,xR,hc,hr,B,vm[iy],wd[iy]
    );
  }

  sm=std::max(vm[0],vm[n-1]);
  sd=std::max(wd[0],wd[n-1]);

  // Refine every local maximum along the HDR boundary.  This avoids
  // branch-selection artefacts when two different talent slices are nearly tied.
  for(int iy=1;iy<n-1;++iy){

    if(vm[iy]>=vm[iy-1] && vm[iy]>=vm[iy+1]){
      double z=refine_theta_max(
        theta_grid[iy-1],theta_grid[iy+1],0,
        b1,b2,p1,p2,rho,x_tie,xL,xR,hc,hr,B
      );
      if(z>sm) sm=z;
    }

    if(wd[iy]>=wd[iy-1] && wd[iy]>=wd[iy+1]){
      double z=refine_theta_max(
        theta_grid[iy-1],theta_grid[iy+1],1,
        b1,b2,p1,p2,rho,x_tie,xL,xR,hc,hr,B
      );
      if(z>sd) sd=z;
    }
  }
}



// [[Rcpp::export]]
List coexistence_scores_cpp(NumericVector beta1,NumericVector beta2,NumericVector phi1,
 NumericVector phi2,NumericVector rho,NumericVector x_tie,NumericVector xL,NumericVector xR,
 NumericVector theta,double hc,double hr,double B,int threads=6){
 int n=beta1.size();NumericMatrix out(n,4);
#ifdef _OPENMP
 omp_set_num_threads(threads);
#pragma omp parallel for schedule(dynamic)
#endif
 for(int k=0;k<n;++k){
  double sm,sd;
  scores_one(beta1[k],beta2[k],phi1[k],phi2[k],rho[k],x_tie,xL,xR,theta,hc,hr,B,sm,sd);
  out(k,0)=sm;out(k,1)=sd;
  scores_one(beta2[k],beta1[k],phi2[k],phi1[k],1/rho[k],x_tie,xL,xR,theta,hc,hr,B,sm,sd);
  out(k,2)=sd;out(k,3)=sm;
 }
 return List::create(_["S2max"]=out(_,0),_["S2div"]=out(_,1),_["S1div"]=out(_,2),_["S1max"]=out(_,3));
}
