# Directly locate diversifier transitions between parameter-grid nodes.
# Fixed moments and responsive equilibria are evaluated at each bisection point.
# The original plot style is unchanged; supplied geometry caches permit redraws.
source('06_prepare_coexistence_deterministic.R')
refine_coexistence_curves<-function(plane=c('beta','phi'),threads=3L){
 plane<-match.arg(plane);input<-paste0('06_coexistence_',plane,'_inputs.rds')
 d<-prepare_coexistence(plane);file<-paste0('06_coexistence_',plane,'_geometry.rds')
 signature<-list(input=unname(tools::md5sum(input)),score=unname(tools::md5sum('06_coexistence_hdr_scores.cpp')),
  curve=unname(tools::md5sum('06_coexistence_curve_evaluations.cpp')),version=if(plane=='phi')3L else 2L)
 if(file.exists(file)){out<-readRDS(file);if(identical(out$signature,signature))return(out)}
 for(pkg in c('sf','Rcpp'))if(!requireNamespace(pkg,quietly=TRUE))stop('Install ',pkg)
 Rcpp::sourceCpp('06_coexistence_hdr_scores.cpp',cacheDir='coexistence_cpp_cache')
 # Inline the include for Rcpp's dependency scanner on Windows paths with spaces.
 hdr<-readLines('06_coexistence_hdr_scores.cpp',warn=FALSE)
 hdr[trimws(hdr)=='// [[Rcpp::export]]']<-'// Internal implementation.'
 extra<-readLines('06_coexistence_curve_evaluations.cpp',warn=FALSE)
 extra<-extra[!grepl('^#include "06_coexistence_hdr_scores.cpp"',extra)]
 Rcpp::sourceCpp(code=paste(c(hdr,extra),collapse='\n'),cacheDir='coexistence_curve_cpp_cache')
 Rcpp::sourceCpp('05_energy_allocation_population_phispace.cpp',cacheDir='energy_phispace_cpp_cache')
 axis<-d$axis;B<-d$B;sd<-d$talent$sdlog;hc<- -sd^2;hr<-sd*sqrt(qchisq(.99,2,ncp=2*sd^2))
 gh<-energy_gh(48,sd);gl<-energy_gl(32);mh<-energy_gh(64,sd);ml<-energy_gl(64)
 u<-seq(0,1,length.out=101);xt<-sort(unique(c(seq(0,B,length.out=181),B*u^2/2,B-B*u^2/2)))
 xl<-seq(0,B/2,length.out=121);xr<-seq(B/2,B,length.out=121);theta<-seq(-pi/2,pi/2,length.out=61)
 fields<-c('S2max','S2div','S1div','S1max');levels<-c(0,1e-4,1e-4,0)
 root<-file.path('coexistence_curve_checkpoints',plane);dir.create(root,recursive=TRUE,showWarnings=FALSE)
 evaluate<-function(x,y,i,j,mode,field=NULL,trial=FALSE){
  b1<-if(plane=='beta')x else rep(d$panel_vals[i],length(x));b2<-if(plane=='beta')y else rep(d$panel_vals[j],length(x))
  f1<-if(plane=='phi')x else rep(d$panel_vals[i],length(x));f2<-if(plane=='phi')y else rep(d$panel_vals[j],length(x))
  m1<-coex_random_moments_cpp(b1,f1,mh$a,mh$w,ml$x,ml$w,B,threads)
  m2<-coex_random_moments_cpp(b2,f2,mh$a,mh$w,ml$x,ml$w,B,threads);rho<-m2/m1
  if(mode=='responsive'){
   initial<-exp(coex_interp(log(d$panels[[i,j]]$rho),x,y,axis))
   if(trial)rho<-initial else{
    # Use the same exchange convention as the saved fields. In self-
    # symmetric panels the lower half is evaluated by reversing labels.
    reverse<-which((i==j & x>y)|(plane=='phi' & i==2 & j==3 & f1>f2))
    cb1<-b1;cb2<-b2;cf1<-f1;cf2<-f2;crho<-rho;cinitial<-initial
    cb1[reverse]<-b2[reverse];cb2[reverse]<-b1[reverse]
    cf1[reverse]<-f2[reverse];cf2[reverse]<-f1[reverse]
    crho[reverse]<-1/rho[reverse];cinitial[reverse]<-1/initial[reverse]
    z<-energy_population_cpp(cb1,cb2,cf1,cf2,crho,cinitial,gh$a,gh$w,gl$x,gl$w,B,sd,nx=101,threads=threads)
    if(any(z$status!=1))stop('Equilibrium failure during direct curve tracing')
    rho<-z$ratio_equilibrium;rho[reverse]<-1/rho[reverse]
   }
   sym<-which(i==j & abs(x-y)<1e-12);rho[sym]<-1
  }
  if(is.null(field))coexistence_scores_cpp(b1,b2,f1,f2,rho,xt,xl,xr,theta,hc,hr,B,threads)
  else coex_selected_scores_cpp(b1,b2,f1,f2,rho,rep(as.integer(field),length(x)),xt,xl,xr,theta,hc,hr,B,threads)
 }
 swap_lines<-function(ls)lapply(ls,function(z){a<-z$x;z$x<-z$y;z$y<-a;z})
 swap_sfg<-function(g){
  if(is.matrix(g))return(structure(g[,c(2,1),drop=FALSE],class=class(g)))
  structure(lapply(g,swap_sfg),class=class(g))
 }
 geometry<-list(fixed=matrix(vector('list',9),3,3),responsive=matrix(vector('list',9),3,3))
 for(i in 1:3)for(j in i:3)for(mode in c('fixed','responsive')){
  p<-d$panels[[i,j]][[mode]]
  lines<-Map(function(z,level)contourLines(axis,axis,z,levels=level),p[fields],levels)
  for(k in c(2,3)){
   if(i==j && k==3){lines[[3]]<-swap_lines(lines[[2]]);next}
   ls<-lines[[k]];if(!length(ls))next
   fn<-file.path(root,sprintf('panel%d_%d_%s_field%d.rds',i,j,mode,k))
   rec<-if(file.exists(fn))tryCatch(readRDS(fn),error=function(e)NULL)else NULL
   if(!is.null(rec)&&!identical(rec$signature,signature)){
    oldsig<-rec$signature;newsig<-signature
    unchanged<-(i!=j)
    if(plane=='phi'&&!is.null(d$compatible_geometry_input)&&identical(oldsig$input,d$compatible_geometry_input)){
     unchanged<- !(i==2&&j==3&&mode=='responsive')
     if(unchanged)oldsig$input<-newsig$input
    }
    oldsig$version<-NULL;newsig$version<-NULL
    if(unchanged&&identical(oldsig,newsig)){rec$signature<-signature;saveRDS(rec,fn)}
   }
   if(!is.null(rec)&&identical(rec$signature,signature)&&isTRUE(rec$complete)){lines[[k]]<-rec$lines;next}
   if(is.null(rec)||!identical(rec$signature,signature)){
    xy<-do.call(rbind,lapply(ls,function(z)cbind(z$x,z$y)))
    ix<-match(round(xy[,1],12),round(axis,12));iy<-match(round(xy[,2],12),round(axis,12))
    vertical<-!is.na(ix);n<-nrow(xy);lo<-hi<-numeric(n)
    fixed<-ifelse(vertical,xy[,1],xy[,2]);v<-ifelse(vertical,xy[,2],xy[,1])
    ind<-pmax(1,pmin(length(axis)-1,findInterval(v,axis)));lo<-axis[ind];hi<-axis[ind+1]
    xlo<-ifelse(vertical,fixed,lo);ylo<-ifelse(vertical,lo,fixed)
    xhi<-ifelse(vertical,fixed,hi);yhi<-ifelse(vertical,hi,fixed)
    flo<-coex_interp(p[[fields[k]]],xlo,ylo,axis)-levels[k]
    fhi<-coex_interp(p[[fields[k]]],xhi,yhi,axis)-levels[k]
    stopifnot(all(flo*fhi<=1e-12))
    rec<-list(signature=signature,iteration=0L,lo=lo,hi=hi,flo=flo,fhi=fhi,lo0=lo,hi0=hi,trial=(mode=='responsive'),vertical=vertical,fixed=fixed,ls=ls)
   }
   if(rec$iteration<9L)for(iter in seq.int(rec$iteration+1L,9L)){
    mid<-(rec$lo+rec$hi)/2;x<-ifelse(rec$vertical,rec$fixed,mid);y<-ifelse(rec$vertical,mid,rec$fixed)
    fm<-evaluate(x,y,i,j,mode,k,trial=isTRUE(rec$trial))-levels[k]
    left<-rec$flo*fm>0
    rec$lo[left]<-mid[left];rec$flo[left]<-fm[left]
    rec$hi[!left]<-mid[!left];rec$fhi[!left]<-fm[!left]
    rec$iteration<-iter;saveRDS(rec,fn)
    message('Direct ',plane,' curve ',i,',',j,' ',mode,' field ',k,': ',iter,'/9 (',length(mid),' vertices)')
   }
   if(isTRUE(rec$trial)){
    # The interpolated equilibrium is only a root-finding predictor. Accept
    # a bracket only after solving the actual equilibrium at both endpoints.
    actual<-function(v,ids){
     xx<-ifelse(rec$vertical[ids],rec$fixed[ids],v)
     yy<-ifelse(rec$vertical[ids],v,rec$fixed[ids])
     evaluate(xx,yy,i,j,mode,k,trial=FALSE)-levels[k]
    }
    ids<-seq_along(rec$lo);alo<-actual(rec$lo,ids);ahi<-actual(rec$hi,ids)
    bad<-which(alo*ahi>0)
    slope<-(rec$fhi-rec$flo)/(rec$hi-rec$lo)
    if(length(bad)){
     centre<-(rec$lo[bad]+rec$hi[bad])/2-(alo[bad]+ahi[bad])/(2*slope[bad])
     radius<-(rec$hi0[bad]-rec$lo0[bad])/1024
     good<-is.finite(centre)&is.finite(slope[bad])&abs(slope[bad])>1e-12
     centre[!good]<-(rec$lo0[bad[!good]]+rec$hi0[bad[!good]])/2
     rec$lo[bad]<-pmax(rec$lo0[bad],pmin(rec$hi0[bad],centre-radius))
     rec$hi[bad]<-pmax(rec$lo0[bad],pmin(rec$hi0[bad],centre+radius))
     alo[bad]<-actual(rec$lo[bad],bad);ahi[bad]<-actual(rec$hi[bad],bad)
    }
    bad<-which(alo*ahi>0)
    if(length(bad)){
     message('Full-model bisection fallback at ',length(bad),' vertices')
     rec$lo[bad]<-rec$lo0[bad];rec$hi[bad]<-rec$hi0[bad]
     alo[bad]<-actual(rec$lo[bad],bad);ahi[bad]<-actual(rec$hi[bad],bad)
     if(any(alo[bad]*ahi[bad]>1e-12))stop('Full-model root bracket failed')
     for(it in 1:9){
      mid<-(rec$lo[bad]+rec$hi[bad])/2;fm<-actual(mid,bad)
      left<-alo[bad]*fm>0
      il<-bad[left];ir<-bad[!left]
      rec$lo[il]<-mid[left];alo[il]<-fm[left]
      rec$hi[ir]<-mid[!left];ahi[ir]<-fm[!left]
     }
    }
    stopifnot(all(alo*ahi<=1e-12))
    rec$flo<-alo;rec$fhi<-ahi;rec$verified<-TRUE
   }
   mid<-(rec$lo+rec$hi)/2;xy<-cbind(ifelse(rec$vertical,rec$fixed,mid),ifelse(rec$vertical,mid,rec$fixed))
   at<-0L;ls<-lapply(rec$ls,function(z){ids<-at+seq_along(z$x);at<<-at+length(z$x);z$x<-xy[ids,1];z$y<-xy[ids,2];z})
   rec$lines<-ls;rec$complete<-TRUE;saveRDS(rec,fn);lines[[k]]<-ls
  }
  low<-min(axis);high<-max(axis)
  square<-sf::st_linestring(rbind(c(low,low),c(high,low),c(high,high),c(low,high),c(low,low)))
  sfg<-list(square)
  for(fam in lines)for(seg in fam)sfg<-c(sfg,list(sf::st_linestring(cbind(seg$x,seg$y))))
  network<-sf::st_node(sf::st_union(sf::st_sfc(sfg)))
  poly<-sf::st_collection_extract(sf::st_polygonize(network),'POLYGON')
  polys<-sf::st_sf(geometry=poly)
  xy<-sf::st_coordinates(suppressWarnings(sf::st_point_on_surface(polys)))
  vals<-as.matrix(as.data.frame(evaluate(xy[,1],xy[,2],i,j,mode)))
  present<-sweep(vals,2,levels,`>`);polys$rid<-rowSums(present*rep(c(1,2,4,8),each=nrow(present)))
  dissolved<-do.call(rbind,lapply(sort(unique(polys$rid)),function(id)sf::st_sf(rid=id,geometry=sf::st_union(polys$geometry[polys$rid==id]))))
  geom<-list(clines=lines,dissolved=dissolved);geometry[[mode]][[i,j]]<-geom
  if(i!=j){
   mirror<-geom;mirror$clines<-lapply(rev(lines),swap_lines)
   sf::st_geometry(mirror$dissolved)<-sf::st_sfc(lapply(sf::st_geometry(dissolved),swap_sfg))
   mirror$dissolved$rid<-vapply(dissolved$rid,function(v)sum(vapply(1:4,function(k)if(bitwAnd(as.integer(v),2^(k-1))>0)2^(4-k)else 0,1)),1)
   geometry[[mode]][[j,i]]<-mirror
  }
 }
 records<-lapply(list.files(root,full.names=TRUE,pattern='[.]rds$'),function(f)tryCatch(readRDS(f),error=function(e)NULL))
 records<-Filter(function(r)!is.null(r)&&identical(r$signature,signature)&&isTRUE(r$complete),records)
 brackets<-lapply(records,function(r)list(lo=r$lo,hi=r$hi,flo=r$flo,fhi=r$fhi))
 out<-list(signature=signature,geometry=geometry,edge_tolerance=max(diff(axis))/1024,brackets=brackets)
 saveRDS(out,file);out
}
