# Run from the scripts folder. Source this file, then change the function inputs.
# In interactive R, sourcing draws the single example in `settings` below.
# Like Script 05: plot in R by default; legend, PNG and PDF are optional.
# Scalar inputs beta1,beta2,phi1,phi2 and vectors beta=c(...),phi=c(...) work.
# Fixed: full compulsory engagement, uniform free-time split, independent talents.
# Responsive: self-consistent expected competitor density with N-1 rivals.
# Talent distribution stays lognormal, median=1, SD=.5. Boredom transfers pool
# contribution to tournament 2 without improving performance there.
# C++ is preferred and reused between calls; R fallback is available if blocked.
# Increase npop and nx/nxpop to check numerical precision, especially at high phi.
#
# explore_boredom(beta=c(.8,.8), phi=c(1.5,2))
# explore_boredom(beta1=.6,beta2=.8,phi1=1,phi2=2,N=15,compulsory=.6)
# explore_boredom(beta=c(.8,.8),phi=c(2,2),mode='fixed',legend=TRUE)
# result <- explore_boredom(beta=c(.8,.8),phi=c(2,2.25),png=TRUE,pdf=TRUE)
# result$summary; result$population$responsive$residual
# plot(result,legend=TRUE) # redraw without recalculating
# saveRDS(result,'my_boredom_example.rds'); plot(readRDS('my_boredom_example.rds'))

settings<-list(beta1=.8,beta2=.8,phi1=1.5,phi2=2,B=50,
 compulsory=.5,N=30,mode='both',alim=c(0,2.5),nalpha=251L,
 npop=8192L,nx=160L,nxpop=80L,max_npop=65536L,
 legend=FALSE,png=FALSE,pdf=FALSE)

source('08_boredom_model.R')
source('08_boredom_examples.R')
if(!exists('.boredom_explorer_cache',envir=.GlobalEnv,inherits=FALSE))
 assign('.boredom_explorer_cache',new.env(parent=emptyenv()),envir=.GlobalEnv)

explore_boredom<-function(beta1=.8,beta2=.8,phi1=1.5,phi2=2,B=50,
 compulsory=.5,N=30,mode=c('both','fixed','responsive'),alim=c(0,2.5),
 nalpha=251L,npop=8192L,nx=160L,nxpop=80L,max_npop=65536L,
 tolerance=1e-5,maxit=120L,legend=FALSE,png=FALSE,pdf=FALSE,
 file=NULL,beta=NULL,phi=NULL,verbose=TRUE){
 mode<-match.arg(mode)
 if(!is.null(beta)){if(!is.numeric(beta)||length(beta)!=2)stop('beta must have two values.');beta1<-beta[1];beta2<-beta[2]}
 if(!is.null(phi)){if(!is.numeric(phi)||length(phi)!=2)stop('phi must have two values.');phi1<-phi[1];phi2<-phi[2]}
 pars<-list(beta1,beta2,phi1,phi2,B,compulsory,N,tolerance)
 if(!all(vapply(pars,function(x)is.numeric(x)&&length(x)==1&&is.finite(x),TRUE)))stop('Model parameters must be finite numeric scalars.')
 if(min(beta1,beta2)<=0||max(beta1,beta2)>=1)stop('Require 0 < beta1,beta2 < 1.')
 if(min(phi1,phi2,B,tolerance)<=0||compulsory<=0||compulsory>=1||N<=1)stop('Require positive phi, B, tolerance; 0 < compulsory < 1; N > 1.')
 if(!is.numeric(alim)||length(alim)!=2||any(!is.finite(alim))||alim[1]<0||diff(alim)<=0)stop('alim must be an increasing nonnegative interval.')
 for(x in list(nalpha,nx,nxpop,npop,max_npop,maxit))if(!is.numeric(x)||length(x)!=1||!is.finite(x)||x!=as.integer(x)||x<3)stop('Resolutions and maxit must be integer scalars of at least 3.')
 if(npop<32||npop%%2||max_npop<npop||max_npop%%2)stop('npop must be even and at least 32; max_npop must be even and at least npop.')
 for(x in list(legend,png,pdf,verbose))if(!is.logical(x)||length(x)!=1||is.na(x))stop('Plot/output flags must be TRUE or FALSE.')
 # Reject numerical overflow rather than silently painting invalid strategies.
 largest_log_weight<-c(phi1,phi2)*log1p(alim[2]*expm1(c(beta1,beta2)*log1p(B)))
 if(any(!is.finite(largest_log_weight))||any(largest_log_weight>log(.Machine$double.xmax)-8))stop('This talent range/phi combination exceeds double precision. Reduce alim or phi.')
 ref_check<-boredom_reference(c(beta1,beta2),c(phi1,phi2),B,compulsory)
 if(any(!is.finite(ref_check))||any(ref_check<=0)||any(!is.finite(ref_check[1]/ref_check[2:3])))stop('The reference moments exceed numerical precision for these parameters.')

 signature<-list(beta=c(beta1,beta2),phi=c(phi1,phi2),B=B,c=compulsory,N=N,mode=mode,
  npop=npop,max_npop=max_npop,nxpop=nxpop,tolerance=tolerance,maxit=maxit,
  model=unname(tools::md5sum(c('08_boredom_model.R','08_boredom_core.cpp','08_boredom_R_fallback.R'))))
 cache<-get('.boredom_explorer_cache',envir=.GlobalEnv)
 reused<-exists('signature',cache,inherits=FALSE)&&identical(cache$signature,signature)
 if(reused){population<-cache$population;den<-boredom_qmc_density(population$order)}else{
  den<-boredom_qmc_density(npop)
  if(mode=='fixed'){
   reference<-boredom_reference(c(beta1,beta2),c(phi1,phi2),B,compulsory)
   first<-boredom_compile()$boredom_density_cpp(den$talent,den$weights,reference,beta1,beta2,phi1,phi2,B,compulsory,N,nxpop,1L,tolerance,FALSE)
   population<-list(beta=c(beta1,beta2),phi=c(phi1,phi2),B=B,compulsory=compulsory,N=N,reference=reference,first=first,order=npop)
  }else{
   population<-boredom_population(c(beta1,beta2),c(phi1,phi2),B,compulsory,N,den,nxpop,tolerance,maxit=maxit)
   while(!population$responsive$converged&&den$order<max_npop){
    den<-boredom_qmc_density(min(max_npop,4L*den$order))
    if(verbose)message('Retrying responsive density with ',den$order,' integration points.')
    population<-boredom_population(c(beta1,beta2),c(phi1,phi2),B,compulsory,N,den,nxpop,tolerance,start=population$responsive$moment,maxit=maxit)
   }
   if(!population$responsive$converged)stop(sprintf('Responsive density did not converge (residual %.3g). Increase npop/max_npop or nxpop; no equilibrium map was drawn.',population$responsive$residual))
  }
  cache$signature<-signature;cache$population<-population
 }
 modes<-switch(mode,both=c('first','responsive'),fixed='first',responsive='responsive')
 alpha<-seq(alim[1],alim[2],length.out=nalpha)
 maps<-setNames(lapply(modes,function(m)boredom_map(population,m,alpha,nx)),modes)
 boundaries<-list();monotone<-vapply(maps,function(q)min(q[-1,,drop=FALSE]-q[-nrow(q),,drop=FALSE])>= -1e-4*B,TRUE)
 for(m in modes)boundaries[[m]]<-if(monotone[m])boredom_boundaries(population,m,alpha,nx)else NULL
 if(any(!monotone))warning('Nonmonotone talent response detected. Those boundaries use the full strategy grid; increase nalpha for more detail.',call.=FALSE)
 equal<-list();core<-boredom_compile();equal_q<-if(compulsory>.5)(1-compulsory)*B else B/2
 for(m in modes){
  moment<-if(m=='first')population$reference else population$responsive$moment
  bg<-(N-1)*c(moment[1],B-moment[1],moment[2:3])
  candidate<-if(monotone[m])list(boundaries[[m]][[2]])else contourLines(alpha,alpha,maps[[m]],levels=equal_q)
  equal[[m]]<-lapply(candidate,function(branch){
   ok<-is.finite(branch$x)&is.finite(branch$y);ids<-which(ok)
   if(length(ids)){
    a1<-branch$x[ids];a2<-branch$y[ids]
    qbest<-core$boredom_map_cpp(cbind(a1,a2),bg,beta1,beta2,phi1,phi2,B,compulsory,nx)
    value<-function(q){w1<-boredom_w(q,a1,beta1,phi1);w2<-boredom_w(pmin((1-compulsory)*B,B-q),a2,beta2,phi2);(bg[1]+q)*w1/(bg[3]+w1)+(bg[2]+B-q)*w2/(bg[4]+w2)}
    ok[ids]<-value(equal_q)>=value(qbest)-1e-7*(1+B)
   }
   branch$x[!ok]<-NA_real_;branch
  })
 }
 summary<-do.call(rbind,lapply(modes,function(m){st<-population[[m]];q<-st$q;v<-st$mapped
  data.frame(condition=if(m=='first')'Fixed E_t: first response'else 'Responsive E_t',productive_1=v[1],productive_2=B-v[1]-v[4],bored_time=v[4],pool_1=v[1],pool_2=B-v[1],actors_bored_pct=100*sum(den$weights[q<compulsory*B-1e-5]),actors_fully_bored_pct=100*sum(den$weights[q<1e-5]))}))
 rownames(summary)<-NULL
 result<-structure(list(population=population,alpha=alpha,maps=maps,boundaries=boundaries,
  monotone=monotone,equal=equal,summary=summary,settings=list(B=B,compulsory=compulsory,N=N,
  beta=c(beta1,beta2),phi=c(phi1,phi2),npop=den$order,nx=nx,nxpop=nxpop,mode=mode,tolerance=tolerance),files=character()),class='boredom_exploration')
 if(is.null(file))file<-paste0('Figure_boredom_explorer_beta_',paste(c(beta1,beta2),collapse='_'),'_phi_',paste(c(phi1,phi2),collapse='_'),'_c_',compulsory,'_N_',N,'_B_',B,'_',mode)
 if(!is.character(file)||length(file)!=1||is.na(file)||!nzchar(file))stop('file must be a nonempty filename or stem.')
 file<-sub('[.](png|pdf)$','',file,ignore.case=TRUE)
 width<-4*length(modes)+if(legend)1.7 else 0;height<-4.8
 if(png){grDevices::png(paste0(file,'.png'),width=width,height=height,units='in',res=600);tryCatch(plot(result,legend=legend),finally=dev.off());result$files<-c(result$files,png=paste0(file,'.png'))}
 if(pdf){grDevices::cairo_pdf(paste0(file,'.pdf'),width=width,height=height,family='Arial');tryCatch(plot(result,legend=legend),finally=dev.off());result$files<-c(result$files,pdf=paste0(file,'.pdf'))}
 if(interactive()||dev.cur()!=1L)plot(result,legend=legend)
 if(verbose){message('Population averages (time units per actor; pool 2 includes bored time):');print(summary,digits=5,row.names=FALSE);if(mode!='fixed')message(sprintf('Responsive residual: %.3g; integration points: %d%s',population$responsive$residual,den$order,if(reused)' (population reused)'else ''))}
 invisible(result)
}

plot.boredom_exploration<-function(x,legend=FALSE,...){
 stopifnot(is.logical(legend),length(legend)==1,!is.na(legend))
 a<-x$alpha;alim<-range(a);p<-x$settings;B<-p$B;c<-p$compulsory;nm<-length(x$maps)
 layout(matrix(seq_len(nm+as.integer(legend)),1),widths=c(rep(1,nm),if(legend).45 else numeric()))
 par(oma=c(0,0,1.7,0),mar=c(4,4,5,1.6),mgp=c(2.5,.7,0),pty='s',cex=1)
 for(j in seq_along(x$maps)){
  mode<-names(x$maps)[j];q<-x$maps[[j]];cl<-x$boundaries[[mode]];cols<-boredom_example_colours(q,B,c)
  if(x$monotone[mode])cols[q>=B-1e-5]<-'#B54343' # rooted polygon supplies the crisp maxxer fill
  if(!x$monotone[mode]&&abs(c-.5)<1e-8)cols[abs(q-B/2)<1e-5]<-'#FF0'
  plot(NA,xlim=alim,ylim=alim,xaxs='i',yaxs='i',xlab=expression(alpha[1]~'(talent in tournament 1)'),ylab=if(j==1)expression(alpha[2]~'(talent in tournament 2)')else '')
  rasterImage(as.raster(t(cols)[ncol(cols):1,,drop=FALSE]),min(a),min(a),max(a),max(a),interpolate=TRUE)
  if(x$monotone[mode]){
   # Rooted polygon edges, rather than pixel-mask outlines.
   edge<-cl[[4]]$x;edge[!is.finite(edge)]<-max(a);edge[q[1,]>=B-1e-5]<-min(a)
   polygon(c(edge,max(a),max(a)),c(a,max(a),min(a)),col='#E06',border=NA)
   edge0<-cl[[1]]$x;edge0[!is.finite(edge0)]<-min(a);edge0[q[nrow(q),]<=1e-5]<-max(a)
   zero_col<-boredom_example_colours(matrix(0,1,1),B,c)[1]
   polygon(c(min(a),edge0,min(a)),c(min(a),a,max(a)),col=zero_col,border=NA)
   if(abs(c-.5)<1e-8){lo<-cl[[2]]$x;hi<-cl[[3]]$x;hi[!is.finite(hi)&abs(q[nrow(q),]-B/2)<1e-4]<-max(a);valid<-which(is.finite(lo)&is.finite(hi)&hi>=lo)
    if(length(valid)>1)for(ii in split(valid,cumsum(c(1,diff(valid)!=1))))if(length(ii)>1)polygon(c(lo[ii],rev(hi[ii])),c(a[ii],rev(a[ii])),col='#FF0',border=NA)
   }
  }
  for(z in boredom_hdr())lines(z[1,],z[2,],col='#EEE',lwd=.8)
  abline(0,1,lty=2,lwd=.8)
  if(x$monotone[mode]){lines(cl[[1]]$x,a,col='black',lty=2,lwd=1);lines(cl[[4]]$x,a,col='#E06',lwd=1)}else for(v in c(1e-5,B-1e-5))for(z in contourLines(a,a,q,levels=v))lines(z$x,z$y,lty=if(v<B/2)2 else 1,col=if(v<B/2)'black'else '#E06')
  if(abs(c-.5)>1e-8||!x$monotone[mode])for(z in x$equal[[mode]])lines(z$x,z$y,col='#FF0',lwd=3)
  box();title(main=if(mode=='first')expression('Fixed '*E[t])else expression('Responsive '*E[t]),line=3.7,cex.main=1.15)
  for(t in 1:2){at<-min(a)+c(.3,.7)[t]*diff(alim);mtext(bquote(beta[.(t)]==.(sprintf('%.2f',p$beta[t]))),3,line=1.7,at=at,cex=1.05,col=c('#E06','#08F')[t]);mtext(bquote(varphi[.(t)]==.(sprintf('%.2f',p$phi[t]))),3,line=.3,at=at,cex=1.05,col=c('#E06','#08F')[t])}
  m<-x$population[[mode]]$mapped;yy<-min(a)+c(0,m[4],m[4]+m[1],B)*diff(alim)/B;x0<-max(a)+.025*diff(alim);x1<-max(a)+.075*diff(alim)
  for(k in 1:3)rect(x0,yy[k],x1,yy[k+1],col=c('black','#E06','#08F')[k],border=NA,xpd=NA)
  rect(x0,min(a),x1,max(a),lwd=.7,xpd=NA)
 }
 mtext(bquote(B==.(B)~~c==.(sprintf('%.2f',c))~~N==.(p$N)),3,outer=TRUE,line=.2,cex=.95)
 if(legend){
  par(mar=c(4,4,5,1.6),pty='m',xpd=NA);plot.new();plot.window(xlim=c(0,1.12),ylim=c(0,1),xaxs='i',yaxs='i')
  gr<-expand.grid(d=seq(0,1,length.out=121),r=seq(0,1,length.out=121));v<-colorRamp(c('#315A9B','white','#B54343'))(gr$r)/255*(1-gr$d)+(85/255)*gr$d;cm<-matrix(rgb(v[,1],v[,2],v[,3]),121,121)
  rasterImage(as.raster(t(cm)[121:1,,drop=FALSE]),0,0,1,1,interpolate=TRUE);box();axis(1,at=c(0,.5,1),cex.axis=.9);axis(2,at=c(0,.5,1),las=1,cex.axis=.9)
  rect(0,.49,.12,.51,col='#FF0');rect(0,.98,.12,1,col='#E06');rect(0,0,.12,.02,col='#08F');rect(1,0,1.1,1,col='black')
  mtext(expression(d/B),1,line=2,cex=.9);mtext('Bored time share',1,line=3,cex=.8);mtext('Tournament-1 share of productive time',2,line=2.4,cex=.8)
 }
 invisible(x)
}

if(interactive())do.call(explore_boredom,settings)



#This will be in the paper
explore_boredom(beta1=.2,beta2=.2,phi1=2,phi2=2,B=50,compulsory=.5,N=30)
explore_boredom(beta1=.8,beta2=.8,phi1=2,phi2=2,B=50,compulsory=.5,N=30)
explore_boredom(beta1=.8,beta2=.6,phi1=1.5,phi2=1,B=50,compulsory=.5,N=30)
explore_boredom(beta1=.8,beta2=.6,phi1=1.5,phi2=1,B=50,compulsory=.5,N=10)


explore_boredom(beta1=.6,beta2=.8,phi1=1,phi2=1.5,B=50,compulsory=.5,N=30)
explore_boredom(beta1=.6,beta2=.8,phi1=1,phi2=1.5,B=50,compulsory=.5,N=10)


