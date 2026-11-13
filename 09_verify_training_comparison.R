source('09_training_comparison.R')
verify_training_comparison<-function(){
 original<-readRDS('09_harvesting_results/examples.rds')
 cmp<-readRDS('09_harvesting_results/training_comparison.rds')
 checks<-list();gapmax<-0
 actors<-as.matrix(expand.grid(c(.05,.2,.5,1,2.5,5),c(.05,.2,.5,1,2.5,5)))
 for(set in names(cmp))for(j in seq_along(cmp[[set]])){
  p<-cmp[[set]][[j]];fine<-training_comparison_population(original[[set]][[j]],nquad=96L,nx=1201L)
  delta<-max(abs(c(p$first[1]-fine$first[1],p$responsive[1]-fine$responsive[1])))/p$B
  stopifnot(delta<.001,abs(p$residual)<2e-6)
  for(mode in c('first','responsive')){
   q<-training_comparison_alloc(p,mode,actors)
   ratio<-if(mode=='first')p$fixed_ratio else p$responsive_ratio
   for(i in seq_len(nrow(actors))){
    value<-function(x)ratio*boredom_w(x,actors[i,1],p$beta[1],p$phi[1])+boredom_w(p$B-x,actors[i,2],p$beta[2],p$phi[2])
    grid<-seq(0,p$B,length.out=2001);v<-value(grid)
    peaks<-which(v[2:2000]>v[1:1999]&v[2:2000]>=v[3:2001])+1L
    best<-max(c(v,vapply(peaks,function(k)optimize(value,grid[c(k-1,k+1)],maximum=TRUE,tol=1e-10)$objective,0.)))
    gap<-(best-value(q[i]))/(1+abs(best));gapmax<-max(gapmax,gap)
    stopifnot(gap<1e-8)
   }
  }
  checks[[paste(set,j)]]<-c(time_share_change=delta,residual=p$residual)
  cat(set,j,'higher-resolution share change=',delta,'\n');flush.console()
 }
 saveRDS(list(checks=checks,max_global_utility_gap=gapmax),'09_harvesting_results/training_comparison_checks.rds')
 cat('Maximum independent global-utility gap:',gapmax,'\n')
 invisible(checks)
}
if(sys.nframe()==0L)verify_training_comparison()
