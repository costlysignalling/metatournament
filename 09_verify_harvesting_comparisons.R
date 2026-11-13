source('09_harvesting_model.R');source('08_finite_pool_model.R')
verify_harvesting_comparisons<-function(){
 examples<-readRDS('09_harvesting_results/examples.rds');out<-list();maxgap<-0
 for(set in names(examples))for(j in seq_along(examples[[set]])){
  p<-examples[[set]][[j]];fine<-harvest_population(p$beta[2],p$phi[1],p$phi[2],both_fixed=p$both_fixed,order=96L,inner=32L,start=p$responsive_ratio)
  err<-max(abs(c(p$first,p$responsive)[c(1,5)]-c(fine$first,fine$responsive)[c(1,5)]))/50
  stopifnot(err<5e-4,p$residual<1e-7)
  actor<-as.matrix(expand.grid(c(.1,.5,1,2.5,5),c(.1,.5,1,2.5,5)))
  for(mode in c('first','responsive')){
   ratio<-if(mode=='first')p$fixed_ratio else p$responsive_ratio;q<-harvest_alloc(p,mode,actor)
   for(k in seq_len(nrow(actor))){a1<-actor[k,1];a2<-actor[k,2]
    fun<-function(x)ratio*x*expm1(p$phi[1]*log1p(a1))+if(p$both_fixed)(50-x)*expm1(p$phi[2]*log1p(a2))else boredom_w(50-x,a2,p$beta[2],p$phi[2])
    xx<-seq(0,50,length.out=1001);vv<-fun(xx);peaks<-which(vv[2:1000]>=vv[1:999]&vv[2:1000]>=vv[3:1001])+1
    ref<-max(c(vv,vapply(peaks,function(i)optimize(fun,xx[c(i-1,i+1)],maximum=TRUE)$objective,0.)))
    gap<-(ref-fun(q[k]))/(1+abs(ref));maxgap<-max(maxgap,gap);stopifnot(gap<1e-8)
   }
  }
  out[[paste(set,j)]]<-c(time_share_difference=err,log_ratio_difference=abs(log(p$responsive_ratio/fine$responsive_ratio)),residual=p$residual)
  cat(set,j,'time share error',err,'\n')
 }
 # Finite-pool best responses against a denser independent allocation grid.
 fp<-readRDS('09_harvesting_results/finite_examples.rds')
 for(p in fp){stopifnot(p$responsive$converged,p$responsive$residual<2e-6)
  actor<-as.matrix(expand.grid(c(.1,.5,1,2.5),c(.1,.5,1,2.5)))
  for(mode in c('first','responsive')){
   q<-finite_pool_alloc(p,mode,actor);m<-if(mode=='first')p$reference else p$responsive$moment;bg<-(p$N-1)*c(m[1],50-m[1],m[2:3])
   for(k in 1:nrow(actor)){
    f<-function(x){w1<-boredom_w(x,actor[k,1],p$beta[1],p$phi[1]);w2<-boredom_w(50-x,actor[k,2],p$beta[2],p$phi[2]);(bg[1]+x)*w1/(bg[3]+w1)+(bg[2]+50-x)*w2/(bg[4]+w2)}
    grid<-seq(0,50,length.out=5001);gap<-(max(f(grid))-f(q[k]))/(1+abs(f(q[k])));stopifnot(gap<1e-8)
   }
  }
 }
 saveRDS(list(examples=out,max_global_utility_gap=maxgap),'09_harvesting_results/precision_checks.rds')
 cat('Largest relative global-utility gap:',maxgap,'\n')
 invisible(out)
}
if(sys.nframe()==0L)verify_harvesting_comparisons()
