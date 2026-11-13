# Integrate rare boredom conditionally instead of counting a handful of QMC
# talent points in its small region. Conditional roots delimit complete and
# partial disengagement; complete disengagement is integrated by its CDF.
source('08_boredom_model.R')
conditional_boredom<-function(p,mode='first',outer=24L,inner=12L){
 if(!exists('.boredom_conditional',.GlobalEnv)){
  env<-new.env();extension<-readLines('08_boredom_conditional.cpp')
  extension<-extension[!grepl('^#include "08_boredom_core.cpp"',extension)]
  Rcpp::sourceCpp(code=paste(c(readLines('08_boredom_core.cpp'),extension),collapse='\n'),env=env)
  assign('.boredom_conditional',env,.GlobalEnv)
 }
 gh<-boredom_gauss(outer);gl<-boredom_gauss(inner,FALSE)
 m<-if(mode=='first')p$reference else p$responsive$moment
 bg<-(p$N-1)*c(m[1],p$B-m[1],m[2:3])
 ans<-get('.boredom_conditional',.GlobalEnv)$boredom_conditional_cpp(bg,p$beta[1],p$beta[2],p$phi[1],p$phi[2],p$B,p$compulsory,
  exp(sqrt(2)*boredom_slog()*gh$x),gh$w,gl$x,2*gl$w,boredom_slog(),80L,22L)
 if(ans$failures)stop('Conditional integration boundary check failed.')
 ans$mean_bored/p$B
}
refine_boredom_level<-function(rows=1:3){
 d<-readRDS('08_boredom_results/energy_c_0.5_N_30.rds');den<-boredom_qmc_density(2048)
 for(j in rows)for(i in 1:3){
  file<-sprintf('08_boredom_results/level_panel_%s_%s.rds',i,j)
  if(file.exists(file))next
  p0<-d$panels[[i,j]];p<-p0;warm<-NULL
  p$moments<-vector('list',length(d$beta)^2)
  for(y in seq_along(d$beta)){
   for(x in seq_along(d$beta)){
    fit<-boredom_population(d$beta[c(x,y)],c(p$phi1,p$phi2),d$B,d$compulsory,d$N,den,80L,1e-5,warm,maxit=80L)
    if(!fit$responsive$converged)for(nhigh in c(8192L,32768L,131072L)){
     fit<-boredom_population(fit$beta,fit$phi,d$B,d$compulsory,d$N,boredom_qmc_density(nhigh),160L,1e-5,fit$responsive$moment,maxit=120L)
     if(fit$responsive$converged)break
    }
    if(!fit$responsive$converged)stop(sprintf('Population unresolved at panel %s,%s cell %s,%s: residual %.4g',i,j,x,y,fit$responsive$residual))
    warm<-fit$responsive$moment
    p$first$s0[x,y]<-conditional_boredom(fit,'first')
    p$responsive$s0[x,y]<-conditional_boredom(fit,'responsive')
    p$moments[[x+(y-1)*length(d$beta)]]<-list(reference=fit$reference,responsive=fit$responsive$moment,n=fit$order,residual=fit$responsive$residual)
   }
   cat('Conditional panel',i,j,'row',y,'\n');flush.console()
  }
  saveRDS(p,file)
 }
 files<-outer(1:3,1:3,Vectorize(function(i,j)sprintf('08_boredom_results/level_panel_%s_%s.rds',i,j)))
 if(all(file.exists(files))){d$panels<-lapply(as.vector(files),readRDS);dim(d$panels)<-c(3,3);d$method<-'conditional talent integration of boredom; QMC population background';saveRDS(d,'08_boredom_results/level_refined_c_0.5_N_30.rds')}
 invisible(NULL)
}
