# Comparison only: replace each fixed-performance tournament by original
# training with beta=1, preserving phi, B, talent distributions and references.
# Uses the existing two-tournament conditional-integration population solver.
source('05_random_allocation_reference.R')
source('05_explorer_population.R')
source('09_strategy_examples.R')
training_comparison_core<-local({core<-NULL;function(){
 if(is.null(core)){
  code<-readLines('05_energy_allocation_population_phispace.cpp',warn=FALSE)
  wrapper<-'
// [[Rcpp::export]]
NumericVector comparison_training_map(NumericMatrix actors,double beta1,double beta2,
 double phi1,double phi2,double B,double ratio,int nx=801){
 ConditionalCell cell{beta1,beta2,phi1,phi2,B,.4,ratio,nx};
 NumericVector out(actors.nrow());
 for(int i=0;i<actors.nrow();i++)out[i]=cell.optimum(0,B,actors(i,0),actors(i,1),nx);
 return out;
}'
  core<<-new.env();Rcpp::sourceCpp(code=paste(c(code,wrapper),collapse='\n'),env=core,
    cacheDir=file.path(tempdir(),'09_training_comparison_cpp'))
 }
 core
}})
training_comparison_population<-function(p,nquad=64L,nx=801L){
 beta<-c(1,if(p$both_fixed)1 else p$beta[2])
 ref<-random_allocation_reference(beta[1],beta[2],p$phi[1],p$phi[2],p$B,nquad=nquad)
 pop<-explorer_population(beta[1],beta[2],p$phi[1],p$phi[2],p$B,ref$s,nquad=nquad,nx=nx)
 list(beta=beta,phi=p$phi,B=p$B,both_fixed=FALSE,
      fixed_ratio=ref$s,responsive_ratio=pop$responsive$s,
      first=unname(unlist(pop$first)[c('E1','E2','M1','M2')]),
      responsive=unname(unlist(pop$responsive)[c('E1','E2','M1','M2')]),
      residual=pop$responsive$residual,order=nquad)
}
training_comparison_alloc<-function(p,mode,actor,nx=801L){
 training_comparison_core()$comparison_training_map(as.matrix(actor),p$beta[1],p$beta[2],p$phi[1],p$phi[2],p$B,
    if(mode=='first')p$fixed_ratio else p$responsive_ratio,nx)
}
prepare_training_comparison_panel<-function(p,mode,a){
 B<-p$B;n<-length(a);alloc<-function(z)training_comparison_alloc(p,mode,z)
 q<-matrix(alloc(as.matrix(expand.grid(a,a))),n)
 boundary<-function(target){
  lo<-rep(min(a),n);hi<-rep(max(a),n)
  qlo<-alloc(cbind(lo,a));qhi<-alloc(cbind(hi,a));valid<-qlo<target&qhi>=target
  if(any(valid))for(k in 1:27){mid<-(lo+hi)/2;above<-alloc(cbind(mid[valid],a[valid]))>=target;ii<-which(valid);hi[ii[above]]<-mid[ii[above]];lo[ii[!above]]<-mid[ii[!above]]}
  out<-(lo+hi)/2;out[!valid&qlo>=target]<-min(a)-1;out[!valid&qhi<target]<-max(a)+1
  out
 }
 lower<-boundary(1e-5);upper<-boundary(B-1e-5)
 # Equal marginal returns: inverse in alpha1, then test global optimality.
 H<-B/2;ratio<-if(mode=='first')p$fixed_ratio else p$responsive_ratio
 marginal<-function(alpha,beta,phi)phi*alpha*beta*(1+H)^(beta-1)*(1+alpha*expm1(beta*log1p(H)))^(phi-1)
 target<-marginal(a,p$beta[2],p$phi[2])/ratio
 lo<-rep(0,n);hi<-rep(max(a),n);valid<-marginal(hi,p$beta[1],p$phi[1])>=target
 for(k in 1:40){mid<-(lo+hi)/2;above<-marginal(mid,p$beta[1],p$phi[1])>target;hi[above]<-mid[above];lo[!above]<-mid[!above]}
 eq<-(lo+hi)/2;eq[!valid]<-NA_real_
 ii<-which(is.finite(eq))
 if(length(ii)){
  xbest<-alloc(cbind(eq[ii],a[ii]))
  value<-function(x)ratio*boredom_w(x,eq[ii],p$beta[1],p$phi[1])+boredom_w(B-x,a[ii],p$beta[2],p$phi[2])
  eq[ii[value(H)<value(xbest)-1e-8*(1+abs(value(xbest)))]]<-NA_real_
 }
 list(q=q,field1=outer(a,upper,'-'),field2= -outer(a,lower,'-'),equal=list(list(x=eq,y=a)))
}
prepare_training_comparison<-function(overwrite=FALSE){
 file<-'09_harvesting_results/training_comparison.rds'
 if(file.exists(file)&&!overwrite)return(readRDS(file))
 original<-readRDS('09_harvesting_results/examples.rds')
 out<-lapply(names(original),function(set)lapply(seq_along(original[[set]]),function(j){message('Training comparison ',set,' ',j);training_comparison_population(original[[set]][[j]])}))
 names(out)<-names(original);saveRDS(out,file);out
}
draw_training_comparison<-function(){
 original<-readRDS('09_harvesting_results/examples.rds');alternative<-prepare_training_comparison()
 for(set in names(original)){
  draw_strategy_examples(original[[set]],paste0('Figure_harvesting_examples_',set))
  draw_strategy_examples(alternative[[set]],paste0('Figure_harvesting_comparison_training_',set),training=TRUE)
 }
 invisible(alternative)
}
