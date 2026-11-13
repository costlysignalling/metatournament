# Rebuild the four main-text examples and their Supplement utility contours.
# Run from the scripts folder; C++ performs the allocation calculations.
source('08_boredom_model.R')
source('08_boredom_examples.R')
prepare_paper_boredom<-function(npop=32768L,nalpha=401L){
 scenarios<-list(list(beta=c(.2,.2),phi=c(2,2),N=30),
  list(beta=c(.8,.8),phi=c(2,2),N=30),
  list(beta=c(.8,.6),phi=c(1.5,1),N=30),
  list(beta=c(.8,.6),phi=c(1.5,1),N=10))
 den<-boredom_qmc_density(npop);alpha<-seq(0,2.5,length.out=nalpha)
 prep<-lapply(seq_along(scenarios),function(j){
  sc<-scenarios[[j]];message('Preparing example ',j)
  p<-boredom_population(sc$beta,sc$phi,50,.5,sc$N,den,160L,2e-6)
  if(!p$responsive$converged)stop('Example ',j,' did not converge: ',p$responsive$residual)
  p$maps<-lapply(c('first','responsive'),function(m)boredom_map(p,m,alpha,240L))
  p$boundaries<-lapply(c('first','responsive'),function(m)boredom_boundaries(p,m,alpha,240L))
  p$utility<-lapply(1:2,function(i){
   q<-p$maps[[i]];a1<-matrix(alpha,length(alpha),length(alpha));a2<-t(a1)
   moment<-if(i==1)p$reference else p$responsive$moment
   bg<-(p$N-1)*c(moment[1],50-moment[1],moment[2:3])
   w1<-boredom_w(q,a1,p$beta[1],p$phi[1]);w2<-boredom_w(pmin(25,50-q),a2,p$beta[2],p$phi[2])
   (bg[1]+q)*w1/(bg[3]+w1)+(bg[2]+50-q)*w2/(bg[4]+w2)
  })
  p
 })
 result<-list(prep=prep,alpha=alpha,B=50,compulsory=.5,N=c(30,30,30,10))
 dir.create('08_boredom_results',showWarnings=FALSE)
 saveRDS(result,'08_boredom_results/paper_examples.rds')
 invisible(result)
}
draw_paper_boredom<-function(){
 d<-readRDS('08_boredom_results/paper_examples.rds')
 plot_boredom_examples(data=d,paper=TRUE,stem='Figure_local_boredom_paper')
 plot_boredom_examples(data=d,contours=TRUE,paper=TRUE,stem='Figure_local_boredom_paper_utility_contours')
}
if(sys.nframe()==0L){prepare_paper_boredom();draw_paper_boredom()}
