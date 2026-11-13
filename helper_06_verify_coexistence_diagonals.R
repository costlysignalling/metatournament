Rcpp::sourceCpp('06_coexistence_hdr_scores.cpp',cacheDir='coexistence_cpp_cache')
B<-50;s<-sqrt(log((1+sqrt(2))/2));hc<- -s^2;hr<-s*sqrt(qchisq(.99,2,ncp=2*s^2))
u<-seq(0,1,length.out=101);xt<-sort(unique(c(seq(0,B,length.out=181),B*u^2/2,B-B*u^2/2)))
fields<-c('S2max','S2div','S1div','S1max')
planes<-commandArgs(trailingOnly=TRUE);if(!length(planes))planes<-c('beta','phi')
for(plane in planes){
 d<-readRDS(paste0('06_coexistence_',plane,'_inputs.rds'))
 for(i in 1:3){
  p<-d$panels[[i,i]];ids<-which(diag(p$direct));a<-d$axis[ids];v<-rep(d$panel_vals[i],length(ids))
  b<-if(plane=='beta')a else v;f<-if(plane=='phi')a else v
  z<-coexistence_scores_cpp(b,b,f,f,rep(1,length(ids)),xt,seq(0,25,length.out=121),seq(25,50,length.out=121),seq(-pi/2,pi/2,length.out=61),hc,hr,B,3)
  error<-max(vapply(fields,function(k)max(abs(diag(p$fixed[[k]])[ids]-z[[k]]),abs(diag(p$responsive[[k]])[ids]-z[[k]])),1))
  stopifnot(error<1e-7)
  cat(plane,'diagonal panel',i,': direct-score error',error,'\n')
 }
}
