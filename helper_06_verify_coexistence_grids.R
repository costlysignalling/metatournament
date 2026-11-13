# Check saved coexistence data independently of the plotting routines.
fields<-c('S2max','S2div','S1div','S1max');levels<-c(0,1e-4,1e-4,0)
planes<-commandArgs(trailingOnly=TRUE);if(!length(planes))planes<-c('beta','phi')
for(plane in planes){
 d<-readRDS(paste0('06_coexistence_',plane,'_inputs.rds'));n<-length(d$axis)
 for(i in 1:3)for(j in 1:3){
  p<-d$panels[[i,j]];q<-d$panels[[j,i]];hit<-matrix(FALSE,n-1,n-1)
  stopifnot(all(is.finite(p$rho)),all(p$rho>0),max(abs(p$rho*t(q$rho)-1))<1e-10)
  for(mode in c('fixed','responsive'))for(k in 1:4){
   f<-p[[mode]][[fields[k]]]
   stopifnot(all(is.finite(f)),max(abs(f-t(q[[mode]][[fields[5-k]]])))<1e-9)
   f<-f-levels[k];z<-list(f[-n,-n],f[-1,-n],f[-n,-1],f[-1,-1]);lo<-Reduce(pmin,z);hi<-Reduce(pmax,z)
   hit<-hit|(lo<=0 & hi>=0 & hi-lo>1e-12)
  }
  direct<-p$direct[-n,-n]&p$direct[-1,-n]&p$direct[-n,-1]&p$direct[-1,-1]
  stopifnot(all(direct[hit]))
  cat(plane,'panel',i,j,':',sum(hit),'boundary cells directly covered\n')
 }
 gr<-expand.grid(x=d$axis,y=d$axis)
 tab<-do.call(rbind,lapply(1:3,function(i)do.call(rbind,lapply(1:3,function(j){
  p<-d$panels[[i,j]]
  do.call(rbind,lapply(c('fixed','responsive'),function(mode){
   z<-gr;names(z)<-if(plane=='beta')c('beta1','beta2')else c('phi1','phi2')
   if(plane=='beta'){z$phi1<-d$panel_vals[i];z$phi2<-d$panel_vals[j]}else{z$beta1<-d$panel_vals[i];z$beta2<-d$panel_vals[j]}
   z$mode<-mode
   for(k in fields)z[[k]]<-as.vector(p[[mode]][[k]])
   z$direct<-as.vector(p$direct);z
  }))
 }))))
 con<-gzfile(paste0('06_coexistence_',plane,'_scores.csv.gz'),'wt');write.csv(tab,con,row.names=FALSE);close(con)
}
cat('All saved-grid finiteness, exchange and direct-boundary-coverage checks passed.\n')
