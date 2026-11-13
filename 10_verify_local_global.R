verify_local_global<-function(){
 tab<-read.csv('10_local_global_results/actor_budget_summaries.csv')
 stopifnot(nrow(tab)==24,all(abs(rowSums(tab[,c('training1','training2','campaign1','campaign2')])-1)<1e-7),
  all(tab[,c('training1','training2','campaign1','campaign2')]>=0))
 d<-readRDS('10_local_global_results/parameter_levels.rds')
 stopifnot(all(is.finite(unlist(d$fields))),all(unlist(d$fields)>=0),all(unlist(d$fields)<=1))
 for(z in d$fields)for(c in 1:3)for(r in 1:3)
  stopifnot(max(abs(z[,,c,r]-t(z[,,r,c])))<1e-10)
 files<-list.files('10_local_global_results/parameter_cells',full.names=TRUE)
 for(file in files){
  x<-readRDS(file);stopifnot(!is.list(x$result),abs(sum(x$result[c('campaign','training1','training2')])-1)<1e-7,
   abs(sum(x$result[c('campaign1','campaign2')])-x$result['campaign'])<1e-7)
 }
 message('Potential-actor budgets, campaign directions, parameter fields and exchange symmetry passed.')
 invisible(TRUE)
}
if(sys.nframe()==0L)verify_local_global()
