# Reproducible higher-resolution parameter maps, with all three efficiencies
# sharing the expensive recipient response profile for each tournament pair.
# source('10_recompute_campaign_parameter_precision.R')
# recompute_campaign_parameter_precision()
source('10_campaign_dense_profile.R')
campaign_precision_tasks<-function(grid_size=25L){
 a<-seq(.1,.9,length.out=grid_size);phis<-c(.66,1,1.5);tasks<-list()
 for(mode in c('local','global'))for(c in 1:3)for(r in 1:c)for(i in seq_along(a))for(j in seq_along(a)){
  if(c==r&&j<i)next
  tasks[[length(tasks)+1L]]<-data.frame(mode=mode,phi1=phis[c],phi2=phis[r],beta1=a[i],beta2=a[j])
 }
 do.call(rbind,tasks)
}
campaign_precision_filename<-function(t,root='10_campaign_precision_results')
 sprintf('%s/cells/%s_p%.2f_%.2f_b%.6f_%.6f.rds',root,t$mode,t$phi1,t$phi2,t$beta1,t$beta2)
recompute_campaign_parameter_precision<-function(worker=1L,workers=1L,grid_size=25L,overwrite=FALSE){
 root<-'10_campaign_precision_results';dir.create(file.path(root,'cells'),recursive=TRUE,showWarnings=FALSE)
 tasks<-campaign_precision_tasks(grid_size)
 for(k in which((seq_len(nrow(tasks))-1L)%%workers==worker-1L)){
  t<-tasks[k,];file<-campaign_precision_filename(t)
  if(file.exists(file)&&!overwrite)next
  message(format(Sys.time()),' ',k,'/',nrow(tasks),' ',basename(file))
  ans<-tryCatch({
   p<-campaign_dense_profile(c(t$beta1,t$beta2),c(t$phi1,t$phi2),N=if(t$mode=='local')30 else 1e60)
   values<-tryCatch(campaign_precision_evaluate(p),error=function(err){
    message('  Escalating to full adaptive profile: ',conditionMessage(err))
    p<<-campaign_precision_profile(c(t$beta1,t$beta2),c(t$phi1,t$phi2),N=if(t$mode=='local')30 else 1e60)
    campaign_precision_evaluate(p)})
   for(v in values)stopifnot(all(is.finite(v)),all(v>=-1e-10),abs(sum(v[c('campaign','training1','training2')])-1)<1e-7,
    abs(v['campaign1']+v['campaign2']-v['campaign'])<1e-10)
   list(task=t,profile=p,values=values,grid_size=grid_size,version=2L)
  },error=function(err)list(error=conditionMessage(err),task=t))
  if(!is.null(ans$error)){saveRDS(ans,paste0(file,'.error'));message('FAILED: ',ans$error);next}
  temporary<-paste0(file,'.',Sys.getpid(),'.tmp');saveRDS(ans,temporary);file.rename(temporary,file)
 }
 message('WORKER COMPLETE')
 invisible(TRUE)
}
if(sys.nframe()==0L)recompute_campaign_parameter_precision()
