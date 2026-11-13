source('10_local_global_figures.R')
source('10_convincing_effort_limits.R')
campaign_entry_curves<-function(d,threshold=1e-5*d$p$B){
 a<-d$map$alpha;n<-length(a);mask<-matrix(abs(d$map$allocation[,3])>threshold,n)
 ij<-which(mask[-1,,drop=FALSE]!=mask[-n,,drop=FALSE],arr.ind=TRUE)
 if(!nrow(ij))return(list())
 lo<-a[ij[,1]];hi<-a[ij[,1]+1];y<-a[ij[,2]];side<-mask[cbind(ij[,1]+1,ij[,2])]
 for(k in 1:12){mid<-(lo+hi)/2;v<-convincing_optimize(d$p,cbind(mid,y))[,3]
  move<-abs(v)>threshold;move<-move==side;hi[move]<-mid[move];lo[!move]<-mid[!move]}
 roots<-lapply(seq_len(n),function(j)(lo+hi)[ij[,2]==j]/2)
 distance<-matrix(0,n,n)
 for(j in seq_len(n))distance[,j]<-if(length(roots[[j]]))
  apply(abs(outer(a,roots[[j]],'-')),1,min)*ifelse(mask[,j],1,-1)else
  ifelse(mask[,j],1,-1)*diff(range(a))
 contourLines(a,a,distance,levels=0)
}
prepare_local_global_figures<-function(png=TRUE,pdf=TRUE){
 dir.create('10_local_global_results',showWarnings=FALSE)
 report<-list()
 limits_file<-'10_local_global_results/effort_limits.rds'
 if(!file.exists(limits_file))campaign_limit_summary()
 limits<-readRDS(limits_file)
 for(family in c('local','global'))for(set in c('symmetric_beta04','symmetric_beta08','asymmetric')){
  message('Preparing ',family,' ',set)
  if(family=='local'){
   data<-readRDS(paste0('10_convincing_results/',set,'.rds'))
   geometry<-readRDS(paste0('10_convincing_results/',set,'_geometry.rds'))
  }else{
   data<-readRDS(paste0('10_comparison_results/popularizing_',set,'.rds'))
   geometry<-readRDS(paste0('10_comparison_results/geometry_popularizing_',set,'.rds'))
  }
  for(j in seq_along(data)){
   fine<-paste0('10_local_global_results/high_resolution/',family,'_',set,'_',j,'.rds')
   if(file.exists(fine))data[[j]]$map<-readRDS(fine)
  }
  entry_file<-paste0('10_local_global_results/entry_',family,'_',set,'.rds')
  if(file.exists(entry_file))entry<-readRDS(entry_file)else{
   entry<-lapply(data,campaign_entry_curves);saveRDS(entry,entry_file)
  }
  for(j in seq_along(data)){
   p<-data[[j]]$p;o<-p$actors;ww<-p$density$weights
   mean_budget<-colSums(ww*cbind(o[,1:2],pmax(o[,3],0),pmax(-o[,3],0)))/p$B
   opposed<-abs(o[,3])>1e-4*p$B&o[,3]*(o[,1]-o[,2])< -1e-5
   report[[paste(family,set,j)]]<-data.frame(family=family,set=set,panel=j,
    training1=mean_budget[1],training2=mean_budget[2],campaign1=mean_budget[3],
    campaign2=mean_budget[4],opposed_training_and_campaign=sum(ww*opposed),
    relative_campaign1=if(sum(mean_budget[3:4])>1e-12)mean_budget[3]/sum(mean_budget[3:4])else NA_real_,
    reallocated_recipient_budget=unname(p$summary['mean_recipient_abs_shift'])/p$B)
  }
  effort_geometry<-if(family=='local')lapply(seq_along(data),function(j){
   file<-paste0('10_local_global_results/effort_geometry_local_',set,'_',j,'.rds')
   if(!file.exists(file))stop('Prepare campaign-direction geometry first: ',file)
   readRDS(file)})else NULL
  draw_local_global_set(set,png=png,pdf=pdf,data=data,geometry=geometry,onset=entry,
   show_entry=FALSE,effort_geometry=effort_geometry,benchmarks=limits[[paste(family,set,sep='_')]],
   file=paste0('Figure_local_global_',family,'_',set),effort_limit=1,response_limit=1)
 }
 ans<-do.call(rbind,report);write.csv(ans,'10_local_global_results/actor_budget_summaries.csv',row.names=FALSE)
 invisible(ans)
}
if(sys.nframe()==0L)prepare_local_global_figures()
