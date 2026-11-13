# Assemble and redraw the three efficiency comparisons from saved calculations.
# Run from the scripts folder. Calculation uses R/C++ and the same lognormal
# talents (median=1, SD=0.5), B=50, and N=30 for the local condition.
# Each map is a first intervention at the non-convincing equilibrium.
source('10_plot_convincing_parameter_levels.R')
assemble_efficiency_fields<-function(efficiency=.1,require_complete=TRUE,save_output=FALSE){
 if(file.exists('10_campaign_precision_results/complete.rds')){
  if(!exists('assemble_campaign_precision_fields',mode='function'))source('10_campaign_precision_maps.R')
  return(assemble_campaign_precision_fields(efficiency,save_output=save_output))
 }
 root<-if(efficiency==.1)'10_local_global_results' else sprintf('10_efficiency_%.2f_results',efficiency)
 a<-seq(.1,.9,length.out=if(efficiency==.1)13L else 9L);n<-length(a);phis<-c(.66,1,1.5)
 panels<-setNames(lapply(1:2,function(k)matrix(vector('list',9),3,3)),c('local','global'))
 get<-function(family,c,r,i,j,axis,kind){
  swap<-FALSE
  if(c<r){z<-c;c<-r;r<-z;z<-i;i<-j;j<-z;swap<-!swap}
  if(c==r&&j<i){z<-i;i<-j;j<-z;swap<-!swap}
  file<-sprintf('%s/%s/%s_p%.2f_%.2f_b%.6f_%.6f.rds',root,kind,family,phis[c],phis[r],axis[i],axis[j])
  if(!file.exists(file)){if(require_complete)stop('Missing calculation: ',file);return(c(NA_real_,NA_real_))}
  d<-readRDS(file);if(is.list(d$result))stop('Unresolved calculation: ',file,': ',d$result$error)
  v<-unname(d$result[c('campaign1','campaign2')]);stopifnot(all(is.finite(v)),all(v>=-1e-12),sum(v)<=1+1e-8)
  if(swap)v<-rev(v);if(c==r&&i==j)v[]<-mean(v);v
 }
 for(f in names(panels))for(c in 1:3)for(r in 1:3){
  C1<-C2<-matrix(NA_real_,n,n);m1<-m2<-matrix(NA_real_,n-1,n-1);mid<-(head(a,-1)+tail(a,-1))/2
  for(i in 1:n)for(j in 1:n){v<-get(f,c,r,i,j,a,'parameter_cells');C1[i,j]<-v[1];C2[i,j]<-v[2]}
  for(i in 1:(n-1))for(j in 1:(n-1)){v<-get(f,c,r,i,j,mid,'parameter_center_cells');m1[i,j]<-v[1];m2[i,j]<-v[2]}
  panels[[f]][[c,r]]<-list(axis=a,C1=C1,C2=C2,center_C1=m1,center_C2=m2,actual_grid=n)
 }
 if(efficiency==.1){z<-readRDS(file.path(root,'linear_beta25.rds'));stopifnot(z$complete)
  panels$global[[2,2]]<-list(axis=z$axis,C1=z$C1,C2=z$C2,actual_grid=length(z$axis))}
 values<-unlist(lapply(panels,function(mm)lapply(mm,function(p)c(as.vector(p$C1+p$C2),as.vector(p$center_C1+p$center_C2)))))
 d<-list(axis=a,phis=phis,panels=panels,range=range(values,na.rm=TRUE),B=50,efficiency=efficiency,local_N=30,schema_version=2L,centres_complete=TRUE)
 d$display_max<-ceiling(max(values,na.rm=TRUE)/.05)*.05
 if(save_output)saveRDS(d,file.path(root,'parameter_levels.rds'));d
}
convincing_efficiency_figures<-function(png=TRUE,pdf=TRUE,rebuild=FALSE){
 # Cosmetic redraws read three compact aggregate files. Set rebuild=TRUE
 # after changing raw calculations or the parameter grid.
 eps<-c(.05,.1,.2)
 datasets<-lapply(eps,function(e){
  root<-if(e==.1)'10_local_global_results' else sprintf('10_efficiency_%.2f_results',e)
  file<-file.path(root,'parameter_levels.rds')
  if(!rebuild&&file.exists(file)){
   cached<-readRDS(file)
   if(cached$schema_version %in% c(2L,3L)&&isTRUE(cached$centres_complete))return(cached)
  }
  assemble_efficiency_fields(e)
 })
 common_max<-max(vapply(datasets,function(d)d$display_max,0.))
 for(k in seq_along(eps)){
  d<-datasets[[k]];d$display_max<-common_max
  root<-if(eps[k]==.1)'10_local_global_results' else sprintf('10_efficiency_%.2f_results',eps[k])
  dir.create(root,recursive=TRUE,showWarnings=FALSE)
  saveRDS(d,file.path(root,'parameter_levels.rds'))
  suffix<-if(eps[k]==.1)'' else sprintf('_epsilon_%.2f',eps[k])
  plot_campaign_levels(png,pdf,data=d,file=paste0('Figure_convincing_local_global_parameter_levels',suffix))
  plot_campaign_direction(png,pdf,data=d,file=paste0('Figure_convincing_campaign_direction',suffix))
 }
 invisible(datasets)
}
if(sys.nframe()==0L)convincing_efficiency_figures()
