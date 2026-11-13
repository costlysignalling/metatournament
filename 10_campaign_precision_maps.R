# Read the genuine parameter-grid calculations; no smoothing across parameters.
source('10_plot_convincing_parameter_levels.R')
assemble_campaign_precision_fields<-function(efficiency=.1,save_output=FALSE){
 meta<-readRDS('10_campaign_precision_results/complete.rds')
 a<-meta$axis;phis<-c(.66,1,1.5);n<-length(a)
 panels<-setNames(lapply(1:2,function(k)matrix(vector('list',9),3,3)),c('local','global'))
 for(mode in names(panels))for(c in 1:3)for(r in 1:3){
  c1<-c2<-matrix(NA_real_,n,n)
  for(i in seq_len(n))for(j in seq_len(n)){
   cc<-c;rr<-r;ii<-i;jj<-j;swap<-FALSE
   if(cc<rr){z<-cc;cc<-rr;rr<-z;z<-ii;ii<-jj;jj<-z;swap<-!swap}
   if(cc==rr&&jj<ii){z<-ii;ii<-jj;jj<-z;swap<-!swap}
   file<-sprintf('10_campaign_precision_results/cells/%s_p%.2f_%.2f_b%.6f_%.6f.rds',mode,phis[cc],phis[rr],a[ii],a[jj])
   cell<-readRDS(file);v<-unname(cell$values[[sprintf('%.2f',efficiency)]][c('campaign1','campaign2')])
   if(swap)v<-rev(v);if(c==r&&i==j)v[]<-mean(v)
   c1[i,j]<-v[1];c2[i,j]<-v[2]
  }
  panels[[mode]][[c,r]]<-list(axis=a,C1=c1,C2=c2,actual_grid=n)
 }
 all<-unlist(lapply(panels,function(m)lapply(m,function(p)as.vector(p$C1+p$C2))))
 d<-list(axis=a,phis=phis,panels=panels,range=range(all),B=50,efficiency=efficiency,local_N=30,
  schema_version=3L,uniform_grid=TRUE,centres_complete=TRUE,
  display_max=ceiling(max(all)/.05)*.05,precision=meta)
 if(save_output){root<-if(efficiency==.1)'10_local_global_results'else sprintf('10_efficiency_%.2f_results',efficiency)
  dir.create(root,recursive=TRUE,showWarnings=FALSE)
  saveRDS(d,file.path(root,'parameter_levels.rds'))}
 d
}
finalize_campaign_parameter_precision<-function(grid_size=25L){
 source('10_recompute_campaign_parameter_precision.R')
 tasks<-campaign_precision_tasks(grid_size);files<-vapply(seq_len(nrow(tasks)),function(k)campaign_precision_filename(tasks[k,]),'')
 if(!all(file.exists(files)))stop('Missing ',sum(!file.exists(files)),' parameter calculations.')
 for(file in files){cell<-readRDS(file)
  stopifnot(is.null(cell$error),length(cell$values)==3L)
  for(v in cell$values)stopifnot(all(is.finite(v)),min(v)>=-1e-10,
   abs(v['campaign1']+v['campaign2']-v['campaign'])<1e-9,
   abs(sum(v[c('campaign','training1','training2')])-1)<1e-7)
 }
 meta<-list(axis=seq(.1,.9,length.out=grid_size),grid_size=grid_size,independent_cells=nrow(tasks),
  efficiencies=c(.05,.1,.2),recipient_order=32768L,focal_order=2048L,pressure_knots=campaign_precision_knots(),
  version=2L,method='dense recipient integration with adaptive rare tails and monotone cubic pressure interpolation')
 saveRDS(meta,'10_campaign_precision_results/complete.rds')
 invisible(meta)
}
