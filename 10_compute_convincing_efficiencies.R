# Recompute the genuine parameter grid, sharing recipient profiles across
# epsilon=.05,.10,.20. Successful checkpoints are reused automatically.
# For a cosmetic redraw, source('10_convincing_efficiency_figures.R') and call
# convincing_efficiency_figures(); no model calculation is needed then.
source('10_recompute_campaign_parameter_precision.R')
compute_convincing_efficiencies<-function(efficiencies=c(.05,.1,.2),overwrite=FALSE,
 export_figures=FALSE,grid_size=25L,workers=1L){
 stopifnot(all(efficiencies%in%c(.05,.1,.2)))
 stopifnot(length(workers)==1L,is.finite(workers),workers>=1L,workers==as.integer(workers))
 if(workers==1L){
  recompute_campaign_parameter_precision(grid_size=grid_size,overwrite=overwrite)
 }else{
  cluster<-parallel::makePSOCKcluster(workers,outfile='')
  on.exit(parallel::stopCluster(cluster),add=TRUE)
  parallel::clusterCall(cluster,function(root){
   setwd(root)
   source('10_recompute_campaign_parameter_precision.R')
   NULL
  },getwd())
  parallel::clusterApply(cluster,seq_len(workers),
   function(worker,workers,grid_size,overwrite){
    recompute_campaign_parameter_precision(worker,workers,grid_size,overwrite)
   },workers,grid_size,overwrite)
 }
 source('10_campaign_precision_maps.R')
 finalize_campaign_parameter_precision(grid_size)
 source('10_convincing_efficiency_figures.R')
 convincing_efficiency_figures(png=export_figures,pdf=export_figures,rebuild=TRUE)
 invisible(TRUE)
}
if(sys.nframe()==0L)compute_convincing_efficiencies()
