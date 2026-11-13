# Shared, approved figure design with common effort/response scales.
source('10_comparison_member_figures.R')
comparison_geometry<-function(data,before=NULL,member=FALSE){
 if(member){
  env<-new.env(parent=environment(convincing_boundary_fields))
  env$convincing_optimize<-comparison_member_predict
  boundary<-convincing_boundary_fields;environment(boundary)<-env
  equal<-comparison_member_equal
 }else{boundary<-convincing_boundary_fields;equal<-convincing_equal_curves}
 list(alpha=lapply(data,function(d)d$map$alpha),
  fields=lapply(seq_along(data),function(j)list(
   if(is.null(before))boundary(data[[j]],1)else before$fields[[j]][[1]],
   boundary(data[[j]],2))),
  equal=lapply(seq_along(data),function(j)list(
   if(is.null(before))equal(data[[j]],1)else before$equal[[j]][[1]],
   equal(data[[j]],2))))
}
comparison_render_set<-function(set='asymmetric',include_members=TRUE){
 original<-readRDS(paste0('10_convincing_results/',set,'.rds'))
 local_geometry<-readRDS(paste0('10_convincing_results/',set,'_geometry.rds'))
 families<-c('local_full','partial_local','global_finite_audience','popularizing','collective_members')
 if(!include_members)families<-families[families!='collective_members']
 zero_geometry<-NULL
 for(family in families){
  message('Render ',family,' ',set)
  if(family=='local_full')data<-original else if(family=='collective_members'){
   path<-paste0('10_comparison_results/collective_members_',set,'.rds')
   equilibria<-lapply(1:4,function(j)readRDS(paste0(
    '10_comparison_results/collective_equilibrium_',set,'_',j,'.rds')))
   stale<-TRUE
   if(file.exists(path)){
    data<-readRDS(path)
    stale<-!all(vapply(1:4,function(j)identical(data[[j]]$p$collective,equilibria[[j]]),TRUE))
   }
   if(stale){
    geometry_file<-paste0('10_comparison_results/geometry_collective_members_',set,'.rds')
    if(file.exists(geometry_file))file.remove(geometry_file)
    data<-lapply(1:4,function(j){
     eq<-equilibria[[j]]
     g<-readRDS(paste0('10_comparison_results/popularizing_',set,'_',j,'.rds'))
     d<-comparison_member_data(g$p,eq,g$map$alpha)
     saveRDS(d,paste0('10_comparison_results/collective_members_',set,'_',j,'.rds'));d
    })
    saveRDS(data,path)
   }
  }else data<-readRDS(paste0('10_comparison_results/',family,'_',set,'.rds'))
  file<-paste0('10_comparison_results/geometry_',family,'_',set,'.rds')
  if(file.exists(file))geometry<-readRDS(file)else{
   if(family=='local_full')geometry<-local_geometry else if(family=='partial_local')
    geometry<-comparison_geometry(data,before=local_geometry)else if(family=='global_finite_audience'){
     geometry<-comparison_geometry(data)
     geometry$fields<-lapply(geometry$fields,function(v)list(v[[1]],v[[1]]))
     geometry$equal<-lapply(geometry$equal,function(v)list(v[[1]],v[[1]]))
   }else geometry<-comparison_geometry(data,before=zero_geometry,member=family=='collective_members')
   saveRDS(geometry,file)
  }
  if(family=='global_finite_audience')zero_geometry<-geometry
  if(family=='collective_members')comparison_draw_members(set,data,png=TRUE,pdf=TRUE,geometry=geometry)else
   draw_convincing_set(set,png=TRUE,pdf=TRUE,data=data,geometry=geometry,
    file=paste0('Figure_comparison_',family,'_',set),effort_limit=1,response_limit=1)
 }
 invisible(TRUE)
}
comparison_redraw_all<-function(){
 for(set in comparison_sets)comparison_render_set(set)
 invisible(TRUE)
}
