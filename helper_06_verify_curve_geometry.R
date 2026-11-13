planes<-commandArgs(trailingOnly=TRUE);if(!length(planes))planes<-c('beta','phi')
for(plane in planes){
 f<-paste0('06_coexistence_',plane,'_geometry.rds');d<-readRDS(f)
 if(is.null(d$brackets)){
  records<-lapply(list.files(file.path('coexistence_curve_checkpoints',plane),full.names=TRUE,pattern='[.]rds$'),function(f)tryCatch(readRDS(f),error=function(e)NULL))
  records<-Filter(function(r)!is.null(r)&&identical(r$signature,d$signature)&&isTRUE(r$complete),records)
  d$brackets<-lapply(records,function(r)list(lo=r$lo,hi=r$hi,flo=r$flo,fhi=r$fhi));saveRDS(d,f)
 }
 for(b in d$brackets){
  stopifnot(all(is.finite(b$flo)),all(is.finite(b$fhi)),max(b$hi-b$lo)<=2.01*d$edge_tolerance,
   all(b$flo*b$fhi<=0|pmin(abs(b$flo),abs(b$fhi))<1e-10))
 }
 input<-readRDS(paste0('06_coexistence_',plane,'_inputs.rds'));area<-diff(range(input$axis))^2
 for(mode in c('fixed','responsive'))for(i in 1:3)for(j in 1:3){
  g<-d$geometry[[mode]][[i,j]]$dissolved
  stopifnot(all(sf::st_is_valid(g)),abs(sum(as.numeric(sf::st_area(g)))-area)<1e-7,
   all(g$rid>=1&g$rid<=15))
 }
 cat(plane,':',sum(vapply(d$brackets,function(b)length(b$lo),1L)),'actual-model root brackets verified; polygons valid and cover the full square\n')
}
