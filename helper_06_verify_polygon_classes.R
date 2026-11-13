planes<-commandArgs(trailingOnly=TRUE);if(!length(planes))planes<-c('beta','phi')
fields<-c('S2max','S2div','S1div','S1max');levels<-c(0,1e-4,1e-4,0)
for(plane in planes){
 d<-readRDS(paste0('06_coexistence_',plane,'_inputs.rds'));g<-readRDS(paste0('06_coexistence_',plane,'_geometry.rds'))
 gr<-expand.grid(x=d$axis,y=d$axis);total<-badtotal<-0L
 for(mode in c('fixed','responsive'))for(i in 1:3)for(j in 1:3){
  p<-d$panels[[i,j]];good<-as.vector(p$direct)&gr$x>min(d$axis)&gr$x<max(d$axis)&gr$y>min(d$axis)&gr$y<max(d$axis)
  ids<-which(good);pts<-sf::st_as_sf(gr[ids,,drop=FALSE],coords=c('x','y'))
  inside<-sf::st_within(pts,g$geometry[[mode]][[i,j]]$dissolved);single<-lengths(inside)==1
  expected<-Reduce(`+`,lapply(1:4,function(k)as.integer(as.vector(p[[mode]][[fields[k]]])[ids]>levels[k])*2^(k-1)))
  actual<-vapply(inside[single],function(z)g$geometry[[mode]][[i,j]]$dissolved$rid[z],1)
  bad<-sum(actual!=expected[single]);total<-total+sum(single);badtotal<-badtotal+bad
  if(bad)cat(plane,mode,i,j,':',bad,'differences\n')
 }
 cat(plane,':',badtotal,'classification differences among',total,'direct-node polygon checks\n')
 stopifnot(badtotal==0)
}
