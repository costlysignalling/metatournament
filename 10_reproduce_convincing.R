# From this folder: source('10_reproduce_convincing.R'); reproduce_convincing()
# Computation is separate from redrawing. Default figures: N=30, B=50,
# effectiveness=.1, lognormal talents with median=1, SD=.5.
source('10_convincing_figures.R')
source('05_asymmetric_examples.R')
reproduce_convincing<-function(n=8192L,npressure=81L,nalpha=NULL,efficiency=.1){
 dir.create('10_convincing_results',showWarnings=FALSE)
 sets<-list(symmetric_beta04=lapply(c(.66,1,1.25,1.5),function(phi)list(beta=c(.4,.4),phi=rep(phi,2))),
  symmetric_beta08=lapply(c(.66,1,1.25,1.5),function(phi)list(beta=c(.8,.8),phi=rep(phi,2))),
  asymmetric=lapply(asymmetric_examples,function(s)list(beta=c(s$beta1,s$beta2),phi=c(s$phi1,s$phi2))))
 summary<-list()
 for(set in names(sets)){
  data<-list()
  for(j in seq_along(sets[[set]])){
   s<-sets[[set]][[j]];d<-NULL
   adaptive<-identical(set,'asymmetric')
   np<-if(adaptive)max(161L,npressure)else npressure
   na<-if(is.null(nalpha))if(adaptive)361L else 181L else nalpha
   cache<-file.path('10_convincing_results',paste0(set,'_',j,'.rds'))
   if(file.exists(cache)){d<-readRDS(cache)
    if(!identical(c(d$p$n,d$p$npressure,length(d$map$alpha),d$p$efficiency),c(n,np,na,efficiency)))stop('Cached settings differ; use a new output folder or remove the specific saved RDS first.')
    if(!isTRUE(all.equal(d$p$beta,s$beta))||!isTRUE(all.equal(d$p$phi,s$phi))||!d$p$version%in%c('energy-convincing-v8','adaptive-recipients-v9'))d<-NULL
    if(!is.null(d)&&(isTRUE(d$p$aligned)||is.null(d$p$group_budget)||!identical(d$p$version,if(adaptive)'adaptive-recipients-v9'else'energy-convincing-v8'))){
     message('Update unrestricted choices ',set,' ',j)
     d$p$aligned<-FALSE;d$p<-if(adaptive)convincing_refine_recipients(d$p,np)else convincing_refresh(d$p)
     d$map<-convincing_map(d$p,seq(0,2.5,length.out=na));saveRDS(d,cache)
    }
   }
   if(is.null(d)){
    message('Calculate ',set,' ',j,'/',length(sets[[set]]))
    s<-sets[[set]][[j]];p<-convincing_prepare(s$beta,s$phi,n=n,npressure=np,efficiency=efficiency,adaptive=adaptive)
    map<-convincing_map(p,seq(0,2.5,length.out=na));d<-list(p=p,map=map);saveRDS(d,cache)
   }
   data[[j]]<-d;summary[[paste(set,j,sep='_')]]<-c(d$p$summary,zero_error=d$p$zero_error,baseline_residual=d$p$baseline$responsive$residual)
  }
  saveRDS(data,file.path('10_convincing_results',paste0(set,'.rds')))
  draw_convincing_set(set,png=TRUE,pdf=TRUE,data=data)
 }
 tab<-do.call(rbind,summary);write.csv(tab,'10_convincing_results/summary.csv')
 invisible(tab)
}
if(sys.nframe()==0L){reproduce_convincing()}
