# Run from the scripts folder. This creates a separate comparison document's
# data and figures; it does not replace Supplement2.Rmd or its original figures.
source('10_comparison_figure_template.R')
source('10_thermodynamic_local_audiences.R')
source('10_popularizing_population.R')
comparison_sets<-c('symmetric_beta04','symmetric_beta08','asymmetric')
comparison_global_prepare<-function(seed,npressure=161L){
 p<-seed;p$N<-1e60;p$baseline$N<-p$N
 p$baseline$recipient_integration<-NULL;p$recipient_integration<-NULL
 p$comparison_family<-'popularizing';p$population_limit<-TRUE
 if(identical(p$beta[1],p$beta[2])&&all(p$phi==1)){
  # This exact reduction avoids two-dimensional quadrature when phi=1.
  p$npressure<-npressure;p$effort<-seq(-p$B,p$B,length.out=npressure)
  p$pressure<-p$efficiency*p$effort
  moments<-lapply(p$pressure,function(delta)
   popularizing_raw_moments_phi1(delta,p$beta[1],p$B))
  p$moments<-t(vapply(moments,function(m)c(m$E[1],m$M,m$E[2]),numeric(4)))
  mid<-(npressure+1)/2
  p$baseline$responsive$moment<-p$moments[mid,1:3]
  p$baseline$responsive$mapped<-p$moments[mid,1:3]
  p$baseline$responsive$residual<-0;p$zero_error<-0
  r<-log(p$density$talent[,1]/p$density$talent[,2])
  p$receiver_population<-vapply(p$pressure,function(delta)
   popularizing_recipient_phi1(r,delta,p$beta[1],p$B),numeric(length(r)))
  p$recipient_integration<-'analytic-phi1'
  p<-convincing_refresh(p)
 }else p<-convincing_refine_recipients(p,npressure)
 p
}
comparison_no_influence<-function(global){
 p<-global;mid<-(length(p$effort)+1)/2
 p$moments<-p$moments[rep(mid,3),,drop=FALSE]
 p$receiver_population<-p$receiver_population[,rep(mid,3),drop=FALSE]
 p$effort<-c(-p$B,0,p$B);p$pressure<-c(0,0,0);p$npressure<-3L
 p$comparison_family<-'global_finite_audience'
 p$adaptive_integration<-NULL
 p<-convincing_refresh(p)
 stopifnot(max(abs(p$actors[,3]))<1e-8)
 p
}
comparison_reproduce_set<-function(set='asymmetric',nalpha=NULL,draw=FALSE){
 dir.create('10_comparison_results',showWarnings=FALSE)
 source_data<-readRDS(file.path('10_convincing_results',paste0(set,'.rds')))
 na<-if(is.null(nalpha))if(set=='asymmetric')361L else 181L else nalpha
 family_data<-setNames(lapply(1:3,function(i)list()),
  c('partial_local','popularizing','global_finite_audience'))
 for(j in seq_along(source_data)){
  seed<-source_data[[j]]$p
  for(family in names(family_data)){
   file<-file.path('10_comparison_results',paste0(family,'_',set,'_',j,'.rds'))
   if(file.exists(file)){
    d<-readRDS(file)
    stopifnot(isTRUE(all.equal(d$p$beta,seed$beta)),isTRUE(all.equal(d$p$phi,seed$phi)),
     d$p$B==seed$B,d$p$efficiency==seed$efficiency,length(d$map$alpha)==na)
   }else{
    message(family,' ',set,' ',j)
    if(family=='partial_local'){
     p<-thermo_local_audience(seed,audience=15);p$comparison_family<-family
     map<-thermo_local_map(p,seq(0,2.5,length.out=na))
    }else if(family=='popularizing'){
     p<-comparison_global_prepare(seed)
     map<-convincing_map(p,seq(0,2.5,length.out=na))
    }else{
     global<-family_data$popularizing[[j]]$p
     p<-comparison_no_influence(global)
     map<-convincing_map(p,seq(0,2.5,length.out=na))
    }
    d<-list(p=p,map=map);saveRDS(d,file)
   }
   family_data[[family]][[j]]<-d
  }
 }
 for(family in names(family_data)){
  data<-family_data[[family]]
  saveRDS(data,file.path('10_comparison_results',paste0(family,'_',set,'.rds')))
  if(draw)draw_convincing_set(set,png=TRUE,pdf=TRUE,data=data,
   file=paste0('Figure_comparison_',family,'_',set))
 }
 invisible(family_data)
}

