source('09_harvesting_model.R')
prepare_harvesting_phispace<-function(nphi=41L,betas=c(.1,.5,.9),order=48L,inner=16L,overwrite=FALSE){
 dir.create('09_harvesting_results',showWarnings=FALSE);axis<-seq(.1,3.1,length.out=nphi)
 radius<-sqrt(qchisq(.99,2,ncp=2*boredom_slog()^2))
 for(b in betas){
  file<-paste0('09_harvesting_results/phi_beta_',b,'.rds')
  if(file.exists(file)&&!overwrite)next
  z<-matrix(NA_real_,nphi,nphi);p<-list(beta2=b,axis=axis,ratio_first=z,ratio_responsive=z,E_first=z,E_responsive=z,residual=z,
   presence_first=replicate(4,z,simplify=FALSE),presence_responsive=replicate(4,z,simplify=FALSE))
  warm<-NULL
  for(j in seq_along(axis)){
   warm<-NULL
   for(i in seq_along(axis)){
    fit<-harvest_population(b,axis[i],axis[j],order=order,inner=inner,start=warm);warm<-fit$responsive_ratio
    p$ratio_first[i,j]<-fit$fixed_ratio;p$ratio_responsive[i,j]<-warm
    p$E_first[i,j]<-fit$first[1];p$E_responsive[i,j]<-fit$responsive[1];p$residual[i,j]<-fit$residual
    en<-harvest_engine(b,axis[i],axis[j],order=12L,inner=6L)
    for(mode in c('first','responsive')){
     score<-harvest_core()$harvesting_presence_cpp(en,if(mode=='first')fit$fixed_ratio else warm,radius)
     for(k in 1:4)p[[paste0('presence_',mode)]][[k]][i,j]<-score[k]
    }
   }
   cat('Harvesting phi beta',b,'row',j,'/',nphi,'\n');flush.console()
  }
  p$presence_method<-'signed-support-v2';saveRDS(p,file)
 }
 invisible(NULL)
}
refresh_harvesting_presence<-function(){
 radius<-sqrt(qchisq(.99,2,ncp=2*boredom_slog()^2))
 for(b in c(.1,.5,.9)){
  file<-paste0('09_harvesting_results/phi_beta_',b,'.rds');p<-readRDS(file)
  if(identical(p$presence_method,'signed-support-v2'))next
  for(j in seq_along(p$axis))for(i in seq_along(p$axis)){
   en<-harvest_engine(b,p$axis[i],p$axis[j],order=8L,inner=4L)
   for(mode in c('first','responsive')){
    sc<-harvest_core()$harvesting_presence_cpp(en,p[[paste0('ratio_',mode)]][i,j],radius)
    for(k in 1:4)p[[paste0('presence_',mode)]][[k]][i,j]<-sc[k]
   }
  }
  p$presence_method<-'signed-support-v2';saveRDS(p,file);message('Refined signed presence: beta2=',b)
 }
}
if(sys.nframe()==0L)prepare_harvesting_phispace()

