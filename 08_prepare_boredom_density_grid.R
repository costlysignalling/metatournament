# Deterministic expected-density integration; no random classroom is selected.
# Splitting phi2 rows into separate R sessions is optional and changes no data.
source('08_boredom_model.R')
prepare_boredom_density_grid<-function(rows=1:3,B=50,compulsory=.5,N=30,nbeta=21L,n=2048L){
 dir.create('08_boredom_results',showWarnings=FALSE)
 den<-boredom_qmc_density(n);beta<-seq(.01,.99,length.out=nbeta);phi<-c(.66,1,1.5)
 tag<-paste0('c_',format(compulsory,trim=TRUE),'_N_',N)
 signature<-list(B=B,compulsory=compulsory,N=N,nbeta=nbeta,n=n,version=1L,core=unname(tools::md5sum('08_boredom_core.cpp')))
 for(j in rows)for(i in 1:3){
  file<-file.path('08_boredom_results',paste0('density_panel_',i,'_',j,'_',tag,'.rds'))
  if(file.exists(file)&&identical(attr(readRDS(file),'signature'),signature))next
  z<-matrix(NA_real_,nbeta,nbeta);p<-list(phi1=phi[i],phi2=phi[j],first=list(s0=z,s1=z,s2=z),responsive=list(s0=z,s1=z,s2=z),status=z,residual=z,order=z)
  st<-NULL
  for(y in seq_along(beta)){
   for(x in seq_along(beta)){
    fit<-boredom_population(c(beta[x],beta[y]),c(phi[i],phi[j]),B,compulsory,N,den,80L,1e-5,st,maxit=70L)
    if(!fit$responsive$converged)for(nhigh in c(8192L,32768L)){
     fit<-boredom_population(c(beta[x],beta[y]),c(phi[i],phi[j]),B,compulsory,N,boredom_qmc_density(nhigh),100L,1e-5,fit$responsive$moment,maxit=100L)
     if(fit$responsive$converged)break
    }
    st<-if(fit$responsive$converged)fit$responsive$moment else NULL
    for(mode in c('first','responsive')){m<-fit[[mode]]$mapped;p[[mode]]$s0[x,y]<-m[4]/B;p[[mode]]$s1[x,y]<-m[1]/B;p[[mode]]$s2[x,y]<-1-(m[4]+m[1])/B}
    p$status[x,y]<-as.integer(fit$responsive$converged);p$residual[x,y]<-fit$responsive$residual;p$order[x,y]<-fit$order
   }
   if(y%%5==0){cat('phi panels',i,j,'row',y,'/',nbeta,'\n');flush.console()}
  }
  attr(p,'signature')<-signature;saveRDS(p,file);cat('Panel',i,j,'ready; unresolved',sum(p$status!=1),'\n');flush.console()
 }
 files<-outer(1:3,1:3,Vectorize(function(i,j)file.path('08_boredom_results',paste0('density_panel_',i,'_',j,'_',tag,'.rds'))))
 if(all(file.exists(files))){panels<-lapply(as.vector(files),readRDS);dim(panels)<-c(3,3)
  saveRDS(list(B=B,compulsory=compulsory,N=N,beta=beta,phi=phi,panels=panels,method='expected-density QMC'),file.path('08_boredom_results',paste0('energy_',tag,'.rds')))
 }
 invisible(NULL)
}
if(sys.nframe()==0L)prepare_boredom_density_grid()
