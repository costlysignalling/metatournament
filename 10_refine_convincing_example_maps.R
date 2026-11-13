# Double the actual talent-grid resolution, reusing all saved grid points.
# This changes the display sampling, not the baseline or population integration.
# source('10_refine_convincing_example_maps.R'); refine_convincing_example_maps()
source('10_convincing_model.R')
convincing_fine_map<-function(d){
 old<-d$map;if(length(old$alpha)>=361L)return(old)
 a<-seq(min(old$alpha),max(old$alpha),length.out=2L*length(old$alpha)-1L)
 n<-length(a);actor<-as.matrix(expand.grid(a,a));out<-matrix(NA_real_,n*n,ncol(old$allocation))
 colnames(out)<-colnames(old$allocation)
 old_index<-as.vector(outer(seq(1L,n,by=2L),seq(0L,n-1L,by=2L)*n,'+'))
 out[old_index,]<-old$allocation
 fresh<-setdiff(seq_len(n*n),old_index)
 message('  Optimizing ',length(fresh),' additional talent combinations')
 out[fresh,]<-convincing_optimize(d$p,actor[fresh,,drop=FALSE])
 keep<-which(d$p$influence_weights>1e-8)
 recipient_abs<-recipient_signed<-numeric(n*n)
 recipient_abs[old_index]<-old$recipient_abs;recipient_signed[old_index]<-old$recipient_signed
 # Batches limit memory while retaining the original recipient optimizer.
 for(begin in seq(1L,length(fresh),by=8192L)){
  ii<-fresh[begin:min(length(fresh),begin+8191L)]
  q<-convincing_compile()$convincing_receivers_cpp(actor[ii,,drop=FALSE],d$p$baseline$responsive$moment,
   d$p$pressure[keep],d$p$beta[1],d$p$beta[2],d$p$phi[1],d$p$phi[2],d$p$B,d$p$N,d$p$nx)
  zero<-which(d$p$pressure[keep]==0)
  if(length(zero))q[,zero]<-out[ii,5]
  change<-q-out[ii,5]
  recipient_abs[ii]<-as.vector(abs(change)%*%d$p$influence_weights[keep])
  recipient_signed[ii]<-as.vector(change%*%d$p$influence_weights[keep])
 }
 stopifnot(all(is.finite(out)),all(recipient_abs>=0),max(recipient_abs)<=d$p$B+1e-7,
  max(abs(rowSums(out[,1:2,drop=FALSE])+abs(out[,3])-d$p$B))<1e-6,
  identical(unname(out[old_index,]),unname(old$allocation)))
 list(alpha=a,allocation=out,recipient_abs=recipient_abs,recipient_signed=recipient_signed)
}
refine_convincing_example_maps<-function(family=c('local','global'),overwrite=FALSE,sets=c('symmetric_beta04','symmetric_beta08','asymmetric')){
 folder<-'10_local_global_results/high_resolution';dir.create(folder,recursive=TRUE,showWarnings=FALSE)
 for(mode in family)for(set in sets){
  data<-if(mode=='local')readRDS(paste0('10_convincing_results/',set,'.rds'))else
   readRDS(paste0('10_comparison_results/popularizing_',set,'.rds'))
  for(j in seq_along(data)){
   file<-file.path(folder,paste0(mode,'_',set,'_',j,'.rds'))
   if(file.exists(file)&&!overwrite)next
   message(format(Sys.time()),' ',mode,' ',set,' panel ',j)
   result<-convincing_fine_map(data[[j]])
   saveRDS(result,file)
   message('  Saved ',length(result$alpha),' x ',length(result$alpha),' map')
  }
 }
 invisible(folder)
}
if(sys.nframe()==0L){
 args<-commandArgs(trailingOnly=TRUE)
 refine_convincing_example_maps(if(length(args))args else c('local','global'))
}


