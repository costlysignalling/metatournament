rows<-list()
for(family in c('local','global'))for(set in c('symmetric_beta04','symmetric_beta08','asymmetric')){
 file<-if(family=='local')paste0('10_convincing_results/',set,'.rds')else paste0('10_comparison_results/popularizing_',set,'.rds')
 data<-readRDS(file)
 for(j in seq_along(data)){
  p<-data[[j]]$p;o<-p$actors;w<-p$density$weights
  q0<-o[,5]/p$B;q<-o[,1]/pmax(1e-10,o[,1]+o[,2]);yes<-abs(o[,3])>1e-4*p$B
  farther<-abs(q-.5)>abs(q0-.5)+1e-6
  rows[[length(rows)+1]]<-data.frame(family=family,set=set,panel=j,
   convince=sum(w*yes),farther_given_convincing=sum(w*yes*farther)/sum(w*yes),
   closer_given_convincing=sum(w*yes*(abs(q-.5)<abs(q0-.5)-1e-6))/sum(w*yes),
   max_change_nonconvincers=if(any(!yes))max(abs(q[!yes]-q0[!yes]))else NA_real_)
 }
}
write.csv(do.call(rbind,rows),'10_local_global_results/polarization_checks.csv',row.names=FALSE)
print(do.call(rbind,rows))
