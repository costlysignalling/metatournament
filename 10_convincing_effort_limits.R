source('10_convincing_model.R')
# Scalar specialist objectives use the prepared recipient response, not a
# regression fitted to the displayed colours.
campaign_scalar_max<-function(fun,B){
 grid<-seq(0,B,length.out=401);value<-fun(grid);best<-which.max(value)
 ans<-c(effort=grid[best],value=value[best])
 peaks<-which(value[2:400]>=value[1:399]&value[2:400]>=value[3:401])+1L
 for(k in peaks){z<-optimize(fun,grid[c(k-1,k+1)],maximum=TRUE,tol=1e-8)
  if(z$objective>ans[2])ans<-c(effort=z$maximum,value=z$objective)}
 if(ans[2]<=value[1]+1e-9*(1+abs(value[1])))ans<-c(effort=0,value=value[1])
 ans
}
campaign_specialist_limits<-function(p){
 bg<-convincing_background_function(p);B<-p$B;local<-p$N<1e30
 out<-list()
 for(t in 1:2){same<-if(t==1)1 else -1
  objective<-function(y,d,low=FALSE){m<-bg(d*y);E<-m[,t];M<-m[,t+2]
   R<-B-y;g<-expm1(p$beta[t]*log1p(R))
   if(local){numerator<-(p$N-1)*E+R+if(d==same)y else 0
    numerator/((p$N-1)*M)*g^if(low)1 else p$phi[t]
   }else E/M*g^if(low)1 else p$phi[t]
  }
  reference<-campaign_scalar_max(function(y){
   m<-bg(same*y);R<-B-y;claim<-boredom_w(R,rep(1,length(y)),p$beta[t],p$phi[t])
   if(local)((p$N-1)*m[,t]+B)*claim/((p$N-1)*m[,t+2]+claim)
   else m[,t]/m[,t+2]*claim
  },B)
  normal<-if(local)reference else campaign_scalar_max(function(y)objective(y,same,FALSE),B)
  cross<-campaign_scalar_max(function(y)objective(y,-same,TRUE),B)
  lowsame<-campaign_scalar_max(function(y)objective(y,same,TRUE),B)
  out[[t]]<-list(tournament=t,normal_signed=same*normal[1]/B,
   low_cross_signed=-same*cross[1]/B,low_cross_preferred=cross[2]>lowsame[2]+1e-8,
   low_same_signed=same*lowsame[1]/B,normal=normal,cross=cross,low_same=lowsame,reference_signed=same*reference[1]/B,reference=reference,
   normal_kind=if(local)'median-talent specialist benchmark'else if(p$phi[t]==1)'exact specialist plateau'else 'large-talent specialist limit')
 }
 out
}
campaign_limit_summary<-function(){
 results<-list();rows<-list()
 for(family in c('local','global'))for(set in c('symmetric_beta04','symmetric_beta08','asymmetric')){
  file<-if(family=='local')paste0('10_convincing_results/',set,'.rds')else paste0('10_comparison_results/popularizing_',set,'.rds')
  data<-readRDS(file);z<-lapply(data,function(d)campaign_specialist_limits(d$p))
  results[[paste(family,set,sep='_')]]<-z
  for(j in seq_along(z))for(t in 1:2)rows[[length(rows)+1]]<-data.frame(family=family,set=set,panel=j,tournament=t,
   solid_marker=unname(z[[j]][[t]]$normal_signed),median_reference=unname(z[[j]][[t]]$reference_signed),low_cross_limit=unname(z[[j]][[t]]$low_cross_signed),
   cross_preferred=z[[j]][[t]]$low_cross_preferred,kind=z[[j]][[t]]$normal_kind)
 }
 saveRDS(results,'10_local_global_results/effort_limits.rds')
 write.csv(do.call(rbind,rows),'10_local_global_results/effort_limits.csv',row.names=FALSE)
 invisible(results)
}
