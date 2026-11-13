source('10_compare_convincing_models.R')
source('10_comparison_collectives.R')
comparison_member_predict<-function(p,actor,nx=80L){
 actor<-as.matrix(actor);eq<-p$collective
 base<-convincing_optimize(p$individual_reference,actor,nx)
 n<-nrow(actor);out<-matrix(0,n,10);shift<-absolute_campaign<-numeric(n)
 probability<-as.vector(t(outer(eq$p1,eq$p2)))
 for(k in seq_along(eq$outcomes)){
  a<-eq$outcomes[[k]];v<-ifelse(actor[,1]>actor[,2],a$v[1],a$v[2])
  R<-p$B-abs(v)
  x<-comparison_price_compile()$comparison_member_training_cpp(actor,
   p$B-abs(a$v[1]),p$B-abs(a$v[2]),a$H,p$beta[1],p$beta[2],p$phi[1],p$phi[2],nx)
  y<-R-x;L<-ifelse(a$M>0,a$E/a$M,0)
  U<-L[1]*boredom_w(x,actor[,1],p$beta[1],p$phi[1])+
   L[2]*boredom_w(y,actor[,2],p$beta[2],p$phi[2])
  out[,1:4]<-out[,1:4]+probability[k]*cbind(x,y,v,U)
  shift<-shift+probability[k]*(abs(x-base[,5])+abs(y-(p$B-base[,5])))/2
  absolute_campaign<-absolute_campaign+probability[k]*abs(v)
 }
 out[,5:6]<-base[,5:6];out[,7]<-out[,4]-out[,6];out[,8]<-out[,6]
 colnames(out)<-colnames(base);attr(out,'member_shift')<-shift
 attr(out,'absolute_campaign')<-absolute_campaign;out
}
comparison_member_equal<-function(d,row){
 if(row==1){d$p<-d$p$individual_reference;return(convincing_equal_curves(d,row))}
 p<-d$p;a<-d$map$alpha;n<-length(a);o<-d$map$allocation
 mask<-matrix(o[,1]>=o[,2],n)
 ij<-which(mask[-1,,drop=FALSE]!=mask[-n,,drop=FALSE],arr.ind=TRUE)
 roots<-vector('list',n)
 if(nrow(ij)){
  lo<-a[ij[,1]];hi<-a[ij[,1]+1];yy<-a[ij[,2]];side<-mask[cbind(ij[,1]+1,ij[,2])]
  for(k in 1:12){mid<-(lo+hi)/2;v<-comparison_member_predict(p,cbind(mid,yy))
   above<-v[,1]>=v[,2];move<-above==side;hi[move]<-mid[move];lo[!move]<-mid[!move]}
  for(j in seq_len(n))roots[[j]]<-(lo+hi)[ij[,2]==j]/2
 }
 nr<-max(lengths(roots));curves<-list()
 if(nr)for(k in seq_len(nr))curves[[k]]<-list(x=vapply(roots,
  function(v)if(length(v)>=k)v[k]else NA_real_,0.),y=a)
 if(p$collective$pure)curves<-lapply(curves,function(cl){
  keep<-which(is.finite(cl$x));if(!length(keep))return(cl)
  xy<-cbind(cl$x[keep],cl$y[keep]);pred<-comparison_member_predict(p,xy)
  state<-p$collective$outcomes[[1]];R<-pred[,1]+pred[,2]
  score<-function(x,y)exp(state$H/2)*boredom_w(x,xy[,1],p$beta[1],p$phi[1])+
   exp(-state$H/2)*boredom_w(y,xy[,2],p$beta[2],p$phi[2])
  valid<-score(R/2,R/2)>=score(pred[,1],pred[,2])-1e-6*(1+p$B)
  cl$x[keep[!valid]]<-NA_real_;cl
 })
 curves
}
comparison_member_data<-function(global,eq,alpha){
 reference<-comparison_no_influence(global)
 p<-reference;p$individual_reference<-reference
 p$collective<-eq;p$comparison_family<-'collective_members'
 probability<-as.vector(t(outer(eq$p1,eq$p2)))
 budget<-energy<-numeric(3)
 for(k in seq_along(eq$outcomes)){
  a<-eq$outcomes[[k]];campaign<-mean(abs(a$v))
  training<-a$E-c(sum(pmax(a$v,0)),sum(pmax(-a$v,0)))/2
  budget<-budget+probability[k]*c(training,campaign)
  energy[1:2]<-energy[1:2]+probability[k]*a$E
 }
 p$group_budget<-p$N*budget;p$group_energy<-p$N*energy[1:2]
 p$summary<-c(mean_convincing=budget[3],group_convincing_share=budget[3]/p$B)
 actor<-as.matrix(expand.grid(alpha,alpha));prediction<-comparison_member_predict(p,actor)
 map<-list(alpha=alpha,allocation=prediction,
  recipient_abs=attr(prediction,'member_shift'),recipient_signed=prediction[,1]-prediction[,5],
  expected_absolute_effort=attr(prediction,'absolute_campaign'))
 list(p=p,map=map)
}
comparison_draw_members<-function(set,data,png=FALSE,pdf=FALSE,geometry=NULL){
 # Reuse the approved renderer and polygon construction. Only its model
 # queries change: members follow the coalition policy, not individual BRs.
 env<-new.env(parent=environment(draw_convincing_set))
 env$convincing_optimize<-comparison_member_predict
 env$convincing_boundary_fields<-convincing_boundary_fields
 environment(env$convincing_boundary_fields)<-env
 env$convincing_equal_curves<-comparison_member_equal
 renderer<-draw_convincing_set;environment(renderer)<-env
 renderer(set,png=png,pdf=pdf,data=data,geometry=geometry,
  file=paste0('Figure_comparison_collective_members_',set),effort_limit=1,response_limit=1,
  row_labels=c('Without campaigns','Members after campaigns','Shared campaign effort','Member allocation change'),
  effort_label='Signed mean campaign time / budget',
  response_label='Absolute training change / budget')
}
comparison_reproduce_members<-function(set='asymmetric'){
 data<-list()
 for(j in 1:4){
  eq<-readRDS(paste0('10_comparison_results/collective_equilibrium_',set,'_',j,'.rds'))
  g<-readRDS(paste0('10_comparison_results/popularizing_',set,'_',j,'.rds'))
  d<-comparison_member_data(g$p,eq,g$map$alpha);data[[j]]<-d
  saveRDS(d,paste0('10_comparison_results/collective_members_',set,'_',j,'.rds'))
 }
 saveRDS(data,paste0('10_comparison_results/collective_members_',set,'.rds'))
 comparison_draw_members(set,data,png=TRUE,pdf=TRUE)
 invisible(data)
}
