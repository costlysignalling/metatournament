# Infinite population of finite, independent local contest groups.
# N in the original solver remains the LOCAL group size. Each group has one
# potential persuader; audience is the number of local competitors reached.
# The total number of groups is unlimited and does not enter the solver.
source('10_convincing_model.R')
thermo_local_audience<-function(reference,audience=reference$N-1){
 stopifnot(is.list(reference),ncol(reference$moments)==4L,
  length(audience)==1,is.finite(audience),audience>=0,
  audience<=reference$N-1,audience==as.integer(audience))
 if(!is.null(reference$audience_size)&&reference$audience_size!=reference$N-1)
  stop('Use a full-audience reference as the first argument.')
 p<-reference;n<-p$N;share<-audience/(n-1);mid<-(length(p$effort)+1L)/2L
 base<-c(p$baseline$responsive$moment,p$B-p$baseline$responsive$moment[1])
 # Mix CLAIMS separately from investments. Applying w() to an averaged
 # investment would incorrectly replace random exposure with partial training.
 p$moments<-share*reference$moments+(1-share)*matrix(base,nrow(reference$moments),4,byrow=TRUE)
 q0<-reference$receiver_population[,mid]
 # This stores expected recipient allocations. For a two-point exposure
 # mixture, |E[X]-X0| equals E[|X-X0|], since the unexposed change is zero.
 p$receiver_population<-share*reference$receiver_population+(1-share)*q0
 p$audience_size<-audience;p$audience_fraction<-share
 p<-convincing_refresh(p)
 p
}
thermo_local_map<-function(p,alpha=seq(0,2.5,length.out=181L)){
 ans<-convincing_map(p,alpha)
 share<-if(is.null(p$audience_fraction))1 else p$audience_fraction
 ans$recipient_abs<-share*ans$recipient_abs;ans$recipient_signed<-share*ans$recipient_signed
 ans
}
thermo_local_summary<-function(p){
 n<-p$N;B<-p$B;s<-p$summary
 data.frame(local_group=n,audience=if(is.null(p$audience_size))n-1 else p$audience_size,
  potential_persuaders_pct=100*unname(s['persuader_fraction']),
  all_actors_convincing_pct=100*unname(s['persuader_fraction'])/n,
  selected_actor_convincing_pct=100*unname(s['mean_convincing'])/B,
  population_budget_convincing_pct=100*unname(s['mean_convincing'])/(n*B),
  tournament1_budget_pct=100*unname(p$group_energy[1])/(n*B),
  tournament2_budget_pct=100*unname(p$group_energy[2])/(n*B))
}
thermo_local_example<-function(audience=15,actor=c(2,.5),
 reference_file='10_convincing_results/symmetric_beta04_2.rds'){
 reference<-readRDS(reference_file)$p
 p<-thermo_local_audience(reference,audience)
 print(convincing_optimize(p,matrix(actor,nrow=1)));print(thermo_local_summary(p))
 invisible(p)
}
plot_thermo_local_audience<-function(data,png=FALSE,pdf=FALSE,
 file='Figure_convincing_thermodynamic_local_audience'){
 table<-data$summary
 draw<-function(){
  par(mfrow=c(1,2),mar=c(3.4,4,2.1,.65),mgp=c(2,.65,0),tcl=-.25,cex=1.1)
  plot(table$audience,table$all_actors_convincing_pct,type='b',pch=16,lwd=1.6,
   xlab='Local competitors reached',ylab='Actors convincing (%)',ylim=c(0,max(.1,table$all_actors_convincing_pct)*1.05))
  mtext('A',3,line=.7,adj=0,font=2,cex=1.15)
  plot(table$audience,table$population_budget_convincing_pct,type='b',pch=16,lwd=1.6,
   xlab='Local competitors reached',ylab='Budget spent convincing (%)',ylim=c(0,max(.01,table$population_budget_convincing_pct)*1.05))
  mtext('B',3,line=.7,adj=0,font=2,cex=1.15)
 }
 if(png){grDevices::png(paste0(file,'.png'),width=7.5,height=3.3,units='in',res=600,type='cairo');tryCatch(draw(),finally=dev.off())}
 if(pdf){cairo_pdf(paste0(file,'.pdf'),width=7.5,height=3.3,family='Arial');tryCatch(draw(),finally=dev.off())}
 if(interactive()||dev.cur()!=1)draw()
 invisible(data)
}
prepare_thermo_local_audiences<-function(audiences=c(0,5,10,15,20,29),
 reference_file='10_convincing_results/symmetric_beta04_2.rds'){
 reference<-readRDS(reference_file)$p;rows<-list();actors<-list()
 for(j in seq_along(audiences)){
  message('Finite personal audience: ',audiences[j])
  p<-thermo_local_audience(reference,audiences[j]);rows[[j]]<-thermo_local_summary(p)
  actors[[j]]<-convincing_optimize(p,rbind(c(2,.5),c(.5,2),c(1,1)))
  stopifnot(abs(sum(p$group_budget)-p$N*p$B)<1e-7,
   max(abs(p$moments[(length(p$effort)+1)/2,]-reference$moments[(length(p$effort)+1)/2,]))<1e-8)
  if(audiences[j]==0)stopifnot(max(abs(p$actors[,3]))<1e-8)
  if(audiences[j]==reference$N-1){
   stopifnot(max(abs(p$actors[,4]-reference$actors[,4]))<1e-6)
  }
 }
 ans<-list(summary=do.call(rbind,rows),actor_examples=actors,
  reference_file=reference_file,beta=reference$beta,phi=reference$phi,B=reference$B,
  local_group=reference$N,efficiency=reference$efficiency)
 saveRDS(ans,'10_thermodynamic_local_results.rds')
 write.csv(ans$summary,'10_thermodynamic_local_summary.csv',row.names=FALSE)
 plot_thermo_local_audience(ans,png=TRUE,pdf=TRUE)
 ans
}
