# source('10_explore_convincing.R')
# explore_convincing(beta1=.4,beta2=.4,phi1=1,phi2=1,N=30)
# Draws in R; exports only with png=TRUE and/or pdf=TRUE.
source('10_convincing_figures.R')
explore_convincing<-function(beta1=.4,beta2=.4,phi1=1,phi2=1,
 B=50,N=30,efficiency=.1,actor=c(2,.5),aligned=FALSE,
 n=8192L,npressure=81L,nalpha=361L,png=FALSE,pdf=FALSE,
 file='Figure_convincing_explorer'){
 p<-convincing_prepare(c(beta1,beta2),c(phi1,phi2),B,N,efficiency,
  n=n,npressure=npressure,aligned=aligned)
 map<-convincing_map(p,seq(0,2.5,length.out=nalpha))
 print(convincing_optimize(p,matrix(actor,nrow=1)));print(p$summary)
 draw_convincing_set(data=list(list(p=p,map=map)),file=file,png=png,pdf=pdf)
 invisible(list(p=p,map=map))
}
