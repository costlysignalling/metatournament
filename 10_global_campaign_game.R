# Additive global broadcasts in a finite population. Campaigns are chosen
# first; the training subgame then responds to all other people's messages.
# Received persuasion enters the objective in BOTH decisions. Own messages
# do not change one's own preference weights. Unweighted material utility
# is reported separately. No unresolved profile is labelled an equilibrium.
source('08_boredom_model.R')
broadcast_compile<-local({e<-NULL;function(){if(is.null(e)){
 e<<-new.env();Rcpp::sourceCpp('10_global_campaign_game.cpp',env=e,
  cacheDir=file.path(tempdir(),'10_global_campaign_game'))};e}})
broadcast_talents<-function(N=30L){
 stopifnot(N>=4,N%%2==0)
 radical<-function(i,b){u<-numeric(length(i));f<-1/b;while(any(i>0)){
  u<-u+f*(i%%b);i<-i%/%b;f<-f/b};u}
 i<-seq_len(N/2)+113L
 a<-exp(boredom_slog()*qnorm(cbind(radical(i,2),radical(i,3))))
 out<-matrix(0,N,2);out[seq(1,N,2),]<-a;out[seq(2,N,2),]<-a[,2:1]
 out
}
broadcast_solve<-function(beta=c(.4,.4),phi=c(1,1),B=50,efficiency=.1,
 N=30L,talent=broadcast_talents(N),start=rep(0,nrow(talent)),
 nx=24L,grid=16L,maxit=100L,tolerance=.005,sequential=TRUE,enforce_symmetry=FALSE){
 stopifnot(length(beta)==2,length(phi)==2,all(beta>0&beta<1),all(phi>0),
  B>0,efficiency>=0,ncol(talent)==2,all(talent>0),length(start)==nrow(talent),
  all(abs(start)<B))
 z<-broadcast_compile()$broadcast_campaign_cpp(talent,beta,phi,start,B,efficiency,
  nx,grid,maxit,tolerance,sequential,enforce_symmetry)
 z$x2<-B-abs(z$campaign)-z$x1
 stopifnot(max(abs(z$x1+z$x2+abs(z$campaign)-B))<1e-8)
 w1<-boredom_w(z$x1,talent[,1],beta[1],phi[1]);w2<-boredom_w(z$x2,talent[,2],beta[2],phi[2])
 E<-c(sum(z$x1+pmax(z$campaign,0)),sum(z$x2+pmax(-z$campaign,0)))
 M<-c(sum(w1),sum(w2));part1<-if(M[1]>0)E[1]*w1/M[1]else rep(0,nrow(talent))
 part2<-if(M[2]>0)E[2]*w2/M[2]else rep(0,nrow(talent))
 z$material<-part1+part2
 z$received_pressure<-efficiency*(sum(z$campaign)-z$campaign)
 z$perceived_scaled<-exp(pmin(z$received_pressure,0))*part1+exp(pmin(-z$received_pressure,0))*part2
 z$E<-E;z$M<-M
 z$budget_shares<-colMeans(cbind(z$x1,z$x2,pmax(z$campaign,0),pmax(-z$campaign,0)))/B
 if(!isTRUE(z$campaign_converged))warning('The campaign/training game is unresolved; this result is not an equilibrium.')
 list(beta=beta,phi=phi,B=B,efficiency=efficiency,N=nrow(talent),talent=talent,
  result=z,settings=list(nx=nx,grid=grid,tolerance=tolerance))
}
