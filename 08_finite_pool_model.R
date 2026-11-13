source('08_boredom_model.R') # shared deterministic talent and quadrature helpers
finite_pool_compile<-local({core<-NULL;function(){
 if(is.null(core)){core<<-new.env();Rcpp::sourceCpp('08_finite_pool_core.cpp',env=core)}
 core
}})
finite_pool_population<-function(N=30,beta=c(.8,.4),phi=c(1,1.5),B=50,n=32768L){
 den<-boredom_qmc_density(n);ref<-boredom_reference(beta,phi,B,compulsory=0)
 core<-finite_pool_compile()
 first<-core$finite_pool_density_cpp(den$talent,den$weights,ref,beta[1],beta[2],phi[1],phi[2],B,.5,N,160L,1L,2e-6,FALSE)
 eq<-core$finite_pool_density_cpp(den$talent,den$weights,first$mapped[1:3],beta[1],beta[2],phi[1],phi[2],B,.5,N,160L,120L,2e-6,TRUE)
 if(!eq$converged)stop('Finite-pool population unresolved for N=',N,': ',eq$residual)
 list(beta=beta,phi=phi,B=B,N=N,reference=ref,first=first,responsive=eq,order=n)
}
finite_pool_alloc<-function(p,mode,actor){
 m<-if(mode=='first')p$reference else p$responsive$moment
 bg<-(p$N-1)*c(m[1],p$B-m[1],m[2:3])
 finite_pool_compile()$finite_pool_map_cpp(as.matrix(actor),bg,p$beta[1],p$beta[2],p$phi[1],p$phi[2],p$B,.5,240L)
}
prepare_finite_pool_examples<-function(){
 dir.create('09_harvesting_results',showWarnings=FALSE)
 p<-lapply(c(100,30,10,5),function(N){message('Finite N=',N);finite_pool_population(N)})
 saveRDS(p,'09_harvesting_results/finite_examples.rds');invisible(p)
}
