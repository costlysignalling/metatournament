Sys.setenv(PKG_CPPFLAGS='-DRCPP_INTERFACE')
analysis_dir<-if(file.exists('known_talents_participation.cpp')) '.' else
  file.path('Etude3_multiple_tournaments_analysis','T10_phi125')
Rcpp::sourceCpp(file.path(analysis_dir,'known_talents_participation.cpp'),
                cacheDir=file.path(analysis_dir,'participation_cpp_cache'))
grid_count<-function(a,beta,phi,n=1000L){
 xx<-seq(0,50,length.out=n+1L)
 value<-rep(-Inf,n+1L);value[1]<-0;support<-integer(n+1L)
 for(aa in a){
  ww<-expm1(phi*log1p(aa*expm1(beta*log1p(xx))))
  nextvalue<-numeric(n+1L);nextsupport<-integer(n+1L)
  for(j in 0:n){
   old<-(j:0)+1L;vals<-value[old]+ww[1:(j+1L)];pick<-which.max(vals)
   nextvalue[j+1L]<-vals[pick]
   nextsupport[j+1L]<-support[old[pick]]+(pick>1L)
  }
  value<-nextvalue;support<-nextsupport
 }
 support[n+1L]
}
cases<-list(
 list(.4,1.25,c(7,6.1,.8,.7,.6,.5,.4,.3,.2,.1)),
 list(.4,1.25,rep(6.1,10)),list(.4,1.25,rep(.8,10)),
 list(.8,1.25,rep(.8,10)),list(.8,1.25,c(1.01,rep(.8,9))),
 list(.8,1.25,c(.8,.81)),list(.8,1.25,seq(.1,.9,length.out=10)),
 list(.4,1.5,rep(6,10)),list(.4,1.5,c(4,3.1,3.05,.9,.8,.7,.6,.5,.4,.3)),
 list(.8,1.5,rep(.01,10)),list(.8,1.5,seq(.02,.11,length.out=10)),
 list(.4,.66,seq(.1,2,length.out=10)),list(.8,1,seq(.1,2,length.out=10)))
for(i in seq_along(cases)){
 beta<-cases[[i]][[1]];phi<-cases[[i]][[2]];a<-cases[[i]][[3]]
 actual<-participation_count(matrix(a,nrow=1),beta,phi)[1]
 reference<-grid_count(a,beta,phi)
 if(actual!=reference)reference<-grid_count(a,beta,phi,2000L)
 stopifnot(actual==reference)
 cat('Independent grid check passed:',i,'beta=',beta,'phi=',phi,
     'T=',length(a),'k=',actual,'\n')
}
cat('All added-parameter and T=10 optimization checks passed.\n')
