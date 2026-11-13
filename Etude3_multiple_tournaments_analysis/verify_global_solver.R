Sys.setenv(PKG_CPPFLAGS='-DRCPP_INTERFACE')
analysis_dir<-if(file.exists('known_talents_participation.cpp')) '.' else 'Etude3_multiple_tournaments_analysis'
Rcpp::sourceCpp(file.path(analysis_dir,'known_talents_participation.cpp'),
                cacheDir=file.path(analysis_dir,'participation_cpp_cache'))

grid_optimum<-function(a,beta,phi,n=1000L){
  B<-50;xx<-seq(0,B,length.out=n+1L)
  value<-rep(-Inf,n+1L);value[1]<-0
  support<-integer(n+1L)
  for(aa in a){
    ww<-expm1(phi*log1p(aa*expm1(beta*log1p(xx))))
    nextval<-numeric(n+1L);nextsupport<-integer(n+1L)
    for(j in 0:n){
      old<-(j:0)+1L
      candidates<-value[old]+ww[1:(j+1L)]
      pick<-which.max(candidates)
      nextval[j+1L]<-candidates[pick]
      nextsupport[j+1L]<-support[old[pick]]+(pick>1L)
    }
    value<-nextval;support<-nextsupport
  }
  support[n+1L]
}
cases<-list(
 c(10,3.2,.5),c(4,3.1,3.05,.9,.8,.6,.3),c(100,50,10,5,4,3.1,2),
 rep(6,7),rep(3.1,7),c(3.01,3.01,.4),c(.3,.5,1),
 c(.01,.01),rep(.01,7),c(.1,.05,.01),c(.49,.48,.2),rep(.02,7),
 c(.02,.03,.04,.05,.06,.07,.08),c(.05,.06),rep(.05,7))
for(i in seq_along(cases)){
 beta<-if(i<=7) .4 else .8
 a<-cases[[i]]
 exact<-participation_count(matrix(a,nrow=1),beta,1.5)[1]
 coarse<-grid_optimum(a,beta,1.5,500)
 fine<-grid_optimum(a,beta,1.5,1000)
 if(exact!=fine){
   finer<-grid_optimum(a,beta,1.5,2000)
   stopifnot(exact==finer)
 }
 cat('Global solver vs independent simplex dynamic program:',
     'case=',i,' T=',length(a),' k=',exact,' grid=',coarse,fine,'\n')
}
cat('All global optimization checks passed.\n')
