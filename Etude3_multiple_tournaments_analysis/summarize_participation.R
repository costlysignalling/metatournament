analysis_dir<-if(file.exists('participation_many.csv')) '.' else 'Etude3_multiple_tournaments_analysis'
d<-read.csv(file.path(analysis_dir,'participation_many.csv'))
stopifnot(nrow(d)==6*sum(2:7)*4*2)
totals<-aggregate(proportion~beta+phi+T+rep+N,d,sum)
stopifnot(max(abs(totals$proportion-1))<1e-10)
means<-aggregate(proportion~beta+phi+T+N+k,d,mean)
fine<-subset(means,N==max(N))
coarse<-subset(means,N==min(N))
comparison<-merge(fine,coarse,by=c('beta','phi','T','k'),suffixes=c('_fine','_coarse'))
cat('Maximum nested-resolution difference, percentage points:',
    100*max(abs(comparison$proportion_fine-comparison$proportion_coarse)),'\n')
shift_range<-aggregate(proportion~beta+phi+T+N+k,subset(d,N==max(N)),function(x) diff(range(x)))
cat('Largest range between shifted rules, percentage points:',100*max(shift_range$proportion),'\n')
analytic_file<-file.path(analysis_dir,'analytic_participation.R')
if(!file.exists(analytic_file))analytic_file<-'../research_known_talents/analytic_participation.R'
source(analytic_file)
for(be in c(.4,.8))for(tt in c(2,3)){
  exact<-participation_phi1(be,50,tt)
  numerical<-subset(fine,beta==be&phi==1&T==tt)
  numerical<-numerical[order(numerical$k),]
  error<-max(abs(exact-numerical$proportion))
  stopifnot(error<0.002)
  cat('Exact phi=1 check, beta=',be,' T=',tt,': max error (percentage points)=',100*error,'\n')
}
write.csv(fine,file.path(analysis_dir,'participation_many_summary.csv'),row.names=FALSE)
write.csv(comparison,file.path(analysis_dir,'participation_many_convergence.csv'),row.names=FALSE)
print(subset(fine,k==T|k==1),row.names=FALSE)
