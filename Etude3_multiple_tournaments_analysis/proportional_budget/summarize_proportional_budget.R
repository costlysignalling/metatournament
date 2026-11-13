analysis_dir<-if(file.exists('participation_many.csv')) '.' else
 file.path('Etude3_multiple_tournaments_analysis','proportional_budget')
d<-read.csv(file.path(analysis_dir,'participation_many.csv'))
stopifnot(nrow(d)==3456,all(d$B==25*d$T))
totals<-aggregate(proportion~beta+phi+T+rep+N,d,sum)
stopifnot(max(abs(totals$proportion-1))<1e-10)
means<-aggregate(proportion~beta+phi+T+N+k+B,d,mean)
fine<-subset(means,N==max(N));coarse<-subset(means,N==min(N))
comparison<-merge(fine,coarse,by=c('beta','phi','T','k','B'),suffixes=c('_fine','_coarse'))
cat('Maximum nested-resolution difference, percentage points:',
    100*max(abs(comparison$proportion_fine-comparison$proportion_coarse)),'\n')
shifts<-aggregate(proportion~beta+phi+T+N+k,subset(d,N==max(N)),function(x)diff(range(x)))
cat('Largest range between shifted rules, percentage points:',100*max(shifts$proportion),'\n')
if(file.exists(file.path(analysis_dir,'participation_fixed_budget_baseline.csv'))) {
 fixed<-read.csv(file.path(analysis_dir,'participation_fixed_budget_baseline.csv'))
 old<-subset(fixed,T==2);restored<-subset(d,T==2)
 key<-c('beta','phi','T','rep','N','k')
 o<-do.call(order,old[key]);r<-do.call(order,restored[key])
 stopifnot(max(abs(old$proportion[o]-restored$proportion[r]))<1e-12)
 cat('All T=2, B=50 probabilities retained unchanged.\n')
}
analytic_file<-'../research_known_talents/analytic_participation.R'
if(!file.exists(analytic_file))analytic_file<-file.path('Etude3_multiple_tournaments_analysis','analytic_participation.R')
if(!file.exists(analytic_file))analytic_file<-'../analytic_participation.R'
source(analytic_file)
for(be in c(.4,.8)){
 exact<-participation_phi1(be,75,3)
 numerical<-subset(fine,beta==be&phi==1&T==3)
 numerical<-numerical[order(numerical$k),]
 error<-max(abs(exact-numerical$proportion))
 stopifnot(error<0.002)
 cat('Exact phi=1, T=3, B=75 check, beta=',be,': max error (percentage points)=',100*error,'\n')
}
write.csv(fine,file.path(analysis_dir,'participation_many_summary.csv'),row.names=FALSE)
write.csv(comparison,file.path(analysis_dir,'participation_many_convergence.csv'),row.names=FALSE)
print(subset(fine,k==T),row.names=FALSE)
