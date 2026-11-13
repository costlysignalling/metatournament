analysis_dir<-if(file.exists('participation_many.csv')) '.' else
 file.path('Etude3_multiple_tournaments_analysis','T10_phi125')
d<-read.csv(file.path(analysis_dir,'participation_many.csv'))
stopifnot(nrow(d)==8*sum(2:10)*4*2)
totals<-aggregate(proportion~beta+phi+T+rep+N,d,sum)
stopifnot(max(abs(totals$proportion-1))<1e-10)
means<-aggregate(proportion~beta+phi+T+N+k,d,mean)
fine<-subset(means,N==max(N));coarse<-subset(means,N==min(N))
comparison<-merge(fine,coarse,by=c('beta','phi','T','k'),suffixes=c('_fine','_coarse'))
cat('Maximum nested-resolution difference, percentage points:',
    100*max(abs(comparison$proportion_fine-comparison$proportion_coarse)),'\n')
shifts<-aggregate(proportion~beta+phi+T+N+k,subset(d,N==max(N)),function(x)diff(range(x)))
cat('Largest range between shifted rules, percentage points:',100*max(shifts$proportion),'\n')
if(file.exists(file.path(analysis_dir,'participation_T7_baseline.csv'))) {
 old<-read.csv(file.path(analysis_dir,'participation_T7_baseline.csv'))
 restored<-subset(d,T<=7&phi!=1.25)
 key<-c('beta','phi','T','rep','N','k')
 o<-do.call(order,old[key]);r<-do.call(order,restored[key])
 stopifnot(max(abs(old$proportion[o]-restored$proportion[r]))<1e-12)
 cat('All previous T<=7 probabilities retained unchanged.\n')
}
write.csv(fine,file.path(analysis_dir,'participation_many_summary.csv'),row.names=FALSE)
write.csv(comparison,file.path(analysis_dir,'participation_many_convergence.csv'),row.names=FALSE)
print(subset(fine,phi==1.25&k==T),row.names=FALSE)
