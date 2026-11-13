Sys.setenv(PKG_CPPFLAGS='-DRCPP_INTERFACE')
analysis_dir<-if(file.exists('known_talents_participation.cpp')) '.' else
  file.path('Etude3_multiple_tournaments_analysis','proportional_budget')
Rcpp::sourceCpp(file.path(analysis_dir,'known_talents_participation.cpp'),
                cacheDir=file.path(analysis_dir,'participation_cpp_cache'))
fixed<-read.csv(file.path(analysis_dir,'participation_fixed_budget_baseline.csv'))
baseline<-subset(fixed,T==2)
stopifnot(max(baseline$N)==262144,nrow(baseline)==128)
new_file<-file.path(analysis_dir,'participation_proportional_increment.csv')
participation_compute(262144L,4L,new_file,proportional=TRUE)
new<-read.csv(new_file)
stopifnot(nrow(new)==3328)
all<-rbind(baseline,new)
stopifnot(nrow(all)==3456,!anyDuplicated(all[c('beta','phi','T','rep','N','k')]))
all$B<-25*all$T
write.csv(all,file.path(analysis_dir,'participation_many.csv'),row.names=FALSE)
cat('All proportional-budget panels completed and verified.\n')
