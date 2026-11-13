if(!requireNamespace('Rcpp',quietly=TRUE)) stop('Install Rcpp before running this script.')
Sys.setenv(PKG_CPPFLAGS='-DRCPP_INTERFACE')
analysis_dir<-if(file.exists('known_talents_participation.cpp')) '.' else 'Etude3_multiple_tournaments_analysis'
Rcpp::sourceCpp(file.path(analysis_dir,'known_talents_participation.cpp'),
                cacheDir=file.path(analysis_dir,'participation_cpp_cache'))
args<-commandArgs(trailingOnly=TRUE)
N<-if(length(args)) as.integer(args[1]) else 65536L
reps<-if(length(args)>1) as.integer(args[2]) else 4L
output<-if(length(args)>2) args[3] else file.path(analysis_dir,'participation_many.csv')
participation_compute(N,reps,output)
stopifnot(file.exists(output),nrow(read.csv(output))==6*sum(2:7)*reps*2)
cat('Completed and verified all participation rows.\n')
