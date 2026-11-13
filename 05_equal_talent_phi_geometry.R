# Compile the analytical crossing classification and optional global-optimum check.
code <- c(readLines('05_energy_allocation_population_phispace.cpp',warn=FALSE),
          readLines('05_equal_talent_phi_geometry.cpp',warn=FALSE))
Rcpp::sourceCpp(code=paste(code,collapse='\n'),cacheDir='05_equal_talent_phi_cpp_cache')
