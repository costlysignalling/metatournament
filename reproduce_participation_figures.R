# Run source("reproduce_participation_figures.R") from the scripts folder.
# Then reproduce_participation_figures() redraws both figures from saved data.
# reproduce_participation_figures(recompute=TRUE) recalculates all probabilities.
# See REPRODUCING_PARTICIPATION_FIGURES.md for requirements and file inventory.
reproduce_participation_figures <- local({
  script_dir <- dirname(normalizePath(sys.frame(1)$ofile, mustWork=TRUE))
  function(recompute=FALSE, figures=c("both","fixed","proportional")) {
    figures <- match.arg(figures)
    old_dir <- getwd()
    on.exit(setwd(old_dir), add=TRUE)
    old_flags <- Sys.getenv("PKG_CPPFLAGS", unset=NA_character_)
    on.exit(if(is.na(old_flags)) Sys.unsetenv("PKG_CPPFLAGS") else
      Sys.setenv(PKG_CPPFLAGS=old_flags), add=TRUE)
    setwd(script_dir)
    selected <- if(figures=="both") c("fixed","proportional") else figures
    for(which in selected) {
      subdir <- if(which=="fixed") "T10_phi125" else "proportional_budget"
      analysis_dir <- file.path(script_dir,"Etude3_multiple_tournaments_analysis",subdir)
      summary_script <- if(which=="fixed") "summarize_T10.R" else "summarize_proportional_budget.R"
      required <- c("known_talents_participation.cpp",summary_script)
      paths <- c(file.path(analysis_dir,required),
        file.path(script_dir,"Etude3_multiple_tournaments_analysis","analytic_participation.R"))
      if(any(!file.exists(paths))) stop("Missing reproduction files: ",paste(paths[!file.exists(paths)],collapse=", "))
      if(recompute) {
        if(!requireNamespace("Rcpp",quietly=TRUE)) stop("Install Rcpp first: install.packages('Rcpp')")
        message("Recalculating every probability for ",which," budget; this can take hours.")
        setwd(analysis_dir)
        Sys.setenv(PKG_CPPFLAGS=paste(if(is.na(old_flags)) "" else old_flags,"-DRCPP_INTERFACE"))
        engine <- new.env(parent=globalenv())
        Rcpp::sourceCpp("known_talents_participation.cpp",env=engine,
          cacheDir="participation_cpp_cache")
        # Full calculation: no previously calculated probabilities are reused.
        if(which=="fixed") engine$participation_compute(262144L,4L,"participation_many.csv",incremental=FALSE) else {
          engine$participation_compute(262144L,4L,"participation_many.csv",proportional=FALSE)
          d <- read.csv("participation_many.csv")
          d$B <- 25*d$T
          write.csv(d,"participation_many.csv",row.names=FALSE)
        }
        sys.source(summary_script,envir=new.env(parent=globalenv()))
      }
      setwd(script_dir)
      plot_script <- if(which=="fixed") "04_known_talents_multiple_tournaments_T10.R" else
        "04_known_talents_multiple_tournaments_proportional_budget.R"
      if(!file.exists(file.path(analysis_dir,"participation_many_summary.csv"))) stop("Missing summary; use recompute=TRUE.")
      # Run in the selected calculation directory so the plot reads its own CSV.
      # Output is copied back to the scripts folder, alongside the manuscript.
      setwd(analysis_dir)
      sys.source(file.path(script_dir,plot_script),envir=new.env(parent=globalenv()))
      png_name <- if(which=="fixed") "Figure_known_talents_multiple_tournaments_T10.png" else
        "Figure_known_talents_multiple_tournaments_proportional_budget.png"
      file.copy(png_name,file.path(script_dir,png_name),overwrite=TRUE)
      message("Updated ",png_name)
    }
    invisible(NULL)
  }
})
