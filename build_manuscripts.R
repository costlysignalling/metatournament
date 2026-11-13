# Rebuild both files together so cross-document numbers stay synchronized.
# Run with Rscript build_manuscripts.R from the scripts directory.
for(package in c('knitr','rmarkdown','bookdown')) {
  if(!requireNamespace(package,quietly=TRUE)) stop('Missing package: ',package)
}
for(doc in c('Paper.Rmd','Supplement.Rmd')) {
  for(format in c('bookdown::word_document2','bookdown::html_document2')) {
    rmarkdown::render(doc,output_format=format,envir=new.env(parent=globalenv()),quiet=FALSE)
  }
}
