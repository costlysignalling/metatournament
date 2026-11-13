# Shared numbering for separately rendered Paper and Supplement.
# Numbering is rebuilt from Supplement.Rmd on every render, never hard-coded.
supplement_outline <- function(path='Supplement.Rmd') {
  lines <- readLines(path,warn=FALSE,encoding='UTF-8')
  inside_chunk <- FALSE
  counts <- integer(6)
  out <- list()
  for(line in lines) {
    if(grepl('^```',line)) { inside_chunk <- !inside_chunk; next }
    if(inside_chunk || !grepl('^#{1,6} ',line)) next
    level <- nchar(sub(' .*','',line))
    counts[level] <- counts[level]+1L
    if(level<6L) counts[(level+1L):6L] <- 0L
    number <- paste(counts[seq_len(level)],collapse='.')
    id <- if(grepl('\\{#[^}]+\\}',line))
      sub('.*\\{#([^} ]+).*','\\1',line) else ''
    out[[length(out)+1L]] <- list(line=line,id=id,number=number)
  }
  out
}
.supplement_sections <- supplement_outline()
# Preserve the existing bookdown configuration, but select the figure prefix
# per manuscript. bookdown 0.48 stores this in its internal options object.
.bookdown_config <- bookdown:::load_config()
if(is.null(.bookdown_config)) .bookdown_config <- list()
.bookdown_config$language$label$fig <-
  if(grepl('^Supplement\\.',basename(knitr::current_input()))) 'Figure S' else 'Figure '
bookdown:::opts$set(config=.bookdown_config)
supplement_ref <- function(id) {
  matches <- Filter(function(x) identical(x$id,id),.supplement_sections)
  if(length(matches)!=1L) stop('Missing or duplicate Supplement ID: ',id)
  ext <- if(knitr::is_html_output()) 'html' else
    if(identical(knitr::pandoc_to(),'docx')) 'docx' else 'Rmd'
  sprintf('[Supplement %s](Supplement.%s#%s)',matches[[1L]]$number,ext,id)
}

# Figure citations in the Supplement use the current Paper figure order.
paper_figure_ref <- function(id,panel='') {
  lines <- readLines('Paper.Rmd',warn=FALSE,encoding='UTF-8')
  figures <- lines[grepl('^```\\{r .*fig\\.cap=',lines)]
  labels <- sub('^```\\{r ([^, ]+).*','\\1',figures)
  number <- match(id,labels)
  if(is.na(number)) stop('Missing Paper figure ID: ',id)
  ext <- if(knitr::is_html_output()) 'html' else
    if(identical(knitr::pandoc_to(),'docx')) 'docx' else 'Rmd'
  sprintf('[Figure%s%s of the main text](Paper.%s#%s-figure)',
          number,panel,ext,id)
}

# Supplement figure citations in Paper follow the current Supplement chunk order.
supplement_figure_ref <- function(id,panel='') {
  lines <- readLines('Supplement.Rmd',warn=FALSE,encoding='UTF-8')
  figures <- lines[grepl('^```\\{r .*fig\\.cap=',lines)]
  labels <- sub('^```\\{r ([^, ]+).*','\\1',figures)
  number <- match(id,labels)
  if(is.na(number)) stop('Missing Supplement figure ID: ',id)
  ext <- if(knitr::is_html_output()) 'html' else
    if(identical(knitr::pandoc_to(),'docx')) 'docx' else 'Rmd'
  sprintf('[Figure S%s%s](Supplement.%s#%s-figure)',number,panel,ext,id)
}

manuscript_document_hook <- function(text) {
  # knitr can supply a vector of document fragments, not one scalar string.
  text <- paste(text,collapse='\n')
  ext <- if(knitr::is_html_output()) 'html' else
    if(identical(knitr::pandoc_to(),'docx')) 'docx' else 'Rmd'
  text <- gsub('(Paper|Supplement)\\.Rmd#',paste0('\\1.',ext,'#'),text)
  if(grepl('^Supplement\\.',basename(knitr::current_input()))) {
    # Literal title prefixes keep numbering identical in Word and HTML.
    # Section numbering is disabled in the output formats to avoid duplication.
    lines <- strsplit(text,'\n',fixed=TRUE)[[1L]]
    for(section in .supplement_sections) {
      at <- which(lines==section$line)
      lines[at] <- sub('^(#+) ',paste0('\\1 ',section$number,' '),section$line)
    }
    text <- paste(lines,collapse='\n')
  }
  text
}
knitr::knit_hooks$set(document=manuscript_document_hook)
