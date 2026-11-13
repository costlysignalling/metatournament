# Self-contained Shiny application. Run shiny::runApp('06B_cube_web').
# The bundled engine is copied from the manuscript's R/C++ source files.
library(shiny)
library(jsonlite)
model <- new.env(parent=globalenv())
for(e in parse('engine/06B_predict_many_tournaments.R'))
  if(is.call(e) && identical(e[[1]],as.name('<-')))eval(e,model)
sys.source('engine/06B_three_tournament_geometry.R',envir=model)
default <- readRDS('default_cube.rds')

ui <- fluidPage(
  tags$head(tags$link(rel='stylesheet',href='cube.css'),
            tags$script(src='cube.js')),
  titlePanel('Three asymmetric tournaments'),
  p('Change the training exponents β and spoils exponents φ, then recalculate. ',
    'Drag either cube to rotate both views. Budget B = 50; independent ',
    'lognormal talents have median 1 and SD 0.5.'),
  div(class='parameters',lapply(1:3,function(t)div(class='tournament',
    h3(paste('Tournament',t),style=paste0('color:',c('#C90056','#0069C7','#28794A')[t])),
    sliderInput(paste0('b',t),HTML(paste0('β<sub>',t,'</sub>')),min=.05,max=.95,
                value=c(.9,.6,.6)[t],step=.01),
    sliderInput(paste0('p',t),HTML(paste0('φ<sub>',t,'</sub>')),min=.1,max=3.1,
                value=c(.8,1.5,1)[t],step=.05)))),
  div(class='actions',actionButton('compute','Recalculate cubes',class='btn-primary'),
      actionButton('restore','Restore manuscript example'),
      selectInput('quality','Calculation',choices=c('Preview'='quick','Finer grid'='fine'),width='180px'),
      downloadButton('snapshot','Download HTML snapshot')),
  div(id='cube-status',role='status'),
  div(class='view-controls',
    tags$label(tags$input(id='specialists',type='checkbox',checked='checked'),' Specialist boundaries'),
    tags$label(tags$input(id='equal',type='checkbox',checked='checked'),' Equal leading investments'),
    tags$label('Surface opacity ',tags$input(id='opacity',type='range',min='.05',max='.65',step='.01',value='.28')),
    tags$button(id='reset',type='button','Reset view')),
  div(class='key',lapply(1:4,function(t)tags$span(
    tags$i(class='sw',style=paste0('background:',c('#EE0066','#0088FF','#40A56A','#D9AE00')[t])),
    c('Tournament 1','Tournament 2','Tournament 3','Equal leading investments')[t]))),
  div(class='cubes',
    tags$section(h2('Fixed equal Eₜ reference'),tags$canvas(id='fixed-cube',`aria-label`='First-response allocation cube')),
    tags$section(h2('Responsive Eₜ'),tags$canvas(id='responsive-cube',`aria-label`='Responsive allocation cube'))),
  p(id='applied-parameters'),
  p('Gold surfaces mark equal positive investments in two jointly leading tournaments; ',
    'each need not receive half the budget. The visible cube spans talents 0–2.5. ',
    'Population calculations include the full talent distribution outside this cube.'),
  tags$details(tags$summary('Calculation details'),verbatimTextOutput('diagnostics'))
)

server <- function(input,output,session) {
  current <- reactiveVal(default)
  send <- function(x)session$sendCustomMessage('cube-data',list(data=x$data,
    parameters=sprintf('Displayed: β = (%s), φ = (%s), B = 50. Population nodes: %s; cube grid: %s³.',
      paste(formatC(x$beta,digits=2,format='f'),collapse=', '),
      paste(formatC(x$phi,digits=2,format='f'),collapse=', '),x$n,x$n_cube),
    beta=x$beta,phi=x$phi))
  session$onFlushed(function()send(isolate(current())),once=TRUE)
  observeEvent(input$restore,{
    for(t in 1:3){updateSliderInput(session,paste0('b',t),value=default$beta[t]);
                 updateSliderInput(session,paste0('p',t),value=default$phi[t])}
    current(default);send(default)
  })
  observeEvent(input$compute_request,{
    beta <- as.numeric(input$compute_request$beta)
    phi <- as.numeric(input$compute_request$phi)
    stopifnot(length(beta)==3,length(phi)==3,all(is.finite(c(beta,phi))),all(beta>0 & beta<1),all(phi>0))
    n <- if(input$compute_request$quality=='fine')16384L else 4096L
    n_cube <- if(input$compute_request$quality=='fine')35L else 21L
    session$sendCustomMessage('cube-status',list(message='Calculating new population and boundaries…'))
    tryCatch(withProgress(message='Calculating allocation cubes',value=0,{
      fixed <- model$many_random_reference(beta,phi,50)
      A <- model$many_lognormal_population(n,3)
      setProgress(.08,detail='Solving the responsive population')
      responsive <- model$many_responsive(A,beta,phi,50,start=fixed$scaler)
      if(!responsive$converged)stop(sprintf('Responsive population did not converge (residual %.3g). No replacement cube is shown.',responsive$residual))
      setProgress(.3,detail='Tracing the allocation surfaces')
      count <- 0
      geometry <- model$many_cube_meshes(beta,phi,50,list(fixed$scaler,responsive$scaler),n_cube=n_cube,
        progress=function(detail){count <<- count+1;setProgress(min(.95,.3+.65*count/8),detail=detail)})
      x <- list(beta=beta,phi=phi,n=n,n_cube=n_cube,
        residual=responsive$residual,max_gap=responsive$max_allocation_gap,
        E=responsive$E,data=model$many_mesh_payload(geometry$meshes))
      current(x);send(x);setProgress(1)
    }),error=function(e){
      session$sendCustomMessage('cube-status',list(message=paste('Calculation stopped:',conditionMessage(e))))
      showNotification(conditionMessage(e),type='error',duration=NULL)
    })
  })
  output$diagnostics <- renderText({x<-current();paste0(
    'Population consistency residual: ',format(x$residual,scientific=TRUE),
    '\nLargest relative actor-utility gap bound: ',format(x$max_gap,scientific=TRUE),
    '\nResponsive mean investments: ',paste(round(x$E,3),collapse=', '),
    '\nInitial allocations: uniform over the budget simplex, independently of talents.',
    '\nPreview and finer modes use 4,096 and 16,384 deterministic population vectors.',
    '\nThe supplied manuscript example uses 65,536 vectors.',
    '\nSurface resolution and population integration are separate numerical approximations.')})
  output$snapshot <- downloadHandler(filename=function()'Three_tournament_cube.html',
    content=function(file){
      x<-current();template<-paste(readLines('www/snapshot.html',warn=FALSE),collapse='\n')
      template<-sub('__CUBE_DATA__',toJSON(x$data,auto_unbox=TRUE,digits=7),template,fixed=TRUE)
      description<-sprintf('β = (%s), φ = (%s), B = 50; %s population nodes.',
        paste(x$beta,collapse=', '),paste(x$phi,collapse=', '),x$n)
      template<-sub('__CUBE_DESCRIPTION__',description,template,fixed=TRUE)
      writeLines(template,file,useBytes=TRUE)
    },contentType='text/html')
}
shinyApp(ui,server)

