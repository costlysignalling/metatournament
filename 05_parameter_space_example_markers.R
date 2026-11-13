# Shared source-figure registry for the six 9+9-panel parameter figures.
# Only examples whose fixed parameters and axis coordinates are visible are drawn.
source('05_asymmetric_examples.R', local=TRUE)
source('05_asymmetric_supplement_examples.R', local=TRUE)
parameter_space_marker_labels <- c(fig3='Fig. 3',linear='Fig. 5 & S7',
                                   fig6='Fig. 6',s14='Fig. S14',s15='Fig. S15')
parameter_space_marker_colours <- c(fig3='#FFD700',linear='#C36BFF',
                                    fig6='#27D3D8',s14='#FF8C00',s15='#39F24A')
parameter_space_example_registry <- function() {
  # Read only the literal plotting parameter assignments, without running plots.
  env<-new.env(parent=baseenv())
  for(expr in as.list(parse('04_known_talents_strategy_map.R'))) {
    if(is.call(expr) && identical(expr[[1]],as.name('<-')) &&
       is.symbol(expr[[2]]) && as.character(expr[[2]]) %in% c('betas','phis'))
      eval(expr,envir=env)
  }
  stopifnot(is.numeric(env$betas),is.numeric(env$phis))
  fig3<-do.call(c,lapply(env$phis,function(p)lapply(env$betas,function(b)
    list(beta1=b,beta2=b,phi1=p,phi2=p))))
  list(fig3=fig3,
       linear=list(list(beta1=.8,beta2=.4,phi1=1,phi2=1)),
       fig6=asymmetric_examples,
       s14=asymmetric_supplement_examples$S14,
       s15=asymmetric_supplement_examples$S15)
}
parameter_space_example_sets <- function(plane=c('beta','phi'),panel_values,
                                        axis_limits) {
  plane <- match.arg(plane)
  fields <- c('beta1','beta2','phi1','phi2')
  lapply(parameter_space_example_registry(),function(scenarios) {
    values <- do.call(rbind,lapply(scenarios,function(sc)unlist(sc[fields],use.names=TRUE)))
    if(plane=='phi')values<-values[,c('phi1','phi2','beta1','beta2'),drop=FALSE]
    present<-function(x)any(abs(x-panel_values)<1e-9)
    keep<-vapply(values[,3],present,logical(1)) & vapply(values[,4],present,logical(1)) &
          values[,1]>=axis_limits[1] & values[,1]<=axis_limits[2] &
          values[,2]>=axis_limits[1] & values[,2]<=axis_limits[2]
    unique(values[keep,,drop=FALSE])
  })
}
parameter_space_example_vectors <- function(plane,panel_values,axis_limits) {
  lapply(parameter_space_example_sets(plane,panel_values,axis_limits),function(z)
    lapply(seq_len(nrow(z)),function(i)z[i,c('beta1','beta2','phi1','phi2')]))
}
parameter_space_visible_sets <- function(examples) {
  names(examples)[vapply(examples,function(z)
    if(is.matrix(z))nrow(z)>0L else length(z)>0L,logical(1))]
}
draw_parameter_example_legend <- function(examples,y,centre=NULL,cex=1.5,
                                           line_label=NULL) {
  shown<-parameter_space_visible_sets(examples)
  if(!length(shown))return(invisible(NULL))
  width<-par('pin')[1];if(is.null(centre))centre<-width/2
  labels<-unname(parameter_space_marker_labels[shown])
  sw<-.15;gap<-.07;item_gap<-.18;prefix_gap<-.14;line_width<-.25
  measure<-function(size) {
    total<-strwidth('Examples',cex=size,units='inches')+prefix_gap+
      sum(strwidth(labels,cex=size,units='inches'))+length(shown)*(sw+gap)+
      max(0,length(shown)-1)*item_gap
    if(!is.null(line_label))total<-total+item_gap+line_width+gap+
      strwidth(line_label,cex=size,units='inches')
    total
  }
  total<-measure(cex)
  if(total>.96*width){cex<-cex*.96*width/total;total<-measure(cex)}
  x<-centre-total/2
  text(x,y,'Examples',adj=c(0,.5),cex=cex)
  x<-x+strwidth('Examples',cex=cex,units='inches')+prefix_gap
  for(i in seq_along(shown)) {
    points(x+sw/2,y,pch=4,cex=1.3,lwd=3.5,col='#333333')
    points(x+sw/2,y,pch=4,cex=1.3,lwd=2.5,col=parameter_space_marker_colours[shown[i]])
    text(x+sw+gap,y,labels[i],adj=c(0,.5),cex=cex)
    x<-x+sw+gap+strwidth(labels[i],cex=cex,units='inches')+item_gap
  }
  if(!is.null(line_label)) {
    segments(x,y,x+line_width,y,col='white',lwd=1.1)
    segments(x,y,x+line_width,y,col='#333333',lwd=1.1,lty=2)
    text(x+line_width+gap,y,line_label,adj=c(0,.5),cex=cex)
  }
}
