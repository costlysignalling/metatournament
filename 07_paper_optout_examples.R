# Four selected opt-out columns for the main Paper. Run from this folder.
# Numerical population solutions and colours match the Supplement figures.
source("07_optout_plot_helpers.R")
scenarios<-list(
  list(beta1=.4,beta2=.4,phi1=1,phi2=1,eta=.5),
  list(beta1=.4,beta2=.4,phi1=1.5,phi2=1.5,eta=.5),
  list(beta1=.8,beta2=.6,phi1=1.5,phi2=1,eta=.5),
  list(beta1=.8,beta2=.6,phi1=1.5,phi2=1,eta=1))
prep<-lapply(scenarios,function(sc){
  cache<-population_cache(sc$beta1,sc$beta2,sc$phi1,sc$phi2,B)
  reference<-benchmark_scaler(cache,B)
  first<-population_state(reference[1],reference[2],cache,sc$eta,B)
  responsive<-responsive_state(cache,sc$eta,B,reference)
  list(pars=sc,cache=cache,k0=reference,fixed=first,eq=responsive)
})
draw_paper_optout_examples<-function(){
  alpha<-seq(alim[1],alim[2],length.out=nalpha)
  hdr_curves<-talent_hdr_curves(median,sd,probs=c(0.50,0.90,0.99))

  # EXACT Figure-07 layout.
  legw<-0.52
  widths<-c(1,1,1,1,legw)

  layout(
    matrix(c(1,2,3,4,9,
             5,6,7,8,9),nrow=2,byrow=TRUE),
    widths=widths,
    heights=c(1,1)
  )

  eta_at<-(sum(widths[1:4])+widths[5]/2)/sum(widths)

  # EXACT Figure-07 outer margins.
  par(oma=c(2,3.5,4.5,0.2))

  for(i in 1:2){
    for(j in 1:4){
      obj<-prep[[j]]
      sc<-obj$pars
      eta<-sc$eta

      if(i==1){
        scaler1<-obj$k0[1]
        scaler2<-obj$k0[2]
        pa<-obj$fixed
      } else {
        scaler1<-obj$eq$scaler[1]
        scaler2<-obj$eq$scaler[2]
        pa<-obj$eq$state
      }

      sm<-strategy_map_optout(
        alpha,B,
        sc$beta1,sc$beta2,sc$phi1,sc$phi2,
        scaler1,scaler2,eta,nx
      )

      cm<-allocation_colours(sm$x0,sm$x1,sm$x2)

      # EXACT Figure-07 panel graphics.
      par(
        mar=c(2.75,2.75,0.5,0.2),
        mgp=c(1.75,0.48,0),
        tcl=-0.23,
        pty="s",
        xpd=FALSE
      )

      plot(
        NA,
        xlim=alim,ylim=alim,
        xaxs="i",yaxs="i",
        axes=FALSE,xlab="",ylab="",asp=1
      )

      # Figure-07 diversifier/self-directed colour logic.
      draw_colour_raster(alpha,cm)

      # Pure maxxers in the same solid Figure-07 colours.
      draw_maxxer_masks(alpha,sm)

      abline(0,1,col="#000",lty=2,lwd=0.9)

      # Reference HDR contours: 50%, 90%, 99%.
      for(hc in hdr_curves){
        lines(hc[1,],hc[2,],col="#EEE",lwd=0.9)
      }

      # Figure-05 logic: exact 50:50 stationary curve, only where globally optimal.
      hc<-half_curve(
        alpha,B,
        sc$beta1,sc$beta2,sc$phi1,sc$phi2,
        scaler1,scaler2,eta,alim,nx
      )

      ii<-which(is.finite(hc$x) & hc$global)

      if(length(ii)>1){
        groups<-split(ii,cumsum(c(1,diff(ii)!=1)))

        for(r in groups){
          if(length(r)>1){
            lines(hc$x[r],hc$y[r],col=coldiag,lwd=3)
          }
        }
      }

      # Equal positive investments with some outside time. Global checks
      # exclude discontinuous switches between unequal allocations.
      equal_partial<-optout_equal_partial(alpha,c(sc$beta1,sc$beta2),
        c(sc$phi1,sc$phi2),c(scaler1,scaler2),eta,B,alim)
      for(branch in 1:2) lines(equal_partial[,branch],alpha,
        col=coldiag,lwd=1,lty=2)

      box()

      # EXACT Figure-07 population bar.
      draw_population_bar(pa$E0,pa$E1,pa$E2,alim)

      if(i==2){
        axis(1,at=seq(0,2.5,0.5),cex.axis=1.2)
        mtext(
          expression(alpha[1]~~"(talent in tournament 1)"),
          side=1,line=2.4,cex=0.95
        )
      }

      if(j==1){
        axis(2,at=seq(0,2.5,0.5),las=1,cex.axis=1.2)
        mtext(
          expression(alpha[2]~~"(talent in tournament 2)"),
          side=2,line=2.4,cex=0.95
        )
      }

      # EXACT Figure-05 two-line coloured parameter header.
      if(i==1){
        tx1<-alim[1]+0.48*diff(alim)
        tx2<-alim[1]+0.88*diff(alim)

        mtext(
          bquote(beta[1]==.(formatC(sc$beta1,format="f",digits=2))),
          side=3,line=1.55,at=tx1,cex=1.12,col=col1_max
        )
        mtext(
          bquote(beta[2]==.(formatC(sc$beta2,format="f",digits=2))),
          side=3,line=1.55,at=tx2,cex=1.12,col=col2_max
        )

        mtext(
          bquote(varphi[1]==.(formatC(sc$phi1,format="f",digits=2))),
          side=3,line=0.15,at=tx1,cex=1.12,col=col1_max
        )
        mtext(
          bquote(varphi[2]==.(formatC(sc$phi2,format="f",digits=2))),
          side=3,line=0.15,at=tx2,cex=1.12,col=col2_max
        )
      }

      if(i==1){
        mtext(bquote(eta==.(formatC(sc$eta,format="f",digits=2))),
          side=3,line=.15,at=alim[1],adj=0,cex=1.0,xpd=NA)
        mtext(LETTERS[j],side=3,line=1.8,at=alim[1],adj=0,cex=1.25,font=2,xpd=NA)
      }

      # Only show a warning in a responsive panel if the active solver did not
      # actually converge. This is diagnostic text, not a strategy label.
      if(i==2 && !obj$eq$status%in%c("active","inactive")){
        text(
          alim[1]+0.04*diff(alim),
          alim[2]-0.04*diff(alim),
          obj$eq$status,
          adj=c(0,1),
          cex=0.62,
          col="white"
        )
      }
    }
  }

  # Row labels follow Figures 05/07.
  mtext(
    expression("Fixed "*E[1]==E[2]),
    side=2,outer=TRUE,at=0.76,line=1.25,cex=1.25
  )

  mtext(
    expression("Responsive "*E[t]),
    side=2,outer=TRUE,at=0.26,line=1.25,cex=1.25
  )

  # EXACT Figure-07 legend panel margins and legend.
  par(mar=c(4.8,4.2,0.95,1.6),pty="m",xpd=FALSE)
  draw_allocation_legend()
}
if(isTRUE(getOption('metatournament.png',TRUE))){
  png('Figure_optout_selected_examples.png',width=3640/300,height=6.6,units='in',res=600)
  tryCatch(draw_paper_optout_examples(),finally=dev.off())
}
if(isTRUE(getOption('metatournament.pdf',TRUE))){
  cairo_pdf('Figure_optout_selected_examples.pdf',width=3640/300,height=6.6,family='Arial')
  tryCatch(draw_paper_optout_examples(),finally=dev.off())
}
if(interactive()) draw_paper_optout_examples()
