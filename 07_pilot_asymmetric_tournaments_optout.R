# A single opt-out talent-space panel, using the current population model.
# Example: pilot_asymmetric_optout(.8,.6,1.5,1,eta=.5,mode='responsive')
pilot_asymmetric_optout<-function(beta1,beta2,phi1,phi2,mode=c('fixed','responsive'),
    eta=.5,B=50,alim=c(0,2.5),nalpha=251L){
 mode<-match.arg(mode)
 h<-new.env(parent=globalenv());sys.source('07_optout_plot_helpers.R',envir=h)
 h$B<-B
 ref<-optout_fixed_scalers(c(beta1,beta2),c(phi1,phi2),B)$scaler
 if(mode=='fixed'){
  scaler<-ref;st<-optout_state(scaler,c(beta1,beta2),c(phi1,phi2),eta,B)
 }else{
  eq<-optout_responsive(c(beta1,beta2),c(phi1,phi2),eta,B);scaler<-eq$scaler;st<-eq$state
 }
 alpha<-seq(alim[1],alim[2],length.out=nalpha)
 sm<-h$strategy_map_optout(alpha,B,beta1,beta2,phi1,phi2,scaler[1],scaler[2],eta)
 par(mar=c(4,4.8,5,1.3),mgp=c(2.4,.7,0),pty='s',xpd=FALSE)
 plot(NA,xlim=alim,ylim=alim,asp=1,xaxs='i',yaxs='i',xlab='',ylab='',axes=FALSE)
 h$draw_colour_raster(alpha,h$allocation_colours(sm$x0,sm$x1,sm$x2))
 h$draw_maxxer_masks(alpha,sm);abline(0,1,lty=2,lwd=.9)
 for(curve in h$talent_hdr_curves())lines(curve[1,],curve[2,],col='#EEE',lwd=.9)
 hc<-h$half_curve(alpha,B,beta1,beta2,phi1,phi2,scaler[1],scaler[2],eta,alim)
 x<-hc$x;x[!hc$global]<-NA_real_;lines(x,hc$y,col='#FF0',lwd=3)
 pp<-optout_equal_partial(alpha,c(beta1,beta2),c(phi1,phi2),scaler,eta,B,alim)
 for(k in 1:2)lines(pp[,k],alpha,col='#FF0',lty=2,lwd=1)
 box();h$draw_population_bar(st$E0,st$E1,st$E2,alim)
 axis(1);axis(2,las=1)
 mtext(expression(alpha[1]~'(talent in tournament 1)'),1,2.4)
 mtext(expression(alpha[2]~'(talent in tournament 2)'),2,2.4)
 f<-function(x)formatC(x,format='f',digits=2)
 for(t in 1:2){
  at<-alim[1]+c(.3,.7)[t]*diff(alim);col<-c('#E06','#08F')[t]
  mtext(bquote(beta[.(t)]==.(f(c(beta1,beta2)[t]))),3,2.7,at=at,col=col)
  mtext(bquote(varphi[.(t)]==.(f(c(phi1,phi2)[t]))),3,1.3,at=at,col=col)
 }
 mtext(bquote(eta==.(f(eta))),3,.05)
 mtext(if(mode=='fixed')expression('Fixed '*E[1]==E[2])else expression('Responsive '*E[t]),2,4)
 invisible(list(allocation=sm,population=st,scaler=scaler))
}
