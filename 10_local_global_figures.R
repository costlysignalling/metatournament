# Allocation and convincing examples used in the Supplement.
# source('10_local_global_prepare_figures.R'); prepare_local_global_figures()
# Replotting reads saved calculations. PNG/PDF are optional, FALSE by default.
source('10_convincing_model.R')
convincing_hdr<-function(){s<-boredom_slog();th<-seq(0,2*pi,length.out=1001)
 lapply(qchisq(c(.5,.9,.99),2,ncp=2*s*s),function(q)exp(-s*s+s*sqrt(q)*rbind(cos(th),sin(th))))}
convincing_polygon<-function(a,f,col){
 iso<-isoband::isobands(a,a,t(f),0,max(1,max(f)+1))
 if(length(iso[[1]]$x))plot(sf::st_sfc(isoband::iso_to_sfg(iso)[[1]]),add=TRUE,col=col,border=NA)
}
# Vector colour bands retain abrupt strategy switches without raster blur.
convincing_colour_bands<-function(a,z,cols,limits,step=.005){
 z[]<-pmax(limits[1],pmin(limits[2],z))
 breaks<-seq(limits[1],limits[2],length.out=ceiling(diff(limits)/step)+1L)
 lo<-head(breaks,-1);hi<-tail(breaks,-1)
 lo[1]<-lo[1]-1e-10;hi[length(hi)]<-hi[length(hi)]+1e-10
 band<-isoband::isobands(a,a,t(z),lo,hi)
 index<-1+round((length(cols)-1)*(((head(breaks,-1)+tail(breaks,-1))/2-limits[1])/diff(limits)))
 shape<-sf::st_sfc(isoband::iso_to_sfg(band))
 plot(shape,add=TRUE,col=cols[pmax(1,pmin(length(cols),index))],border=NA)
}
convincing_resample_field<-function(z,old,new){
 if(identical(old,new))return(z)
 tmp<-vapply(seq_along(old),function(j)approx(old,z[,j],new,rule=2)$y,numeric(length(new)))
 t(vapply(seq_along(new),function(i)approx(old,tmp[i,],new,rule=2)$y,numeric(length(new))))
}
convincing_strategy_legend<-function(cols,show_equal=TRUE,compact=FALSE){
 # Approved S14 legend: same colours, wording, ticks and endpoint labels.
 par(mar=c(2.5,.55,.5,1.60),pty='m',xpd=NA);plot.new()
 plot.window(xlim=c(-.08,1),ylim=c(0,1),xaxs='i',yaxs='i')
 zz<-seq(0,1,length.out=801)
 for(k in 1:800)rect(.34,zz[k],.64,zz[k+1],col=cols[k],border=NA)
 rect(.34,0,.64,1,border='black',lwd=.8);wd<-.006
 rect(.34,-wd,.64,wd,col='#08F',lwd=.8,border=1)
 if(show_equal)rect(.34,.5-wd,.64,.5+wd,col='#FF0',lwd=.8,border=1)
 rect(.34,1-wd,.64,1+wd,col='#E06',lwd=.8,border=1)
 ticks<-seq(0,1,.25);segments(.64,ticks,.75,ticks,lwd=.8)
 text(.88,ticks,sprintf('%.2f',ticks),adj=c(0,.5),cex=1.2)
 label<-if(compact)'Training share in tournament 1' else 'Relative productive time in tournament 1'
 text(.16,.5,label,srt=90,cex=if(compact)1.25 else 1.4)
 text(.49,1.035,'1maxxers',cex=1.5,col='#E06')
 text(.49,-.035,'2maxxers',cex=1.5,col='#08F')
}
convincing_scale<-function(cols,limits,label){
 par(mar=c(2.75,.55,.5,1.6),pty='m',xpd=NA);plot.new()
 plot.window(xlim=c(-.08,1),ylim=limits,xaxs='i',yaxs='i')
 z<-seq(limits[1],limits[2],length.out=length(cols)+1)
 for(k in seq_along(cols))rect(.34,z[k],.64,z[k+1],col=cols[k],border=NA)
 rect(.34,limits[1],.64,limits[2],lwd=.8)
 ticks<-seq(limits[1],limits[2],length.out=5);segments(.64,ticks,.75,ticks,lwd=.8)
 text(.88,ticks,sprintf('%.2f',ticks),adj=c(0,.5),cex=1.2)
 text(.16,mean(limits),label,srt=90,cex=1.25)
}
convincing_equal_curves<-function(d,row){
 p<-d$p;a<-d$map$alpha;o<-d$map$allocation
 background<-convincing_background_function(p)
 q<-if(row==1)o[,5]/p$B else o[,1]/pmax(1e-8,o[,1]+o[,2])
 n<-length(a);mask<-matrix(q,n)>=.5
 ij<-which(mask[-1,,drop=FALSE]!=mask[-n,,drop=FALSE],arr.ind=TRUE)
 roots<-vector('list',n)
 if(nrow(ij)){
  lo<-a[ij[,1]];hi<-a[ij[,1]+1];yy<-a[ij[,2]];side<-mask[cbind(ij[,1]+1,ij[,2])]
  for(k in 1:12){mid<-(lo+hi)/2
   if(row==1)above<-finite_pool_alloc(p$baseline,'responsive',cbind(mid,yy))>=p$B/2 else{
    vv<-convincing_optimize(p,cbind(mid,yy));above<-vv[,1]>=vv[,2]}
   movehi<-above==side;hi[movehi]<-mid[movehi];lo[!movehi]<-mid[!movehi]
  }
  for(j in seq_len(n))roots[[j]]<-(lo+hi)[ij[,2]==j]/2
 }
 curves<-list();nr<-max(lengths(roots))
 if(nr)for(k in seq_len(nr))curves[[k]]<-list(x=vapply(roots,function(v)if(length(v)>=k)v[k]else NA_real_,0.),y=a)
 # A jump from one global maximum to another is not an optimal 50:50 split.
 # Check the equal-allocation candidate against the unconstrained global search.
 lapply(curves,function(cl){
  keep<-which(is.finite(cl$x));if(!length(keep))return(cl)
  xy<-cbind(cl$x[keep],cl$y[keep]);opt<-convincing_optimize(p,xy)
  eqvalue<-function(alpha,c){
   m<-background(c)
   x<-(p$B-abs(c))/2;w1<-boredom_w(x,alpha[1],p$beta[1],p$phi[1]);w2<-boredom_w(x,alpha[2],p$beta[2],p$phi[2])
   den1<-(p$N-1)*m[,3]+w1;den2<-(p$N-1)*m[,4]+w2
   ifelse(den1>0,((p$N-1)*m[,1]+x+pmax(c,0))*w1/pmax(den1,1e-300),0)+
    ifelse(den2>0,((p$N-1)*m[,2]+x+pmax(-c,0))*w2/pmax(den2,1e-300),0)
  }
  best<-vapply(seq_len(nrow(xy)),function(i){
   if(row==1)return(eqvalue(xy[i,],0))
   cc<-p$effort
   if(isTRUE(p$aligned))cc<-cc[cc*(2*opt[i,5]-p$B)>= -1e-7]
   v<-eqvalue(xy[i,],cc);ans<-max(v)
   kk<-which(v[2:(length(v)-1)]>=v[1:(length(v)-2)]&v[2:(length(v)-1)]>=v[3:length(v)])+1
   for(k in kk)ans<-max(ans,optimize(function(c)eqvalue(xy[i,],c),cc[c(k-1,k+1)],maximum=TRUE,tol=1e-7)$objective)
   ans
  },0.)
  valid<-best>=opt[,if(row==1)6 else 4]-1e-5*(1+p$B)
  cl$x[keep[!valid]]<-NA_real_;cl
 })
}
convincing_boundary_fields<-function(d,row){
 p<-d$p;a<-d$map$alpha;o<-d$map$allocation;n<-length(a)
 q<-matrix(if(row==1)o[,5]/p$B else o[,1]/pmax(1e-8,o[,1]+o[,2]),n)
 allocate<-if(row==1)function(actor)finite_pool_alloc(p$baseline,'responsive',actor)/p$B else function(actor){
  v<-convincing_optimize(p,actor);v[,1]/pmax(1e-8,v[,1]+v[,2])}
 fields<-lapply(c(1e-5,1-1e-5),function(target){
  mask<-q>=target;ij<-which(mask[-1,,drop=FALSE]!=mask[-n,,drop=FALSE],arr.ind=TRUE)
  roots<-vector('list',n)
  if(nrow(ij)){
   lo<-a[ij[,1]];hi<-a[ij[,1]+1];y<-a[ij[,2]];side<-mask[cbind(ij[,1]+1,ij[,2])]
   for(k in 1:10){mid<-(lo+hi)/2;above<-allocate(cbind(mid,y))>=target;movehi<-above==side;hi[movehi]<-mid[movehi];lo[!movehi]<-mid[!movehi]}
   for(j in seq_len(n))roots[[j]]<-(lo+hi)[ij[,2]==j]/2
  }
  field<-matrix(0,n,n)
  for(j in seq_len(n))field[,j]<-(if(length(roots[[j]]))apply(abs(outer(a,roots[[j]],'-')),1,min)else rep(diff(range(a)),n))*ifelse(mask[,j],1,-1)
  field
 })
 # Ordered scalar fields define disjoint specialist polygons on the same mesh.
 # In particular, two independent contour reconstructions cannot overlap at
 # a direct jump between the two pure strategies.
 list(pmax(fields[[1]],fields[[2]]),pmin(fields[[1]],fields[[2]]))
}
draw_local_global_set<-function(set='asymmetric',png=FALSE,pdf=FALSE,
 file=paste0('Figure_convincing_',set),data=readRDS(file.path('10_convincing_results',paste0(set,'.rds'))),
 sidebars=TRUE,geometry=NULL,
 row_labels=c('Without convincing','With convincing','Convincing effort','Recipient response'),
 effort_label='Signed convincing time / budget',
 response_label='Mean absolute allocation change / budget',
 effort_limit=NULL,response_limit=NULL,onset=NULL,
 entry_lty="26",entry_alpha=.40,entry_lwd=1.1,show_entry=FALSE,
 effort_geometry=NULL,benchmarks=NULL,
 limit_colours=c('#222222','#D89000','#009E73','#8C4CA6'),
 show_maxxer_contours=FALSE,show_equal=TRUE,budget_legend_cex=1.05,
 rows=1:4,column_labels=NULL){
 stopifnot(identical(rows,1:4)||identical(rows,2:3),
  is.null(column_labels)||length(column_labels)==length(data))
 stopifnot(requireNamespace('isoband',quietly=TRUE),requireNamespace('sf',quietly=TRUE))
 nc<-length(data);cols<-colorRampPalette(c('#315A9B','white','#B54343'))(801)
 efforts<-colorRampPalette(c('#08F','white','#E06'))(801);greys<-gray(seq(1,.06,length.out=801))
 cmax<-max(.05,ceiling(max(vapply(data,function(d)max(abs(d$map$allocation[,3]))/d$p$B,0.))*20)/20)
 rmax<-if(4L %in% rows)max(.05,ceiling(max(vapply(data,function(d)max(d$map$recipient_abs)/d$p$B,0.))*20)/20) else 1
 if(!is.null(effort_limit))cmax<-effort_limit
 if(!is.null(response_limit))rmax<-response_limit
 if(is.null(geometry)){
  equal<-lapply(data,function(d)lapply(1:2,function(r)convincing_equal_curves(d,r)))
  fields<-lapply(data,function(d)lapply(1:2,function(r)convincing_boundary_fields(d,r)))
 }else{
  equal<-geometry$equal
  fields<-lapply(seq_along(data),function(j)lapply(geometry$fields[[j]],function(row)
   lapply(row,convincing_resample_field,old=geometry$alpha[[j]],new=data[[j]]$map$alpha)))
 }
 if((png||pdf)&&identical(file,paste0('Figure_convincing_',set))&&set%in%c('symmetric_beta04','symmetric_beta08','asymmetric'))
  saveRDS(list(fields=fields,equal=equal,alpha=lapply(data,function(d)d$map$alpha)),file.path('10_convincing_results',paste0(set,'_geometry.rds')))
 draw<-function(){
  nr<-length(rows)
  lay<-matrix(seq_len(nr*nc),nr,nc,byrow=TRUE)
  legend_ids<-if(nr==4L)c(1,1,2,3)else 1:2
  layout(cbind(lay,nr*nc+legend_ids),widths=c(rep(1,nc),.46))
  par(oma=c(3.5,3.8,3,0))
  labels<-row_labels
  for(r in rows)for(j in seq_len(nc)){
   d<-data[[j]];p<-d$p;a<-d$map$alpha;o<-d$map$allocation;n<-length(a)
   z<-switch(r,o[,5]/p$B,o[,1]/pmax(1e-8,o[,1]+o[,2]),o[,3]/p$B,d$map$recipient_abs/p$B)
   z<-matrix(z,n,n)
   par(mar=c(2.75,2.75,.5,.2),mgp=c(1.75,.48,0),tcl=-.23,pty='s',xpd=FALSE)
   plot(NA,xlim=range(a),ylim=range(a),xaxs='i',yaxs='i',axes=FALSE,xlab='',ylab='',asp=1)
   pal<-if(r<=2)cols else if(r==3)efforts else greys
   ix<-if(r<=2)1+round(800*z)else if(r==3)1+round(400*(1+z/cmax))else 1+round(800*z/rmax)
   cm<-matrix(pal[pmax(1,pmin(801,ix))],n,n)
   if(r==3&&!is.null(effort_geometry[[j]])){
    rect(min(a),min(a),max(a),max(a),col='white',border=NA)
    gg<-effort_geometry[[j]]
    ixg<-1+round(400*(1+((gg$colour_index-401)/400)/cmax))
    if(length(gg$geometry))plot(gg$geometry,col=efforts[pmax(1,pmin(801,ixg))],border=NA,add=TRUE)
   }else convincing_colour_bands(a,z,pal,if(r<=2)c(0,1)else if(r==3)c(-cmax,cmax)else c(0,rmax))
   if(r<=2){
    # Vector polygons cover pure specialists; no pixel outlines.
    convincing_polygon(a,-fields[[j]][[r]][[1]],'#08F');convincing_polygon(a,fields[[j]][[r]][[2]],'#E06')
   }
   abline(0,1,lty=2,lwd=.65)
   for(h in convincing_hdr())lines(h[1,],h[2,],col='#EEE',lwd=.8)
   if(r<=2){
    if(show_maxxer_contours)for(f in fields[[j]][[r]])for(cl in contourLines(a,a,f,levels=0))lines(cl$x,cl$y,lwd=.75,col='black')
    eql<-equal[[j]][[r]]
    if(show_equal)for(cl in eql)lines(cl$x,cl$y,col='black',lwd=3.5)
    if(show_equal)for(cl in eql)lines(cl$x,cl$y,col='#FF0',lwd=2.3)
   }
   box()
   if(show_entry&&r==2&&!is.null(onset))for(cl in onset[[j]])lines(cl$x,cl$y,lty=entry_lty,lwd=entry_lwd,col=adjustcolor('black',alpha.f=entry_alpha))
   if(sidebars&&r<=2){
    budget<-if(r==1)c(p$baseline$responsive$moment[1],p$B-p$baseline$responsive$moment[1],0,0)else
     colSums(p$density$weights*cbind(p$actors[,1:2],pmax(p$actors[,3],0),pmax(-p$actors[,3],0)))
    stopifnot(abs(sum(budget)-p$B)<1e-6)
    ends<-c(0,cumsum(budget)/p$B)*max(a)
    for(k in 1:4)rect(2.54,ends[k],2.69,ends[k+1],col=c('#E06','#08F','#8F004D','#004E85')[k],border=NA,xpd=NA)
    rect(2.54,0,2.69,max(a),lwd=.7,xpd=NA)
   }
   if(sidebars&&r==3){
    campaign<-colSums(p$density$weights*cbind(pmax(p$actors[,3],0),pmax(-p$actors[,3],0)))
    rect(2.54,0,2.69,max(a),col='white',border=NA,xpd=NA)
    if(sum(campaign)>1e-10){
     ends<-c(0,cumsum(campaign)/sum(campaign))*max(a)
     for(k in 1:2)rect(2.54,ends[k],2.69,ends[k+1],col=c('#8F004D','#004E85')[k],border=NA,xpd=NA)
    }
    rect(2.54,0,2.69,max(a),lwd=.7,xpd=NA)
   }
   if(sidebars&&r==4){
    moved<-unname(p$summary['mean_recipient_abs_shift'])/p$B
    stopifnot(is.finite(moved),moved>=0,moved<=1+1e-8)
    rect(2.54,0,2.69,max(a),col='white',border=NA,xpd=NA)
    rect(2.54,0,2.69,max(a)*moved,col='black',border=NA,xpd=NA)
    rect(2.54,0,2.69,max(a),lwd=.7,xpd=NA)
   }
   if(r==2&&j==1){
    budget_y<-if(nr==2L)grconvertY(par('usr')[3],from='user',to='ndc')-1.35*par('csi')/par('din')[2] else
     grconvertY(par('usr')[4],from='user',to='ndc')+1.1*par('csi')/par('din')[2]
    budget_left<-grconvertX(par('usr')[1],from='user',to='ndc')
   }
   if(r==2&&j==nc)budget_right<-grconvertX(par('usr')[2],from='user',to='ndc')
   if(r==3&&!is.null(benchmarks))mtext(LETTERS[j],3,line=.55,at=0,adj=0,cex=1.25,font=2,xpd=NA,col=limit_colours[j])
   if(r==tail(rows,1)){axis(1,at=seq(0,2.5,.5),cex.axis=1.2);mtext(expression(alpha[1]~~'(talent in tournament 1)'),1,line=2.4,cex=.95)}
   if(j==1){axis(2,at=seq(0,2.5,.5),las=1,cex.axis=1.2);mtext(expression(alpha[2]~~'(talent in tournament 2)'),2,line=2.4,cex=.95)}
   if(r==rows[1]){
    if(!is.null(column_labels))mtext(column_labels[j],3,line=.15,at=0,adj=0,cex=.95,xpd=NA)
    mtext(LETTERS[j],3,line=1.8,at=0,adj=0,cex=1.25,font=2,xpd=NA,col='black')
    for(t in 1:2){at<-2.5*c(.48,.88)[t];cc<-c('#E06','#08F')[t]
     mtext(bquote(beta[.(t)]==.(sprintf('%.2f',p$beta[t]))),3,line=1.55,at=at,cex=1.12,col=cc)
     mtext(bquote(varphi[.(t)]==.(sprintf('%.2f',p$phi[t]))),3,line=.15,at=at,cex=1.12,col=cc)
    }
   }
  }
  for(k in seq_along(rows))mtext(labels[rows[k]],2,outer=TRUE,at=1-(k-.5)/nr,line=1.8,cex=1.25)
  convincing_strategy_legend(cols,show_equal=show_equal,compact=nr==2L)
  convincing_scale(efforts,c(-cmax,cmax),effort_label)
  if(!is.null(benchmarks))for(j in seq_len(nc))for(t in 1:2){
   b<-benchmarks[[j]][[t]];local<-data[[j]]$p$N<1e30
   value<-if(local)b$reference_signed else b$normal_signed
   x0<-.27+(j-1)*.48/nc;x1<-.27+j*.48/nc
   if(abs(value)<=cmax)segments(x0,value,x1,value,col=limit_colours[j],lwd=2)
   if(local&&b$low_cross_preferred&&abs(b$low_cross_signed)<=cmax)
    segments(x0,b$low_cross_signed,x1,b$low_cross_signed,col=limit_colours[j],lwd=1.8,lty='22')
  }
  if(4L %in% rows)convincing_scale(greys,c(0,rmax),response_label)
  budget_labels<-c('Mean actor budget:','Training 1','Training 2','Convincing for 1','Convincing for 2')
  budget_cols<-c('black','#E06','#08F','#8F004D','#004E85')
  budget_gp<-grid::gpar(fontsize=par('ps')*budget_legend_cex,fontfamily=par('family'))
  widths<-vapply(budget_labels,function(label)grid::convertWidth(grid::grobWidth(grid::textGrob(label,gp=budget_gp)),'npc',valueOnly=TRUE),0.)
  gap<-.25/par('din')[1];total<-sum(widths)+gap*4
  left<-(budget_left+budget_right-total)/2
  x<-left+cumsum(c(0,head(widths,-1)+gap))+widths/2
  for(k in seq_along(budget_labels))grid::grid.text(budget_labels[k],x=x[k],y=budget_y,
   gp=grid::gpar(col=budget_cols[k],fontsize=par('ps')*budget_legend_cex,fontfamily=par('family')))

 }
 width<-(3580/300)*(nc+.46)/4.46;height<-1740/300*2-length(setdiff(1:4,rows))*2.35
 if(png){grDevices::png(paste0(file,'.png'),width=width,height=height,units='in',res=600,type='cairo');tryCatch(draw(),finally=dev.off())}
 if(pdf){cairo_pdf(paste0(file,'.pdf'),width=width,height=height,family='Arial');tryCatch(draw(),finally=dev.off())}
 if(interactive()||dev.cur()!=1)draw();invisible(data)
}

