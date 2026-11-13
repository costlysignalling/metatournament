source('10_convincing_model.R')
campaign_signed_distance<-function(a,mask,roots){
 n<-length(a);field<-matrix(0,n,n)
 for(j in seq_len(n))field[,j]<-(if(length(roots[[j]]))apply(abs(outer(a,roots[[j]],'-')),1,min)else rep(diff(range(a)),n))*ifelse(mask[,j],1,-1)
 field
}
campaign_curve_roots<-function(a,curves){
 lapply(a,function(y){v<-numeric()
  for(cl in curves){if(length(cl$x)<2)next
   x0<-head(cl$x,-1);x1<-tail(cl$x,-1);y0<-head(cl$y,-1);y1<-tail(cl$y,-1)
   keep<-is.finite(x0+x1+y0+y1)&y1!=y0&pmin(y0,y1)<=y&pmax(y0,y1)>=y
   v<-c(v,x0[keep]+(y-y0[keep])/(y1[keep]-y0[keep])*(x1[keep]-x0[keep]))
  }
  unique(v)
 })
}
campaign_phase_extension<-function(a,value,phase){
 n<-length(a);out<-matrix(NA_real_,n,n)
 good<-which(colSums(phase)>0)
 for(j in good){keep<-which(phase[,j]);out[,j]<-if(length(keep)==1)value[keep,j]else approx(a[keep],value[keep,j],xout=a,rule=2)$y}
 if(!length(good))return(matrix(0,n,n))
 for(j in setdiff(seq_len(n),good))out[,j]<-out[,good[which.min(abs(good-j))]]
 out
}
campaign_isopoly<-function(a,field){
 z<-isoband::isobands(a,a,t(field),0,max(1,max(field)+1))
 sf::st_make_valid(sf::st_sfc(isoband::iso_to_sfg(z)[[1]]))
}

campaign_refine_edges<-function(p,a,c,horizontal=TRUE,iterations=12L,nx=64L){
 z<-if(horizontal)c else t(c);phase<-z>=0;n<-length(a)
 ij<-which(phase[-1,,drop=FALSE]!=phase[-n,,drop=FALSE],arr.ind=TRUE)
 roots<-matrix(NA_real_,n-1,n)
 if(!nrow(ij))return(list(roots=roots,refined=0L,unresolved=0L))
 roots[ij]<-(a[ij[,1]]+a[ij[,1]+1])/2
 pick<-which(abs(z[ij])>1e-4&abs(z[cbind(ij[,1]+1,ij[,2])])>1e-4)
 if(!length(pick))return(list(roots=roots,refined=0L,unresolved=0L))
 lo<-a[ij[pick,1]];hi<-a[ij[pick,1]+1];fixed<-a[ij[pick,2]]
 choice<-function(v)convincing_optimize(p,if(horizontal)cbind(v,fixed)else cbind(fixed,v),nx=nx)[,3]>=0
 left<-choice(lo);side<-choice(hi);valid<-left!=side;bad<-sum(!valid)
 pick<-pick[valid];lo<-lo[valid];hi<-hi[valid];fixed<-fixed[valid];side<-side[valid]
 for(k in if(length(pick))seq_len(iterations)else integer()){
  mid<-(lo+hi)/2;move<-choice(mid)==side;hi[move]<-mid[move];lo[!move]<-mid[!move]
 }
 roots[ij[pick,,drop=FALSE]]<-(lo+hi)/2
 list(roots=roots,refined=length(pick),unresolved=bad)
}
campaign_warped_poly<-function(a,mask,xroots,yroots){
 n<-length(a);step<-a[2]-a[1]
 bands<-isoband::isobands(a,a,t(ifelse(mask,1,-1)),0,2)
 b<-bands[[1]];if(!length(b$x))return(sf::st_sfc(sf::st_geometrycollection()))
 xindex<-(b$x-a[1])/step;yindex<-(b$y-a[1])/step
 hx<-abs(xindex-round(xindex-.5)-.5)<1e-7 & abs(yindex-round(yindex))<1e-7
 ix<-round(xindex-.5)+1L;jx<-round(yindex)+1L
 keep<-which(hx&ix>=1&ix<n&jx>=1&jx<=n)
 if(length(keep)){v<-xroots[cbind(ix[keep],jx[keep])];ok<-is.finite(v);b$x[keep[ok]]<-v[ok]}
 hy<-abs(yindex-round(yindex-.5)-.5)<1e-7 & abs(xindex-round(xindex))<1e-7
 iy<-round(yindex-.5)+1L;jy<-round(xindex)+1L
 keep<-which(hy&iy>=1&iy<n&jy>=1&jy<=n)
 if(length(keep)){v<-yroots[cbind(iy[keep],jy[keep])];ok<-is.finite(v);b$y[keep[ok]]<-v[ok]}
 bands[[1]]<-b;sf::st_make_valid(sf::st_sfc(isoband::iso_to_sfg(bands)[[1]]))
}
campaign_entry_edge_roots<-function(a,mask,curves){
 n<-length(a);roots<-matrix(NA_real_,n-1,n);all<-campaign_curve_roots(a,curves)
 ij<-which(mask[-1,,drop=FALSE]!=mask[-n,,drop=FALSE],arr.ind=TRUE)
 if(nrow(ij))for(k in seq_len(nrow(ij))){i<-ij[k,1];j<-ij[k,2]
  v<-all[[j]];v<-v[v>=a[i]-1e-8&v<=a[i+1]+1e-8]
  roots[i,j]<-if(length(v))v[which.min(abs(v-mean(a[c(i,i+1)])))]else mean(a[c(i,i+1)])
 }
 roots
}

prepare_campaign_effort_geometry<-function(d,onset,iterations=12L,nx=64L,band_width=.005){
 a<-d$map$alpha;n<-length(a);o<-d$map$allocation;p<-d$p
 c<-matrix(o[,3]/p$B,n,n);positive<-c>=0;active<-abs(c)>1e-5
 xr<-campaign_refine_edges(p,a,c,TRUE,iterations,nx)
 yr<-campaign_refine_edges(p,a,c,FALSE,iterations,nx)
 ax<-campaign_entry_edge_roots(a,active,onset)
 reversed<-lapply(onset,function(cl)list(x=cl$y,y=cl$x))
 ay<-campaign_entry_edge_roots(a,t(active),reversed)
 active_poly<-campaign_warped_poly(a,active,ax,ay)
 positive_poly<-campaign_warped_poly(a,positive,xr$roots,yr$roots)
 negative_poly<-campaign_warped_poly(a,!positive,xr$roots,yr$roots)
 shapes<-list();colours<-integer();breaks<-seq(0,1+band_width/2,by=band_width)
 for(side in c(-1,1)){
  phase<-if(side==1)positive else !positive
  if(!any(phase&active))next
  region<-suppressWarnings(sf::st_intersection(if(side==1)positive_poly else negative_poly,active_poly))
  if(all(sf::st_is_empty(region)))next
  magnitude<-campaign_phase_extension(a,abs(c),phase&active)
  iso<-isoband::isobands(a,a,t(magnitude),head(breaks,-1),tail(breaks,-1))
  g<-sf::st_sfc(isoband::iso_to_sfg(iso));keep<-!sf::st_is_empty(g)
  bands<-sf::st_sf(bin=which(keep),geometry=sf::st_make_valid(g[keep]))
  clipped<-suppressWarnings(sf::st_intersection(bands,region))
  if(nrow(clipped)){
   mag<-(breaks[clipped$bin]+breaks[clipped$bin+1])/2
   shapes[[length(shapes)+1]]<-sf::st_geometry(clipped)
   colours<-c(colours,1+round(400*(1+side*mag)))
  }
 }
 geometry<-if(length(shapes))do.call("c",shapes)else sf::st_sfc()
 list(geometry=geometry,colour_index=colours,alpha=a,version=2L,refined_brackets=xr$refined+yr$refined,unresolved_brackets=xr$unresolved+yr$unresolved,
  max_bracket_width=diff(a)[1]/2^iterations,band_width=band_width,
  phase_geometry=positive_poly)
}
prepare_campaign_effort_geometries<-function(overwrite=FALSE){
 for(set in c('symmetric_beta04','symmetric_beta08','asymmetric')){
  data<-readRDS(paste0('10_convincing_results/',set,'.rds'))
  entry<-readRDS(paste0('10_local_global_results/entry_local_',set,'.rds'))
  for(j in seq_along(data)){
   file<-paste0('10_local_global_results/effort_geometry_local_',set,'_',j,'.rds')
   if(!overwrite&&file.exists(file)&&identical(readRDS(file)$version,2L))next
   message('Refine campaign directions ',set,' ',j)
   g<-prepare_campaign_effort_geometry(data[[j]],entry[[j]])
   saveRDS(g,file);message('  refined ',g$refined_brackets,' brackets; polygons ',length(g$geometry))
  }
 }
 invisible(TRUE)
}
