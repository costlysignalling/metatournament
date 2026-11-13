# Script 06: beta-space strategy coexistence within the 99% talent HDR.
# Fixed E1=E2: actors respond to uniform random allocations, with M1 and M2
# calculated separately. Responsive E_t: continuous-population equilibrium.
# Geometry uses global optima and positive-area HDR intersections, not QMC counts.
# Supplied score and geometry caches permit redraws without a compiler.
# Score cache: 06_coexistence_beta_inputs.rds.
# Delete that file to rebuild with 06_prepare_coexistence_deterministic.R;
# the numerical calculation is checkpointed and requires Rcpp/Rtools.
if(!requireNamespace('sf',quietly=TRUE))stop("Package 'sf' is required.")
source('06_prepare_coexistence_deterministic.R')
coex<-prepare_coexistence('beta')
source('06_refine_coexistence_curves.R')
refined_geometry<-refine_coexistence_curves('beta')$geometry
B<-coex$B;div_tol<-coex$levels[2]
phi_vals<-coex$panel_vals;betas<-coex$axis

source('05_parameter_space_example_markers.R')
example_sets<-parameter_space_example_vectors('beta',phi_vals,range(betas))
example_set_cols<-parameter_space_marker_colours
example_set_labs<-parameter_space_marker_labels

example_pch<-4
example_cex<-1.30
example_lwd<-2.5


bits<-c(1L,2L,4L,8L)
levels_presence<-c(0,div_tol,div_tol,0)

# Region fills. These correspond to the 4-bit strategy-presence code.
fill_map<-c(
  "#F2F2F2", # none / background
  "#AFC7EE", # 2max only
  "#7FAEE6", # 2-lean only
  "#A9D2FF", # 2max + 2-lean
  "#E6B0B0", # 1-lean only
  "#D8C7C5", # 2max + 1-lean
  "#D9D9D9", # both leaners
  "#CFE2F8", # 2max + both leaners
  "#F2A9C2", # 1max only
  "#CBB6E7", # 2max + 1max
  "#E4C0D6", # 2-lean + 1max
  "#C2C7F1", # 2max + 2-lean + 1max
  "#F6C2CF", # 1-lean + 1max
  "#E8C1D6", # 2max + 1-lean + 1max
  "#F0D7D7", # both leaners + 1max
  "#CFE8C7"  # all four
)

cat_cols<-c("#08F","#6F8FCF","#D27A7A","#E06")
cat_labs<-c(
  "2maxxers",
  "2-leaning diversifiers",
  "1-leaning diversifiers",
  "1maxxers"
)

line_cols<-c("#005BBB","#3E7CB1","#C06C84","#D0006F")
line_lty<-c(1,2,2,1)
line_labs<-c(
  "2maxxer presence boundary",
  "2-leaning presence boundary",
  "1-leaning presence boundary",
  "1maxxer presence boundary"
)

normalize_example_sets<-function(example_sets,phi_vals,tol=1e-9){
  
  if(length(example_sets)==0) return(data.frame())
  
  one_row<-function(v,set_name){
    vv<-v
    scope<-if(is.null(names(vv)) || !("scope" %in% names(vv))){
      "both"
    } else {
      as.character(vv[["scope"]])
    }
    
    if(is.null(names(vv))) stop("Each example combination must be named.")
    
    need<-c("beta1","beta2","phi1","phi2")
    miss<-setdiff(need,names(vv))
    if(length(miss)){
      stop("Missing fields in example combination: ",
           paste(miss,collapse=", "))
    }
    
    data.frame(
      beta1=as.numeric(vv[["beta1"]]),
      beta2=as.numeric(vv[["beta2"]]),
      phi1 =as.numeric(vv[["phi1"]]),
      phi2 =as.numeric(vv[["phi2"]]),
      scope=scope,
      set  =set_name,
      stringsAsFactors=FALSE
    )
  }
  
  rows<-list()
  
  for(set_name in names(example_sets)){
    sets_i<-example_sets[[set_name]]
    if(length(sets_i)==0) next
    
    rows_i<-lapply(sets_i,one_row,set_name=set_name)
    rows<-c(rows,rows_i)
  }
  
  if(length(rows)==0) return(data.frame())
  
  ex<-do.call(rbind,rows)
  
  ex$col<-match(round(ex$phi1,10),round(phi_vals,10))
  ex$row<-match(round(ex$phi2,10),round(phi_vals,10))
  
  bad<-is.na(ex$col) | is.na(ex$row)
  if(any(bad)){
    warning(
      "Dropping example combinations whose phi values are not in phi_vals: ",
      paste(
        apply(ex[bad,c("phi1","phi2"),drop=FALSE],1,paste,collapse=","),
        collapse="; "
      )
    )
    ex<-ex[!bad,,drop=FALSE]
  }
  
  if(nrow(ex)){
    ex$scope<-tolower(ex$scope)
    ex$scope[!(ex$scope %in% c("fixed","responsive","both"))]<-"both"
    ex$col_marker<-unname(example_set_cols[ex$set])
    ex$lab_marker<-unname(example_set_labs[ex$set])
  }
  
  rownames(ex)<-NULL
  ex
}

example_df<-normalize_example_sets(example_sets,phi_vals)

if(nrow(example_df)){
  message("\nExample combinations that will be plotted:")
  print(
    example_df[,c("set","beta1","beta2","phi1","phi2","row","col","scope")],
    row.names=FALSE
  )
}

# Exchange tournaments 1 <-> 2, including beta axes.
swap_panel<-function(x){
  y<-list(
    mode=x$mode,
    p1=x$p2,p2=x$p1,
    S2max=t(x$S1max),
    S2div=t(x$S1div),
    S1div=t(x$S2div),
    S1max=t(x$S2max)
  )

  if(!is.null(x$rho)){
    y$rho<-1/t(x$rho)
    y$E1share<-1-t(x$E1share)
  }

  y
}

# r indexes phi2 (row); c indexes phi1 (column).
load_panel<-function(mode,r,c){
 p<-coex$panels[[c,r]]
 out<-c(list(mode=mode,p1=phi_vals[c],p2=phi_vals[r]),p[[mode]])
 if(mode=='responsive'){out$rho<-p$rho;out$E1share<-p$E1share}
 out
}

interp_grid<-function(z,x,y){
  i<-pmax(1L,pmin(length(betas)-1L,findInterval(x,betas)))
  j<-pmax(1L,pmin(length(betas)-1L,findInterval(y,betas)))

  x0<-betas[i];x1<-betas[i+1L]
  y0<-betas[j];y1<-betas[j+1L]
  tx<-(x-x0)/(x1-x0)
  ty<-(y-y0)/(y1-y0)

  (1-tx)*(1-ty)*z[cbind(i,j)] +
    tx*(1-ty)*z[cbind(i+1L,j)] +
    (1-tx)*ty*z[cbind(i,j+1L)] +
    tx*ty*z[cbind(i+1L,j+1L)]
}

contour_family<-function(z,level){
  contourLines(x=betas,y=betas,z=z,levels=level)
}

line_to_sfg<-function(x){
  xy<-cbind(x$x,x$y)

  # Snap endpoints that are numerically indistinguishable from the square edge.
  xy[abs(xy[,1]-.01)<1e-10,1]<-.01
  xy[abs(xy[,1]-.99)<1e-10,1]<-.99
  xy[abs(xy[,2]-.01)<1e-10,2]<-.01
  xy[abs(xy[,2]-.99)<1e-10,2]<-.99

  sf::st_linestring(xy)
}

make_panel_geometry<-function(panel){
  fields<-list(panel$S2max,panel$S2div,panel$S1div,panel$S1max)
  clines<-Map(contour_family,fields,levels_presence)

  square<-sf::st_linestring(
    rbind(c(.01,.01),c(.99,.01),c(.99,.99),c(.01,.99),c(.01,.01))
  )

  all_sfg<-list(square)
  for(fam in clines){
    if(length(fam)) all_sfg<-c(all_sfg,lapply(fam,line_to_sfg))
  }

  network<-sf::st_sfc(all_sfg)
  network<-sf::st_node(sf::st_union(network))
  polygeom<-sf::st_collection_extract(sf::st_polygonize(network),"POLYGON")

  if(length(polygeom)==0) stop("Polygonization failed.")

  polys<-sf::st_sf(geometry=polygeom)
  rp<-sf::st_point_on_surface(polys)
  xy<-sf::st_coordinates(rp)

  vals<-vapply(
    fields,
    function(z) interp_grid(z,xy[,1],xy[,2]),
    numeric(nrow(xy))
  )

  present<-sweep(vals,2,levels_presence,`>`)
  polys$rid<-rowSums(present * rep(bits,each=nrow(present)))

  ids<-sort(unique(polys$rid))
  dissolved<-do.call(
    rbind,
    lapply(ids,function(id){
      sf::st_sf(
        rid=id,
        geometry=sf::st_union(polys$geometry[polys$rid==id])
      )
    })
  )

  list(clines=clines,dissolved=dissolved)
}

# Approximate pole of inaccessibility: visually central point inside each region.
visual_centre_one<-function(poly,n=31,refine=5){
  geom<-sf::st_geometry(poly)
  bb<-sf::st_bbox(geom)

  xmin<-unname(bb["xmin"]);xmax<-unname(bb["xmax"])
  ymin<-unname(bb["ymin"]);ymax<-unname(bb["ymax"])

  best<-sf::st_coordinates(sf::st_point_on_surface(sf::st_geometry(poly)))[1,1:2]

  for(k in seq_len(refine)){
    xs<-seq(xmin,xmax,length.out=n)
    ys<-seq(ymin,ymax,length.out=n)
    gr<-expand.grid(x=xs,y=ys)

    pts<-sf::st_as_sf(gr,coords=c("x","y"),crs=sf::st_crs(poly))
    inside<-lengths(sf::st_within(pts,geom))>0

    if(any(inside)){
      pin<-pts[inside,,drop=FALSE]
      d<-as.numeric(sf::st_distance(pin,sf::st_boundary(geom)))
      imax<-which.max(d)
      if(length(imax) && is.finite(d[imax])){
        best<-sf::st_coordinates(pin[imax,])[1,1:2]
      }
    }

    dx<-(xmax-xmin)/(n-1)
    dy<-(ymax-ymin)/(n-1)

    xmin<-max(unname(bb["xmin"]),best[1]-2*dx)
    xmax<-min(unname(bb["xmax"]),best[1]+2*dx)
    ymin<-max(unname(bb["ymin"]),best[2]-2*dy)
    ymax<-min(unname(bb["ymax"]),best[2]+2*dy)
  }

  best
}

region_centres<-function(gg,n=31,refine=5){
  if(nrow(gg)==0) return(matrix(numeric(0),ncol=2))

  t(vapply(
    seq_len(nrow(gg)),
    function(i) visual_centre_one(gg[i,,drop=FALSE],n=n,refine=refine),
    numeric(2)
  ))
}

draw_panel<-function(panel,r,c,block){
  
  geom<-refined_geometry[[block]][[c,r]]
  
  show_x<-r==3
  show_y<-c==1
  
  phi_cex  <- 1.18
  axis_cex <- 1.02
  lab_cex  <- 1.18
  
  par(
    mar=c(1.50,1.50,0.5,0.5),
    mgp=c(1.35,.36,0),
    tcl=-.24,
    xaxs="i",yaxs="i",
    pty="s",
    xpd=FALSE
  )
  
  plot(
    NA,
    xlim=c(.01,.99),ylim=c(.01,.99),asp=1,
    xlab="",ylab="",axes=FALSE
  )
  
  # filled coexistence regions
  for(i in seq_len(nrow(geom$dissolved))){
    plot(
      sf::st_geometry(geom$dissolved[i,]),
      add=TRUE,
      col=fill_map[geom$dissolved$rid[i]+1L],
      border=NA
    )
  }
  
  # presence boundaries
  for(k in seq_along(geom$clines)){
    if(length(geom$clines[[k]])){
      for(seg in geom$clines[[k]]){
        lines(
          seg$x,seg$y,
          col=line_cols[k],
          lty=line_lty[k],
          lwd=1.55
        )
      }
    }
  }
  
  # strategy-symbol clusters
  for(i in seq_len(nrow(geom$dissolved))){
    gg<-suppressWarnings(sf::st_cast(geom$dissolved[i,],"POLYGON"))
    if(nrow(gg)==0) next
    
    aa<-as.numeric(sf::st_area(gg))
    gg<-gg[aa>.006,,drop=FALSE]
    if(nrow(gg)==0) next
    
    pts<-region_centres(gg)
    rid<-geom$dissolved$rid[i]
    
    for(j in seq_len(nrow(pts))){
      cx<-pts[j,1]
      cy<-pts[j,2]
      d<-.021
      
      pos<-rbind(
        c(cx,   cy+d),
        c(cx-d, cy),
        c(cx+d, cy),
        c(cx,   cy-d)
      )
      
      for(h in seq_along(bits)){
        if(bitwAnd(rid,bits[h])>0){
          points(
            pos[h,1],pos[h,2],
            pch=21,bg=cat_cols[h],col="black",
            lwd=.85,cex=.92
          )
        }
      }
    }
  }
  
  # example parameter combinations
  # example parameter combinations
  if(nrow(example_df)){
    keep<-example_df$row==r &
      example_df$col==c &
      (example_df$scope=="both" | example_df$scope==block)
    
    ex<-example_df[keep,,drop=FALSE]
    
    if(nrow(ex)){
      for(ii in seq_len(nrow(ex))){
        # A thin black underlay makes both yellow and turquoise crosses
        # visible over any region fill or boundary line.
        points(
          ex$beta1[ii],ex$beta2[ii],
          pch=example_pch,
          col="black",
          lwd=example_lwd+1.0,
          cex=example_cex
        )
        points(
          ex$beta1[ii],ex$beta2[ii],
          pch=example_pch,
          col=ex$col_marker[ii],
          lwd=example_lwd,
          cex=example_cex
        )
      }
    }
  }
  
  at<-c(.2,.4,.6,.8)
  
  # x axis
  if(show_x){
    axis(1,at=at,cex.axis=axis_cex,lwd=.8,lwd.ticks=.8)
    mtext(expression(beta[1]),side=1,line=1.90,cex=lab_cex*0.9,xpd=NA)
  }
  
  # y axis
  if(show_y){
    axis(2,las=0,at=at,cex.axis=axis_cex,lwd=.8,lwd.ticks=.8)
    mtext(expression(beta[2]),side=2,line=1.4,cex=lab_cex*0.9,las=1)
    
    usr<-par("usr")
    text(
      x=usr[1]-.36,
      y=mean(usr[3:4]),
      labels=bquote(varphi[2]==.(formatC(phi_vals[r],digits=2,format="f"))),
      srt=90,
      xpd=NA,
      cex=phi_cex*1.5
    )
  }
  
  # top phi1 labels
  if(r==1){
    mtext(
      bquote(varphi[1]==.(formatC(phi_vals[c],digits=2,format="f"))),
      side=3,line=.22,cex=phi_cex
    )
  }
  
  box(lwd=.8)
}


# Horizontal legend measured in PHYSICAL inches rather than user coordinates.
# This makes the spacing independent of the layout panel and output dimensions.
draw_point_legend_row<-function(labels,pch,bg,col,lwd,
                                y,cex=1.5,pt.cex=1.30,
                                point_text_gap=.10,item_gap=.24,
                                centre){
  
  text_w<-strwidth(labels,cex=cex,units="inches")
  total_w<-sum(text_w)+length(labels)*point_text_gap+
    (length(labels)-1)*item_gap
  
  x<-centre-total_w/2
  
  for(i in seq_along(labels)){
    
    # Crosses get the same black underlay as in the panels.
    if(pch[i]==example_pch){
      points(
        x,y,
        pch=pch[i],
        col="black",
        cex=pt.cex,
        lwd=lwd[i]+1
      )
    }
    
    points(
      x,y,
      pch=pch[i],
      bg=bg[i],
      col=col[i],
      cex=pt.cex,
      lwd=lwd[i]
    )
    
    tx<-x+point_text_gap
    
    text(
      tx,y,
      labels[i],
      adj=c(0,.5),
      cex=cex
    )
    
    x<-tx+text_w[i]
    
    if(i<length(labels)){
      x<-x+item_gap
    }
  }
}


draw_line_legend_row<-function(labels,col,lty,lwd,
                               y,cex=1.5,
                               seg_width=.28,
                               line_text_gap=.10,item_gap=.24,
                               centre){

  text_w<-strwidth(labels,cex=cex,units="inches")
  total_w<-length(labels)*(seg_width+line_text_gap)+
    sum(text_w)+(length(labels)-1)*item_gap

  x<-centre-total_w/2

  for(i in seq_along(labels)){
    segments(
      x,y,
      x+seg_width,y,
      col=col[i],
      lty=lty[i],
      lwd=lwd[i]
    )

    tx<-x+seg_width+line_text_gap
    text(tx,y,labels[i],adj=c(0,.5),cex=cex)

    x<-tx+text_w[i]
    if(i<length(labels)) x<-x+item_gap
  }
}


draw_common_legend<-function(){

  par(
    mar=c(.03,.03,.03,.03),
    xpd=NA,pty="m"
  )

  plot.new()

  # Use physical inches as the x-coordinate system. strwidth(...,
  # units="inches") is therefore directly commensurate with the plotting
  # coordinates and the whole row can be centred exactly.
  legend_width<-par("pin")[1]
  plot.window(
    xlim=c(0,legend_width),
    ylim=c(0,1),
    xaxs="i",yaxs="i"
  )

  centre<-legend_width/2

  # Separate rows keep the strategy symbols and source-figure examples readable.
  draw_point_legend_row(labels=cat_labs,pch=rep(21,4),bg=cat_cols,
    col=rep('black',4),lwd=rep(.95,4),y=.49,cex=1.5,pt.cex=1.30,
    point_text_gap=.10,item_gap=.24,centre=centre)
  draw_parameter_example_legend(example_sets,y=.80,centre=centre)

  # Second row: four presence-boundary lines.
  draw_line_legend_row(
    labels=line_labs,
    col=line_cols,
    lty=line_lty,
    lwd=rep(2.1,length(line_labs)),
    y=.15,
    cex=1.5,
    seg_width=.28,
    line_text_gap=.10,
    item_gap=.24,
    centre=centre
  )
}

plot_combined<-function(){
  
  lay<-rbind(
    c( 1, 2, 3, 0,10,11,12),
    c( 4, 5, 6, 0,13,14,15),
    c( 7, 8, 9, 0,16,17,18),
    c(19,19,19,19,19,19,19)
  )
  
  layout(
    lay,
    widths =c(1,1,1,.32,1,1,1),
    heights=c(1,1,1,.95),
    respect=FALSE
  )
  
  par(
    oma=c(.45,3.35,6.55,.45)
  )
  
  # fixed block
  for(r in 1:3){
    for(c in 1:3){
      draw_panel(
        load_panel("fixed",r,c),
        r,c,"fixed"
      )
    }
  }
  
  # responsive block
  for(r in 1:3){
    for(c in 1:3){
      draw_panel(
        load_panel("responsive",r,c),
        r,c,"responsive"
      )
    }
  }
  
  # common legend
  draw_common_legend()
  
  # outer headings
  mtext(
    expression("Fixed "*E[1]==E[2]),
    side=3,outer=TRUE,
    at=.225,line=3.95,
    cex=1.55,font=2
  )
  mtext(
    "What people prefer to do if tournament interests of others are equal",
    side=3,outer=TRUE,
    at=.225,line=2.25,
    cex=1.1,xpd=NA
  )
  
  mtext(
    expression("Responsive "*E[t]),
    side=3,outer=TRUE,
    at=.775,line=3.95,
    cex=1.55,font=2
  )
  mtext(
    "What people end up doing in unequal-interest equilibrium",
    side=3,outer=TRUE,
    at=.775,line=2.25,
    cex=1.1
  )
}



# Identical physical layout; 600 dpi PNG and matching PDF.
figure_width<-3500/300
reference_height<-figure_width*1900/3100
reference_row_height<-(reference_height-7*.132)/3.56
figure_height<-reference_height+(.95-.56)*reference_row_height
png('Figure_strategy_coexistence_fixed_vs_responsive.png',width=figure_width,height=figure_height,units='in',res=600)
tryCatch(plot_combined(),finally=dev.off())
# Retain the original script's alternate output name for existing links.
invisible(file.copy('Figure_strategy_coexistence_fixed_vs_responsive.png','strategy_coexistence_fixed_vs_responsive.png',overwrite=TRUE))
cairo_pdf('Figure_strategy_coexistence_fixed_vs_responsive.pdf',width=figure_width,height=figure_height,family='Arial')
tryCatch(plot_combined(),finally=dev.off())
if(interactive())plot_combined()
