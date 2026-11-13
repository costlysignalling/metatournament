# contour2.R
#
# A reproducible, pure-R variant of graphics::contour() with one extra
# argument, label.shift, for method = "flattest".
#
# label.shift is a fraction of the visible connected contour length in the
# current plotting panel.  For example, label.shift = 0.10 moves the label
# 10% of that visible contour length forward along the contour as the target;
# the nearest valid base-style label gap is used. Negative values move it in
# the opposite direction.
#
# The flattest-gap search and gap-drawing logic are adapted from the R 4.6.1
# graphics contour implementation (src/library/graphics/src/plot3d.c).
# R is distributed under GPL-2-or-later; this derived file is likewise
# distributed under GPL-2-or-later.

.contour2_to_inches <- function(p) {
  cbind(
    x = grconvertX(p[, 1L], from = "user", to = "inches"),
    y = grconvertY(p[, 2L], from = "user", to = "inches")
  )
}

.contour2_arc <- function(p) {
  pi <- .contour2_to_inches(p)
  if (nrow(pi) < 2L)
    return(list(p.in = pi, s = 0, length = 0))
  ds <- sqrt(rowSums((pi[-1L, , drop = FALSE] -
                        pi[-nrow(pi), , drop = FALSE])^2))
  s <- c(0, cumsum(ds))
  list(p.in = pi, s = s, length = s[length(s)])
}

# Liang-Barsky clipping of one line segment to a rectangular plot region.
.contour2_clip_segment <- function(a, b, rect) {
  dx <- b[1L] - a[1L]
  dy <- b[2L] - a[2L]
  p <- c(-dx, dx, -dy, dy)
  q <- c(a[1L] - rect[1L], rect[2L] - a[1L],
         a[2L] - rect[3L], rect[4L] - a[2L])
  t0 <- 0
  t1 <- 1
  
  for (k in 1:4) {
    if (abs(p[k]) <= .Machine$double.eps) {
      if (q[k] < 0)
        return(NULL)
    } else {
      r <- q[k] / p[k]
      if (p[k] < 0) {
        if (r > t1)
          return(NULL)
        if (r > t0)
          t0 <- r
      } else {
        if (r < t0)
          return(NULL)
        if (r < t1)
          t1 <- r
      }
    }
  }
  
  list(
    a = a + t0 * c(dx, dy),
    b = a + t1 * c(dx, dy)
  )
}

# Clip a contour polyline to the current panel.  A contour that leaves and
# re-enters the panel may become several visible connected runs.
.contour2_clip_polyline <- function(p, rect) {
  if (nrow(p) < 2L)
    return(list())
  
  out <- list()
  current <- NULL
  tol <- 1e-10 * max(abs(rect), diff(rect[1:2]), diff(rect[3:4]), 1)
  
  flush <- function() {
    if (!is.null(current) && nrow(current) >= 2L)
      out[[length(out) + 1L]] <<- current
    current <<- NULL
  }
  
  for (i in seq_len(nrow(p) - 1L)) {
    z <- .contour2_clip_segment(p[i, ], p[i + 1L, ], rect)
    if (is.null(z)) {
      flush()
      next
    }
    
    if (is.null(current)) {
      current <- rbind(z$a, z$b)
    } else {
      last <- current[nrow(current), ]
      if (sqrt(sum((last - z$a)^2)) <= tol) {
        if (sqrt(sum((last - z$b)^2)) > tol)
          current <- rbind(current, z$b)
      } else {
        flush()
        current <- rbind(z$a, z$b)
      }
    }
  }
  flush()
  out
}

# Port of R's FindCorners() geometry.  The returned 4x2 matrix is used for
# rejecting labels outside the panel and labels intersecting earlier labels.
.contour2_find_corners <- function(width, height, p0, p1, p0.in, p1.in) {
  delta <- height / width
  
  # The original C code converts the x displacement to inches and then
  # converts that amount through the y scale, and vice versa.
  yoff <- yinch((p1.in[1L] - p0.in[1L]) * delta, warn.log = FALSE)
  xoff <- xinch((p1.in[2L] - p0.in[2L]) * delta, warn.log = FALSE)
  
  rbind(
    c(p0[1L] + xoff, p0[2L] - yoff),
    c(p0[1L] - xoff, p0[2L] + yoff),
    c(p1[1L] - xoff, p1[2L] + yoff),
    c(p1[1L] + xoff, p1[2L] - yoff)
  )
}

.contour2_label_inside_panel <- function(corners, rect) {
  all(corners[, 1L] >= rect[1L] & corners[, 1L] <= rect[2L] &
        corners[, 2L] >= rect[3L] & corners[, 2L] <= rect[4L])
}

# Same edge/edge test used by R's TestLabelIntersection().
.contour2_labels_intersect <- function(a, b) {
  for (i in 1:4) {
    i2 <- if (i == 4L) 1L else i + 1L
    Ax <- a[i, 1L]; Ay <- a[i, 2L]
    Bx <- a[i2, 1L]; By <- a[i2, 2L]
    
    for (j in 1:4) {
      j2 <- if (j == 4L) 1L else j + 1L
      ax <- b[j, 1L]; ay <- b[j, 2L]
      bx <- b[j2, 1L]; by <- b[j2, 2L]
      
      dom <- Bx * by - Bx * ay - Ax * by + Ax * ay -
        bx * By + bx * Ay + ax * By - ax * Ay
      
      if (dom == 0) {
        result1 <- -1
        result2 <- -1
      } else {
        result1 <- (bx * Ay - ax * Ay - ay * bx - Ax * by +
                      Ax * ay + by * ax) / dom
        
        if (bx - ax == 0) {
          if (by - ay == 0)
            result2 <- -1
          else
            result2 <- (Ay + (By - Ay) * result1 - ay) / (by - ay)
        } else {
          result2 <- (Ax + (Bx - Ax) * result1 - ax) / (bx - ax)
        }
      }
      
      if (result1 >= 0 && result1 <= 1 &&
          result2 >= 0 && result2 <= 1)
        return(TRUE)
    }
  }
  FALSE
}

.contour2_hits_old_label <- function(corners, label.boxes) {
  if (!length(label.boxes))
    return(FALSE)
  any(vapply(label.boxes,
             function(z) .contour2_labels_intersect(corners, z),
             logical(1)))
}

.contour2_edge_distance <- function(p, rect) {
  min(p[1L] - rect[1L], rect[2L] - p[1L],
      p[2L] - rect[3L], rect[4L] - p[2L])
}

# Enumerate the same sort of label-sized gaps considered by method="flattest".
# Variance deliberately follows the R 4.6.1 C implementation, including its
# cumulative avgGradient update.
.contour2_candidates <- function(p, label.width, label.height,
                                 label.boxes, rect) {
  ns <- nrow(p)
  if (ns < 3L)
    return(list())
  
  ai <- .contour2_arc(p)
  pin <- ai$p.in
  s <- ai$s
  ans <- list()
  
  # R's C code uses iii=0,... and jjj < ns-1.  In 1-based R indexing this
  # means a candidate start up to ns-2 and an endpoint no later than ns-1.
  if (ns < 4L)
    return(ans)
  
  for (iii in seq_len(ns - 2L)) {
    distance.sum <- 0
    avg.gradient <- 0
    square.sum <- 0
    n <- 0L
    jjj <- iii + 1L
    
    while (jjj <= ns - 1L && distance.sum < label.width) {
      d <- pin[jjj, ] - pin[iii, ]
      distance.sum <- sqrt(sum(d^2))
      
      dx <- p[jjj, 1L] - p[jjj - 1L, 1L]
      dy <- p[jjj, 2L] - p[jjj - 1L, 2L]
      if (dx == 0)
        dx <- 1
      avg.gradient <- avg.gradient + dy / dx
      square.sum <- square.sum + avg.gradient * avg.gradient
      
      jjj <- jjj + 1L
      n <- n + 1L
    }
    
    if (distance.sum < label.width)
      break
    
    jj <- iii + n
    if (jj > ns)
      next
    
    corners <- .contour2_find_corners(
      label.width, label.height,
      p[iii, ], p[jj, ], pin[iii, ], pin[jj, ]
    )
    
    if (!.contour2_label_inside_panel(corners, rect))
      next
    if (.contour2_hits_old_label(corners, label.boxes))
      next
    
    variance <- (square.sum - avg.gradient * avg.gradient / n) / n
    avg.gradient <- avg.gradient / n  # retained for correspondence to C
    
    dxy <- sqrt(sum((pin[jj, ] - pin[iii, ])^2))
    if (!is.finite(dxy) || dxy <= 0)
      next
    
    # R fills any excess gap from the side farther from the plot edge,
    # leaving the actual text-sized segment at the nearer edge.
    closest.end <- .contour2_edge_distance(p[jj, ], rect) <=
      .contour2_edge_distance(p[iii, ], rect)
    frac <- min(1, label.width / dxy)
    
    centre.s <- if (closest.end) {
      s[jj] - 0.5 * frac * (s[jj] - s[iii])
    } else {
      s[iii] + 0.5 * frac * (s[jj] - s[iii])
    }
    
    ans[[length(ans) + 1L]] <- list(
      i = iii,
      j = jj,
      range = n,
      variance = variance,
      corners = corners,
      centre.s = centre.s,
      dxy = dxy,
      closest.end = closest.end
    )
  }
  
  attr(ans, "arc") <- ai
  ans
}

.contour2_draw_run <- function(p, candidate = NULL, label = NULL,
                               label.width = NULL, labcex = 0.6,
                               col = par("fg"), lty = par("lty"),
                               lwd = par("lwd"), vfont = NULL) {
  ns <- nrow(p)
  
  if (is.null(candidate)) {
    if (ns >= 2L)
      lines(p[, 1L], p[, 2L], col = col, lty = lty, lwd = lwd)
    return(invisible())
  }
  
  i <- candidate$i
  j <- candidate$j
  
  # Draw the two contour pieces, leaving a real gap.  Nothing is painted
  # over the contour, so this cannot create a white/background-colour nugget.
  if (i > 1L)
    lines(p[seq_len(i), 1L], p[seq_len(i), 2L],
          col = col, lty = lty, lwd = lwd)
  if (j < ns)
    lines(p[j:ns, 1L], p[j:ns, 2L],
          col = col, lty = lty, lwd = lwd)
  
  p0 <- p[i, ]
  p1 <- p[j, ]
  frac <- min(1, label.width / candidate$dxy)
  
  if (candidate$closest.end) {
    xstart <- p1 - (p1 - p0) * frac
    if (frac < 1)
      segments(p0[1L], p0[2L], xstart[1L], xstart[2L],
               col = col, lty = lty, lwd = lwd)
    u <- xstart
    v <- p1
  } else {
    xstart <- p0 + (p1 - p0) * frac
    if (frac < 1)
      segments(xstart[1L], xstart[2L], p1[1L], p1[2L],
               col = col, lty = lty, lwd = lwd)
    u <- p0
    v <- xstart
  }
  
  # R orients the text from left to right where possible.
  if (u[1L] > v[1L]) {
    tmp <- u
    u <- v
    v <- tmp
  }
  
  ui <- c(grconvertX(u[1L], "user", "inches"),
          grconvertY(u[2L], "user", "inches"))
  vi <- c(grconvertX(v[1L], "user", "inches"),
          grconvertY(v[2L], "user", "inches"))
  
  # R 4.6.1 literally uses 180/3.14 here; retaining that makes this port
  # track the base implementation more closely than using 180/pi.
  angle <- (180 / 3.14) * atan2(vi[2L] - ui[2L], vi[1L] - ui[1L])
  
  txt.args <- list(
    x = u[1L], y = u[2L], labels = label,
    adj = c(0, 0.5), srt = angle,
    cex = labcex, col = col
  )
  if (!is.null(vfont))
    txt.args$vfont <- vfont
  do.call(text, txt.args)
  
  invisible()
}

#' Base contour plot with shifted flattest labels
#'
#' This is intentionally close to graphics::contour.default(), with one new
#' argument: label.shift.  The custom code is used only when labels are drawn,
#' method="flattest", and label.shift is non-zero.  Otherwise the call is
#' delegated to graphics::contour.default().
contour2 <- function(x = seq(0, 1, length.out = nrow(z)),
                     y = seq(0, 1, length.out = ncol(z)),
                     z,
                     nlevels = 10,
                     levels = pretty(zlim, nlevels),
                     labels = NULL,
                     xlim = range(x, finite = TRUE),
                     ylim = range(y, finite = TRUE),
                     zlim = range(z, finite = TRUE),
                     labcex = 0.6,
                     drawlabels = TRUE,
                     method = "flattest",
                     label.shift = 0,
                     vfont,
                     axes = TRUE,
                     frame.plot = axes,
                     col = par("fg"),
                     lty = par("lty"),
                     lwd = par("lwd"),
                     add = FALSE,
                     ...) {
  
  if (missing(z)) {
    if (!missing(x)) {
      if (is.list(x)) {
        z <- x$z
        y <- x$y
        x <- x$x
      } else {
        z <- x
        x <- seq.int(0, 1, length.out = nrow(z))
      }
    } else {
      stop("no 'z' matrix specified")
    }
  } else if (is.list(x)) {
    y <- x$y
    x <- x$x
  }
  
  if (any(diff(x) <= 0) || any(diff(y) <= 0))
    stop("increasing 'x' and 'y' values expected")
  if (!is.matrix(z) || nrow(z) <= 1L || ncol(z) <= 1L)
    stop("no proper 'z' matrix specified")
  
  method.no <- pmatch(method[1L], c("simple", "edge", "flattest"))
  if (is.na(method.no))
    stop("invalid 'method'")
  method <- c("simple", "edge", "flattest")[method.no]
  
  have.vfont <- !missing(vfont)
  
  # Exact base-R fallback whenever the new feature is not needed.
  if (!isTRUE(drawlabels) || method != "flattest" ||
      length(label.shift) == 0L || all(label.shift == 0)) {
    args <- list(
      x = x, y = y, z = z,
      nlevels = nlevels, levels = levels, labels = labels,
      xlim = xlim, ylim = ylim, zlim = zlim,
      labcex = labcex, drawlabels = drawlabels, method = method,
      axes = axes, frame.plot = frame.plot,
      col = col, lty = lty, lwd = lwd, add = add
    )
    if (have.vfont)
      args$vfont <- vfont
    args <- c(args, list(...))
    return(do.call(graphics::contour.default, args))
  }
  
  if (!is.numeric(label.shift) || any(!is.finite(label.shift)))
    stop("'label.shift' must contain finite numeric values")
  
  if (!is.null(labels)) {
    labels <- as.character(labels)
    if (!length(labels))
      stop("'labels' is length zero. Use 'drawlabels = FALSE' to suppress labels.")
  }
  
  # Reproduce contour.default() plot setup.
  if (!add) {
    localPlotWindow <- function(xlim, ylim, ..., main, sub, xlab, ylab,
                                outer, line)
      plot.window(xlim, ylim, ...)
    localTitle <- function(..., log) title(...)
    plot.new()
    localPlotWindow(xlim, ylim, ...)
    localTitle(...)
  }
  
  usr <- par("usr")
  rect <- c(min(usr[1:2]), max(usr[1:2]), min(usr[3:4]), max(usr[3:4]))
  
  cls <- grDevices::contourLines(x, y, z, levels = levels)
  
  # Base contour formats each level separately; format() is the closest
  # public-R equivalent of the internal labelformat() used in C.
  if (is.null(labels))
    label.text <- vapply(levels, function(v) format(v, trim = TRUE), character(1))
  else
    label.text <- labels
  
  # Recycling follows graphics::contour().
  if (!length(col)) col <- par("fg")
  if (!length(lty)) lty <- par("lty")
  if (!length(lwd)) lwd <- par("lwd")
  label.shift <- rep(label.shift, length.out = length(levels))
  
  label.boxes <- list()
  
  for (cl in cls) {
    k <- which.min(abs(levels - cl$level))
    if (!length(k) || !is.finite(k))
      next
    
    this.col <- col[(k - 1L) %% length(col) + 1L]
    if (is.na(this.col))
      this.col <- par("fg")
    this.lty <- lty[(k - 1L) %% length(lty) + 1L]
    if (is.na(this.lty))
      this.lty <- par("lty")
    this.lwd <- lwd[(k - 1L) %% length(lwd) + 1L]
    if (!is.finite(this.lwd))
      this.lwd <- par("lwd")
    
    lab <- label.text[(k - 1L) %% length(label.text) + 1L]
    buffer <- paste0(" ", lab, " ")
    
    sw.args <- list(s = buffer, units = "inches", cex = labcex)
    sh.args <- list(s = buffer, units = "inches", cex = labcex)
    if (have.vfont) {
      sw.args$vfont <- vfont
      sh.args$vfont <- vfont
    }
    label.width <- do.call(strwidth, sw.args)
    label.height <- do.call(strheight, sh.args)
    
    p <- cbind(cl$x, cl$y)
    runs <- .contour2_clip_polyline(p, rect)
    if (!length(runs))
      next
    
    # Find the flattest valid candidate over all visible connected runs.
    best <- NULL
    best.run <- NA_integer_
    run.candidates <- vector("list", length(runs))
    
    for (rr in seq_along(runs)) {
      cand <- .contour2_candidates(
        runs[[rr]], label.width, label.height,
        label.boxes, rect
      )
      run.candidates[[rr]] <- cand
      if (!length(cand))
        next
      
      vv <- vapply(cand, function(z) z$variance, numeric(1))
      ii <- which.min(vv)
      if (is.null(best) || vv[ii] < best$variance) {
        best <- cand[[ii]]
        best.run <- rr
      }
    }
    
    # No valid label location: draw all visible contour pieces normally.
    if (is.null(best)) {
      for (rr in seq_along(runs))
        .contour2_draw_run(runs[[rr]], col = this.col,
                           lty = this.lty, lwd = this.lwd)
      next
    }
    
    # Shift from the original flattest label position by a fraction of the
    # physical arc length of the same visible connected contour run.
    cand <- run.candidates[[best.run]]
    ai <- attr(cand, "arc")
    L <- ai$length
    target <- best$centre.s + label.shift[k] * L
    
    # Closed contours wrap naturally; open contours clamp to the visible run.
    closed <- nrow(runs[[best.run]]) > 2L &&
      sqrt(sum((runs[[best.run]][1L, ] -
                  runs[[best.run]][nrow(runs[[best.run]]), ])^2)) <=
      1e-10 * max(abs(rect), 1)
    if (closed && L > 0) {
      target <- target %% L
      dist.to.target <- vapply(cand, function(z) {
        d <- abs(z$centre.s - target)
        min(d, L - d)
      }, numeric(1))
    } else {
      target <- max(0, min(L, target))
      dist.to.target <- vapply(cand,
                               function(z) abs(z$centre.s - target),
                               numeric(1))
    }
    
    shifted <- cand[[which.min(dist.to.target)]]
    
    # Draw all runs, with one genuine label gap on the selected run.
    for (rr in seq_along(runs)) {
      if (rr == best.run) {
        .contour2_draw_run(
          runs[[rr]], shifted,
          label = buffer,
          label.width = label.width,
          labcex = labcex,
          col = this.col, lty = this.lty, lwd = this.lwd,
          vfont = if (have.vfont) vfont else NULL
        )
        label.boxes[[length(label.boxes) + 1L]] <- shifted$corners
      } else {
        .contour2_draw_run(runs[[rr]], col = this.col,
                           lty = this.lty, lwd = this.lwd)
      }
    }
  }
  
  if (!add) {
    localAxis <- function(..., col, bg, pch, cex, lty, lwd) Axis(...)
    localBox <- function(..., col, bg, pch, cex, lty, lwd) box(...)
    if (axes) {
      localAxis(x, side = 1, ...)
      localAxis(y, side = 2, ...)
    }
    if (frame.plot)
      localBox(...)
  }
  
  invisible()
}
