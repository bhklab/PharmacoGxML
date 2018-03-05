myScatterPlot <-
  function(x, y, method=c("plain", "transparent", "smooth"), transparency=0.10, smooth.pch=".", pch=16, minp=50, col=blues9[7], smooth.col=c("white", blues9), ...) {
    require(grDevices) || stop("Library grDevices is not available!")
    method <- match.arg(method)
    if (length(col) != length(x)) {
      col <- rep(col, length.out=length(x))
    }
    ccix <- complete.cases(x, y)
    x <- x[ccix]
    y <- y[ccix]
    col <- col[ccix]

    if (sum(ccix) < minp) {
      ## too few points, no transparency, no smoothing
      if (sum(ccix) > 0) { rr <- plot(x=x, y=y, col=col, pch=pch, ...) } else { rr <- plot(x=x, y=y, col=col, pch=pch, ...) }
    } else {
      ## enough data points
      switch(method,
             "plain"={
               rr <- plot(x=x, y=y, col=col, pch=pch, ...)
             },
             "transparent"={
               myrgb <- sapply(col, grDevices::col2rgb, alpha=FALSE) / 255
               myrgb <- apply(myrgb, 2, function (x, transparency) {
                 return (rgb(red=x[1], green=x[2], blue=x[3], alpha=transparency, maxColorValue=1))
               }, transparency=transparency)
               rr <- plot(x=x, y=y, pch=pch, col=myrgb, ...)
             },
             "smooth"={
               rr <- smoothScatter(x=x, y=y, col="lightgray", colramp=colorRampPalette(smooth.col), pch=smooth.pch, ...)
             }
      )
    }

    invisible(rr)
  }
