gw_nsprcomp <- function (data, elocat, vars, bw, k = 2, kernel, adaptive = TRUE,
                         p = 2, theta = 0, longlat = FALSE, dMat = NULL, n.obs = NA,
                         n.iter = 1, ncomp = k, nneg = TRUE, localcenter = TRUE, localscale = FALSE,...)
{
  requireNamespace("GWmodel")
  requireNamespace("nsprcomp")
  w_nsprcomp <- function(x, wt, ncomp, nneg = nneg, localcenter = localcenter, localscale = localscale, ...) {
    wt_x <- x * wt
    nsprcomp(wt_x, ncomp, nneg = nneg, localcenter = localcenter,  localscale = localscale, ...)
    }

  if (is(data, "Spatial")) {
    p4s <- proj4string(data)
    dp.locat <- coordinates(data)
  }
  else if (is(data, "data.frame") && (!missing(dMat)))
    data <- data
  else stop("Given data must be a Spatial*DataFrame or data.frame object")
  if (missing(elocat)) {
    ep.given <- FALSE
    elocat <- coordinates(data)
  }
  else {
    ep.given <- TRUE
    if (is(elocat, "Spatial")) {
      espdf <- elocat
      elocat <- coordinates(espdf)
    }
    else if (is.numeric(elocat) && dim(elocat)[2] == 2)
      elocat <- elocat
    else {
      warning("Output loactions are not packed in a Spatial object,and it has to be a two-column numeric vector")
      elocat <- dp.locat
    }
  }
  data <- as(data, "data.frame")
  dp.n <- nrow(data)
  ep.n <- nrow(elocat)
  if (missing(dMat)) {
    DM.given <- FALSE
    DM1.given <- FALSE
    if (dp.n + ep.n <= 10000) {
      dMat <- gw.dist(dp.locat = dp.locat, rp.locat = elocat,
                      p = p, theta = theta, longlat = longlat)
      DM.given <- TRUE
    }
  }
  else {
    DM.given <- TRUE
    DM1.given <- TRUE
    dim.dMat <- dim(dMat)
    if (dim.dMat[1] != dp.n || dim.dMat[2] != ep.n)
      stop("Dimensions of dMat are not correct")
  }
  if (missing(vars))
    stop("Variables input error")
  if (missing(bw) || bw <= 0)
    stop("Bandwidth is not specified incorrectly")
  len.var <- length(vars)
  col.nm <- colnames(data)
  var.idx <- match(vars, col.nm)[!is.na(match(vars, col.nm))]
  if (length(var.idx) == 0)
    stop("Variables input doesn't match with data")
  x <- data[, var.idx]
  x <- as.matrix(x)
  var.nms <- colnames(x)
  var.n <- ncol(x)
  if (len.var > var.n)
    warning("Invalid variables have been specified, please check them again!")
  load <- array(NA, c(ep.n, var.n, k))
  resid_sqsum <- rep(NA, nrow(elocat))
  s <- matrix(NA, ep.n, k)
  sdev <- matrix(NA, ep.n, var.n)
  score.all <- matrix(NA, ep.n, k)
  for (i in 1:ep.n) {
    if (DM.given)
      dist.vi <- dMat[, i]
    else {
      if (ep.given)
        dist.vi <- gw.dist(dp.locat, elocat, focus = i,
                           p, theta, longlat)
      else dist.vi <- gw.dist(dp.locat, focus = i, p = p,
                              theta = theta, longlat = longlat)
    }
    wt <- gw.weight(dist.vi, bw, kernel, adaptive)
    use <- wt > 0
    wt <- wt[use]
    if (length(wt) <= 5) {
      expr <- paste("Too small bandwidth at location: ",
                    i)
      warning(paste(expr, "and the results can't be given there.",
                    sep = ", "))
      next
    }
    temp <- w_nsprcomp(x = x[use, ], wt, ncomp = k, k = ncol(x),
                       nneg, localcenter, localscale)
    load[i, , ] <- matrix(temp$rotation, ncol = k, nrow = var.n)
    score.all[use, ] <- temp$x[, 1:k]
    s[i, ] <- score.all[i, ]
    sdev[i, ] <- temp$sdev
  }
  if (!is.null(rownames(x)))
    dimnames(load)[[1]] <- rownames(x)
  if (!is.null(colnames(x)))
    dimnames(load)[[2]] <- colnames(x)
  dimnames(load)[[3]] <- paste("PC", 1:k, sep = "")
  list(loadings = load, score = s, sdev = sdev)
}
