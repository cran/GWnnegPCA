#' Geographically Weighted Non-negative Principal Component Analysis
#'
#' @description
#' Implementation of geographically weighted non-negative principal component analysis,
#' which consists of the fusion of GWPCA and sparse non-negative PCA.
#'
#' @param data An sf object containing the spatial data and attributes for analysis
#' @param elocat Two-column numeric array or sf object for providing evaluation locations
#' @param vars Character vector of variable names to be used in the analysis
#' @param bw Bandwidth used in the weighting function
#' @param k The number of retained components (default: 2)
#' @param kernel Kernel function type: "gaussian", "exponential", "bisquare", "tricube", or "boxcar"
#' @param adaptive If TRUE, calculate adaptive kernel (default: TRUE)
#' @param p Power of the Minkowski distance (default: 2)
#' @param theta Angle in radians to rotate coordinate system (default: 0)
#' @param longlat If TRUE, great circle distances will be calculated (default: FALSE)
#' @param geodisic_measure Method for geodesic distance calculation (default: "cheap")
#' @param dMat Pre-specified distance matrix (default: NULL)
#' @param n.obs Number of observations for correlation matrix (default: NA)
#' @param n.iter Number of bootstrap iterations (default: 1)
#' @param ncomp Number of principal components to compute (default: k)
#' @param nneg If TRUE, constrain loadings to be non-negative (default: TRUE)
#' @param localcenter If TRUE, center local weighted x (default: TRUE)
#' @param localscale If TRUE, scale local weighted x (default: FALSE)
#' @param ... Additional arguments passed to methods
#'
#' @return A list containing:
#' \item{loadings}{The localized loadings matrix}
#' \item{score}{The PC score matrix from the localized non-negative PCA}
#' \item{sdev}{The localized standard deviation vector of the principal components}
#'
#' @importFrom methods is
#' @importFrom sf st_centroid st_crs st_coordinates st_geometry st_drop_geometry st_read
#' @importFrom nsprcomp nsprcomp
#' @importFrom geodist geodist
#'
#' @export
#'
#' @examples
#' # Read North Carolina SIDS data from sf package
#' nc <- sf::st_read(system.file("shape/nc.shp", package="sf"), quiet = TRUE)
#'
#' # Scale selected variables for analysis
#' vars_to_use <- c("SID74", "NWBIR74", "BIR74")
#' Data.scaled <- scale(as.matrix(sf::st_drop_geometry(nc[, vars_to_use])))
#'
#' # Create sf object with scaled data
#' nc_scaled <- nc
#' nc_scaled[vars_to_use] <- Data.scaled
#'
#' gwnnegpca_ans <- gw_nsprcomp(
#'   data = nc_scaled,
#'   vars = vars_to_use,
#'   bw = 0.25,
#'   k = 3,
#'   longlat = TRUE,
#'   kernel = "bisquare",
#'   adaptive = TRUE,
#'   nneg = TRUE,
#'   geodisic_measure = "geodesic"
#' )
gw_nsprcomp <-
  function (data,
            elocat,
            vars,
            bw,
            k = 2,
            kernel = "gaussian",
            adaptive = TRUE,
            p = 2,
            theta = 0,
            longlat = FALSE,
            geodisic_measure = "cheap",
            dMat = NULL,
            n.obs = NA,
            n.iter = 1,
            ncomp = k,
            nneg = TRUE,
            localcenter = TRUE,
            localscale = FALSE,
            ...)
  {
    # Validate kernel type
    valid_kernels <- c("gaussian", "exponential", "bisquare", "tricube", "boxcar")
    if (!kernel %in% valid_kernels) {
      stop(sprintf("Invalid kernel type. Must be one of: %s",
                  paste(valid_kernels, collapse = ", ")))
    }

    # Internal function to calculate Euclidean distance matrix
    distmat <- function(dp.locat, elocat) {
      if (!is.matrix(dp.locat)) dp.locat <- as.matrix(dp.locat)
      if (!is.matrix(elocat)) elocat <- as.matrix(elocat)

      n <- nrow(dp.locat)
      m <- nrow(elocat)
      dMat <- matrix(0, n, m)

      for (i in 1:m) {
        dMat[, i] <- sqrt((dp.locat[,1] - elocat[i,1])^2 +
                         (dp.locat[,2] - elocat[i,2])^2)
      }
      return(dMat)
    }

    # Internal function to calculate kernel weights
    weight_func <- function(type, adapt, dist_vec, bw) {
      dist_vec <- as.double(dist_vec)

      if (adapt) {
        bw_size <- as.integer(length(dist_vec) * bw)
        bw_dist <- as.double(sort(dist_vec)[bw_size])

        if (type == "gaussian") {
          weight <- exp((-0.5) * ((dist_vec ^ 2) / (bw_dist ^ 2)))
        } else if (type == "exponential") {
          weight <- exp((-1) * dist_vec / bw_dist)
        } else if (type == "bisquare") {
          weight <- ifelse((dist_vec > bw_dist), 0, (1 - (dist_vec ^ 2 / bw_dist ^ 2)) ^ 2)
        } else if (type == "tricube") {
          weight <- ifelse((dist_vec > bw_dist), 0, (1 - (dist_vec ^ 3 / bw_dist ^ 3)) ^ 3)
        } else if (type == "boxcar") {
          weight <- ifelse((dist_vec > bw_dist), 0, 1)
        }
      } else {
        if (type == "gaussian") {
          weight <- exp((-0.5) * ((dist_vec ^ 2) / (bw ^ 2)))
        } else if (type == "exponential") {
          weight <- exp((-1) * dist_vec / bw)
        } else if (type == "bisquare") {
          weight <- ifelse((dist_vec > bw), 0, (1 - (dist_vec ^ 2 / bw ^ 2)) ^ 2)
        } else if (type == "tricube") {
          weight <- ifelse((dist_vec > bw), 0, (1 - (dist_vec ^ 3 / bw ^ 3)) ^ 3)
        } else if (type == "boxcar") {
          weight <- ifelse((dist_vec > bw), 0, 1)
        }
      }
      return(weight)
    }

    # Internal function to perform weighted non-negative PCA
    w_nsprcomp <- function(x,
                          wt,
                          ncomp,
                          nneg = nneg,
                          localcenter = localcenter,
                          localscale = localscale,
                          ...) {
      wt_x <- x * wt
      nsprcomp::nsprcomp(
        wt_x,
        ncomp,
        nneg = nneg,
        localcenter = localcenter,
        localscale = localscale,
        ...
      )
    }

    if (methods::is(data, "sf")) {
      p4s <- sf::st_crs(data)
      dp.locat <- sf::st_coordinates(sf::st_centroid(data))
      names(dp.locat) <- c("longitude", "latitude")
    }
    else if (methods::is(data, "data.frame") && (!missing(dMat)))
      data <- data
    else
      stop("Given data must be an sf object or data.frame object")

    if (missing(elocat)) {
      ep.given <- FALSE
      elocat <- sf::st_coordinates(sf::st_centroid(data))
      names(elocat) <- c("longitude", "latitude")
    }
    else {
      ep.given <- TRUE
      if (methods::is(elocat, "sf")) {
        espdf <- elocat
        elocat <- sf::st_coordinates(sf::st_centroid(espdf))
        names(elocat) <- c("longitude", "latitude")
      }
      else if (is.numeric(elocat) && dim(elocat)[2] == 2)
        elocat <- elocat
      else {
        warning("Output locations are not packed in a Spatial object, and it has to be a two-column numeric vector")
        elocat <- dp.locat
      }
    }

    if (methods::is(data, "sf")) {
      data <- sf::st_drop_geometry(data)
    }

    dp.n <- nrow(data)
    ep.n <- nrow(elocat)

    if (missing(dMat)) {
      DM.given <- FALSE
      DM1.given <- FALSE

      if (longlat) {
        dMat <- geodist::geodist(dp.locat,
                                elocat,
                                measure = geodisic_measure)
      } else {
        dMat <- distmat(dp.locat, elocat)
      }

      DM.given <- TRUE
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
    sdev <- matrix(NA, ep.n, k)
    score.all <- matrix(NA, ep.n, k)

    for (i in 1:ep.n) {
      dist.vi <- dMat[, i]

      wt <- weight_func(
        type = kernel,
        adapt = adaptive,
        dist_vec = dist.vi,
        bw = bw
      )

      use <- wt > 0
      wt <- wt[use]
      if (length(wt) <= 5) {
        expr <- paste("Too small bandwidth at location: ", i)
        warning(paste(expr, "and the results can't be given there.",
                     sep = ", "))
        next
      }

      temp <- w_nsprcomp(
        x = x[use,],
        wt = wt,
        ncomp = k,
        k = ncol(x),
        nneg = nneg,
        localcenter = localcenter,
        localscale = localscale
      )

      load[i, ,] <- matrix(temp$rotation, ncol = k, nrow = var.n)
      score.all[use,] <- (as.matrix(x[use, ]) %*% as.matrix(temp$rotation))[, 1:k]
      s[i, ] <- score.all[i, ]
      sdev[i,] <- temp$sdev
    }

    if (!is.null(rownames(x)))
      dimnames(load)[[1]] <- rownames(x)
    if (!is.null(colnames(x)))
      dimnames(load)[[2]] <- colnames(x)
    dimnames(load)[[3]] <- paste("PC", 1:k, sep = "")

    list(loadings = load,
         score = s,
         sdev = sdev)
  }
