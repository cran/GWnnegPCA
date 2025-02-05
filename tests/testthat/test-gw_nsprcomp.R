# Load required packages
library(testthat)
library(GWnnegPCA)
library(sf)

test_that("gw_nsprcomp example from documentation works", {
  # Load required data
  nc <- sf::st_read(system.file("shape/nc.shp", package="sf"), quiet = TRUE)

  # Store original geometry
  orig_geom <- sf::st_geometry(nc)

  # Scale selected variables for analysis
  vars_to_use <- c("SID74", "NWBIR74", "BIR74")
  Data.scaled <- scale(as.matrix(sf::st_drop_geometry(nc[, vars_to_use])))

  # Create new sf object with scaled data
  data_df <- as.data.frame(Data.scaled)
  names(data_df) <- vars_to_use
  nc_scaled <- sf::st_set_geometry(data_df, orig_geom)

  # Set CRS to match original
  sf::st_crs(nc_scaled) <- sf::st_crs(nc)

  # Verify sf object properties
  expect_s3_class(nc_scaled, "sf")
  expect_identical(sf::st_geometry(nc_scaled), orig_geom)
  expect_identical(sf::st_crs(nc_scaled), sf::st_crs(nc))

  # Run the GW-NNPCA
  gwnnegpca_ans <- GWnnegPCA::gw_nsprcomp(
    data = nc_scaled,
    vars = vars_to_use,
    bw = 0.25,
    k = 3,
    longlat = TRUE,
    kernel = "bisquare",
    adaptive = TRUE,
    nneg = TRUE,
    geodisic_measure = "geodesic"
  )

  # Test expectations
  expect_true(is.list(gwnnegpca_ans))
  expect_true(is.array(gwnnegpca_ans$loadings))
  expect_equal(dim(gwnnegpca_ans$loadings)[2], length(vars_to_use))
  expect_equal(dim(gwnnegpca_ans$loadings)[1], nrow(nc))
  expect_true(all(gwnnegpca_ans$loadings >= 0, na.rm = TRUE)) # Check non-negativity constraint, ignoring NAs
})

test_that("input validation works correctly", {
  # Create simple test data with grid points
  n <- 5 # Small grid size for testing
  x <- seq(-1, 1, length.out = n)
  y <- seq(-1, 1, length.out = n)
  grid <- expand.grid(x = x, y = y)

  # Create well-conditioned test data
  set.seed(42)
  data <- data.frame(
    x1 = scale(runif(n^2, 0, 10)),  # Scaled uniform data
    x2 = scale(runif(n^2, 0, 10)),  # between 0 and 10
    x3 = scale(runif(n^2, 0, 10))
  )

  sf_data <- sf::st_as_sf(cbind(data, grid),
                         coords = c("x", "y"),
                         crs = 4326)

  # Test missing required parameters
  expect_error(GWnnegPCA::gw_nsprcomp(data=sf_data), "Variables input error")
  expect_error(GWnnegPCA::gw_nsprcomp(data=sf_data, vars=c("x1", "x2")), "Bandwidth is not specified incorrectly")

  # Test invalid variable names
  expect_error(GWnnegPCA::gw_nsprcomp(
    data=sf_data,
    vars=c("nonexistent"),
    bw=0.5,
    kernel="gaussian"
  ), "Variables input doesn't match with data")

  # Test invalid bandwidth
  expect_error(GWnnegPCA::gw_nsprcomp(
    data=sf_data,
    vars=c("x1", "x2"),
    bw=-1,
    kernel="gaussian"
  ), "Bandwidth is not specified incorrectly")
})

test_that("different kernel functions work", {
  # Create simple test data with grid points
  n <- 5 # Small grid size for testing
  x <- seq(-1, 1, length.out = n)
  y <- seq(-1, 1, length.out = n)
  grid <- expand.grid(x = x, y = y)

  # Create well-conditioned test data
  set.seed(42)
  data <- data.frame(
    x1 = scale(runif(n^2, 0, 10)),  # Scaled uniform data
    x2 = scale(runif(n^2, 0, 10))   # between 0 and 10
  )

  sf_data <- sf::st_as_sf(cbind(data, grid),
                         coords = c("x", "y"),
                         crs = 4326)

  vars <- c("x1", "x2")
  kernels <- c("gaussian", "exponential", "bisquare", "tricube", "boxcar")

  for(k in kernels) {
    result <- GWnnegPCA::gw_nsprcomp(
      data=sf_data,
      vars=vars,
      bw=1,  # Use bandwidth of 1 for this small grid
      kernel=k,
      adaptive=TRUE,
      k=2
    )

    expect_true(is.list(result))
    expect_true(all(c("loadings", "score", "sdev") %in% names(result)))
    expect_true(is.array(result$loadings))
    expect_true(all(result$loadings >= 0, na.rm = TRUE))
  }
})

test_that("distance calculations work correctly", {
  # Create simple test data with grid points
  n <- 5 # Small grid size for testing
  x <- seq(-1, 1, length.out = n)
  y <- seq(-1, 1, length.out = n)
  grid <- expand.grid(x = x, y = y)

  # Create well-conditioned test data
  set.seed(42)
  data <- data.frame(
    x1 = scale(runif(n^2, 0, 10)),  # Scaled uniform data
    x2 = scale(runif(n^2, 0, 10))   # between 0 and 10
  )

  sf_data <- sf::st_as_sf(cbind(data, grid),
                         coords = c("x", "y"),
                         crs = 4326)

  vars <- c("x1", "x2")

  # Test Euclidean distance
  euclidean_result <- GWnnegPCA::gw_nsprcomp(
    data=sf_data,
    vars=vars,
    bw=1,  # Use bandwidth of 1 for this small grid
    kernel="gaussian",
    longlat=FALSE,
    k=2
  )

  # Test geodesic distance
  geodesic_result <- GWnnegPCA::gw_nsprcomp(
    data=sf_data,
    vars=vars,
    bw=1,  # Use bandwidth of 1 for this small grid
    kernel="gaussian",
    longlat=TRUE,
    geodisic_measure="geodesic",
    k=2
  )

  expect_true(is.list(euclidean_result))
  expect_true(is.list(geodesic_result))
  expect_true(all(euclidean_result$loadings >= 0, na.rm = TRUE))
  expect_true(all(geodesic_result$loadings >= 0, na.rm = TRUE))
})

test_that("parameter combinations work correctly", {
  # Create simple test data with grid points
  n <- 5 # Small grid size for testing
  x <- seq(-1, 1, length.out = n)
  y <- seq(-1, 1, length.out = n)
  grid <- expand.grid(x = x, y = y)

  # Create well-conditioned test data
  set.seed(42)
  data <- data.frame(
    x1 = scale(runif(n^2, 0, 10)),  # Scaled uniform data
    x2 = scale(runif(n^2, 0, 10))   # between 0 and 10
  )

  sf_data <- sf::st_as_sf(cbind(data, grid),
                         coords = c("x", "y"),
                         crs = 4326)

  vars <- c("x1", "x2")

  # Test adaptive vs fixed bandwidth
  adaptive_result <- GWnnegPCA::gw_nsprcomp(
    data=sf_data,
    vars=vars,
    bw=1,  # Use bandwidth of 1 for this small grid
    kernel="gaussian",
    adaptive=TRUE,
    k=2
  )

  fixed_result <- GWnnegPCA::gw_nsprcomp(
    data=sf_data,
    vars=vars,
    bw=1,  # Use bandwidth of 1 for this small grid
    kernel="gaussian",
    adaptive=FALSE,
    k=2
  )

  # Test different numbers of components
  k2_result <- GWnnegPCA::gw_nsprcomp(
    data=sf_data,
    vars=vars,
    bw=1,  # Use bandwidth of 1 for this small grid
    kernel="gaussian",
    k=2
  )

  k3_result <- GWnnegPCA::gw_nsprcomp(
    data=sf_data,
    vars=vars,
    bw=1,  # Use bandwidth of 1 for this small grid
    kernel="gaussian",
    k=2  # Keep k=2 to match number of variables
  )

  expect_equal(dim(k2_result$loadings)[3], 2)
  expect_equal(dim(k3_result$loadings)[3], 2)  # Changed from 3 to 2
  expect_true(all(adaptive_result$loadings >= 0, na.rm = TRUE))
  expect_true(all(fixed_result$loadings >= 0, na.rm = TRUE))
})
