test_that("smoothclust validates parameter values", {
  input <- matrix(1:12, nrow = 3)
  spatial_coords <- cbind(0:3, 0)
  
  expect_error(smoothclust(input, spatial_coords = spatial_coords, bandwidth = 0))
  expect_error(smoothclust(input, spatial_coords = spatial_coords, bandwidth = c(0.1, 0.2)))
  expect_error(smoothclust(input, spatial_coords = spatial_coords, k = 1.5))
  expect_error(smoothclust(input, spatial_coords = spatial_coords, truncate = 1))
  expect_error(smoothclust(input, spatial_coords = spatial_coords, truncate = NA_real_))
  expect_error(smoothclust(input, spatial_coords = spatial_coords, n_threads = 0))
  expect_error(smoothclust(input, spatial_coords = spatial_coords, n_threads = 1.5))
})

test_that("smoothclust validates preprocessing flags", {
  input <- matrix(1:12, nrow = 3)
  spatial_coords <- cbind(0:3, 0)

  expect_error(smoothclust(input, spatial_coords = spatial_coords, normalize = 1))
  expect_error(smoothclust(input, spatial_coords = spatial_coords, normalize = NA))
  expect_error(smoothclust(input, spatial_coords = spatial_coords,
                           normalize = c(TRUE, FALSE)))
  expect_error(smoothclust(input, spatial_coords = spatial_coords, log = 1))
  expect_error(smoothclust(input, spatial_coords = spatial_coords, log = NA))
  expect_error(smoothclust(input, spatial_coords = spatial_coords,
                           log = c(TRUE, FALSE)))
})

test_that("smoothclust validates size factors", {
  input <- matrix(1:12, nrow = 3,
                  dimnames = list(NULL, paste0("spot", 1:4)))
  spatial_coords <- cbind(0:3, 0)

  expect_error(smoothclust(input, spatial_coords = spatial_coords,
                           size_factors = 1:3))
  expect_error(smoothclust(input, spatial_coords = spatial_coords,
                           size_factors = c(1, 1, 0, 1)))
  expect_error(smoothclust(input, spatial_coords = spatial_coords,
                           size_factors = c(1, 1, -1, 1)))
  expect_error(smoothclust(input, spatial_coords = spatial_coords,
                           size_factors = c(1, 1, NA_real_, 1)))
  expect_error(smoothclust(input, spatial_coords = spatial_coords,
                           size_factors = c(1, 1, Inf, 1)))
  expect_error(smoothclust(input, spatial_coords = spatial_coords,
                           size_factors = letters[1:4]))
  expect_error(smoothclust(input, spatial_coords = spatial_coords,
                           size_factors = matrix(1, nrow = 4)))
  expect_error(smoothclust(input, spatial_coords = spatial_coords,
                           normalize = FALSE, size_factors = rep(1, 4)))

  duplicated_names <- setNames(rep(1, 4), c("spot1", "spot1", "spot3", "spot4"))
  missing_name <- setNames(rep(1, 4), c("spot1", "spot2", "spot3", "other"))
  expect_error(smoothclust(input, spatial_coords = spatial_coords,
                           size_factors = duplicated_names))
  expect_error(smoothclust(input, spatial_coords = spatial_coords,
                           size_factors = missing_name))

  unlabelled_input <- unname(input)
  named_factors <- setNames(rep(1, 4), paste0("spot", 1:4))
  expect_error(smoothclust(unlabelled_input, spatial_coords = spatial_coords,
                           size_factors = named_factors))
})

test_that("smoothclust rejects invalid values for requested preprocessing", {
  spatial_coords <- cbind(0:3, 0)
  zero_library <- matrix(c(1, 2, 3, 0, 0, 0, 4, 5, 6, 7, 8, 9), nrow = 3)
  negative <- matrix(1:12, nrow = 3)
  negative[1, 1] <- -1
  non_finite <- matrix(as.numeric(1:12), nrow = 3)
  non_finite[1, 1] <- Inf

  expect_error(
    smoothclust(zero_library, spatial_coords = spatial_coords),
    "positive sums"
  )
  expect_no_error(smoothclust(
    zero_library,
    spatial_coords = spatial_coords,
    size_factors = rep(1, ncol(zero_library))
  ))
  expect_error(smoothclust(negative, spatial_coords = spatial_coords))
  expect_error(smoothclust(negative, spatial_coords = spatial_coords,
                           normalize = FALSE, log = TRUE))
  expect_error(smoothclust(non_finite, spatial_coords = spatial_coords))

  # Pure smoothing continues to support arbitrary finite numeric values.
  expect_no_error(smoothclust(negative, spatial_coords = spatial_coords,
                              normalize = FALSE, log = FALSE))
})

test_that("smoothclust validates log pseudocounts and output assay names", {
  input <- matrix(1:12, nrow = 3)
  spatial_coords <- cbind(0:3, 0)

  for (pseudo_count in list(0, -1, NA_real_, Inf, c(1, 2), "1")) {
    expect_error(smoothclust(input, spatial_coords = spatial_coords,
                             pseudo_count = pseudo_count))
  }

  for (output_assay_name in list(NA_character_, "", c("a", "b"), 1)) {
    expect_error(smoothclust(input, spatial_coords = spatial_coords,
                             output_assay_name = output_assay_name))
  }
})

test_that("smoothclust safeguards assays that are already transformed", {
  input <- matrix(1:8, nrow = 2,
                  dimnames = list(NULL, paste0("spot", 1:4)))
  spatial_coords <- cbind(x = 0:3, y = 0)
  rownames(spatial_coords) <- colnames(input)
  spe <- SpatialExperiment::SpatialExperiment(
    assays = list(logcounts = log2(input + 1), normcounts = input),
    spatialCoords = spatial_coords
  )

  expect_error(smoothclust(spe, assay_name = "logcounts", bandwidth = 0.4))
  expect_no_error(smoothclust(spe, assay_name = "logcounts", bandwidth = 0.4,
                              normalize = FALSE, log = FALSE))
  expect_error(smoothclust(spe, assay_name = "normcounts", bandwidth = 0.4,
                           normalize = TRUE, log = FALSE))
  expect_no_error(smoothclust(spe, assay_name = "normcounts", bandwidth = 0.4,
                              normalize = FALSE, log = TRUE))
})

test_that("smoothclust validates spatial coordinate compatibility", {
  input <- matrix(1:12, nrow = 3)
  spatial_coords <- cbind(0:3, 0)
  
  expect_error(smoothclust(input, spatial_coords = spatial_coords[-1, ]))
  
  spatial_coords_bad <- spatial_coords
  spatial_coords_bad[1, 1] <- Inf
  expect_error(smoothclust(input, spatial_coords = spatial_coords_bad))
})

test_that("smoothclust validates k when using knn smoothing", {
  input <- matrix(1:12, nrow = 3)
  spatial_coords <- cbind(0:3, 0)
  
  expect_error(smoothclust(input, spatial_coords = spatial_coords, method = "knn", k = 4))
})

test_that("boundary_density validates inputs", {
  spatial_coords <- cbind(0:3, 0)
  labels <- c(1, 1, 2, 2)
  
  expect_error(boundary_density(spatial_coords[-1, ], labels, k = 1))
  
  spatial_coords_bad <- spatial_coords
  spatial_coords_bad[1, 1] <- NA_real_
  expect_error(boundary_density(spatial_coords_bad, labels, k = 1))
  
  expect_error(boundary_density(spatial_coords, labels, k = 0))
  expect_error(boundary_density(spatial_coords, labels, k = 1.5))
  expect_error(boundary_density(spatial_coords, labels, k = 4))
  expect_error(boundary_density(spatial_coords, labels, n_threads = 0))
  expect_error(boundary_density(spatial_coords, labels, n_threads = Inf))
  expect_error(boundary_density(spatial_coords, labels, adjust = "bad"))
  expect_error(boundary_density(spatial_coords, labels, 
                                adjust = "permutation", n_permutations = 1))
  expect_error(boundary_density(spatial_coords, labels, 
                                adjust = "permutation", n_permutations = 1.5))
  expect_error(boundary_density(spatial_coords, labels, 
                                adjust = "permutation", seed = NA_real_))
  expect_error(boundary_density(spatial_coords, c(1, 1, 2, NA), k = 1))
  expect_error(boundary_density(spatial_coords, list(1, 1, 2, 2), k = 1))
  expect_error(boundary_density(spatial_coords, 
                                matrix(c(1, 1, 2, 2), ncol = 1), k = 1))
})
