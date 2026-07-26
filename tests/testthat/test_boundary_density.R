test_that("results object has correct structure", {
  spatial_coords <- cbind(c(0, 1, 3, 6), 0)
  labels <- c("A", "B", "A", "A")
  
  res <- boundary_density(spatial_coords, labels, k = 2)
  
  expect_is(res, "list")
  expect_named(res, c("n_discordant", "n_neighbors", "local_boundary_density", 
                      "mean_discordant", "boundary_density"))
})

test_that("boundary_density returns analytic composition adjustment", {
  spatial_coords <- cbind(c(0, 1, 3, 6), 0)
  labels <- c(1, 2, 1, 1)
  
  res <- boundary_density(spatial_coords, labels, k = 2, adjust = "analytic")
  
  expect_named(res, c("n_discordant", "n_neighbors",
                      "local_boundary_density", "mean_discordant",
                      "boundary_density", "expected_boundary_density",
                      "relative_boundary_density"))
  expect_equal(res$boundary_density, 0.625)
  expect_equal(res$expected_boundary_density, 0.5)
  expect_equal(res$relative_boundary_density, 1.25)
})

test_that("boundary_density supports character and factor labels", {
  spatial_coords <- cbind(c(0, 1, 3, 6), 0)
  labels_chr <- c("A", "B", "A", "A")
  labels_fac <- factor(labels_chr)
  
  res_chr <- boundary_density(spatial_coords, labels_chr, k = 2, 
                              adjust = "analytic")
  res_fac <- boundary_density(spatial_coords, labels_fac, k = 2, 
                              adjust = "analytic")
  
  expect_equal(res_chr$n_discordant, c(1, 2, 1, 1))
  expect_equal(res_fac$n_discordant, res_chr$n_discordant)
  expect_equal(res_fac$boundary_density, res_chr$boundary_density)
  expect_equal(res_fac$expected_boundary_density, 
               res_chr$expected_boundary_density)
})

test_that("boundary_density handles one-label composition references", {
  spatial_coords <- cbind(c(0, 1, 3, 6), 0)
  labels <- c(1, 1, 1, 1)
  
  analytic <- boundary_density(spatial_coords, labels, k = 2, 
                               adjust = "analytic")
  
  expect_equal(analytic$boundary_density, 0)
  expect_equal(analytic$expected_boundary_density, 0)
  expect_true(is.na(analytic$relative_boundary_density))
})
