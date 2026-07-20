set.seed(123)
input <- matrix(rpois(60, lambda = 10), nrow = 10)
spatial_coords <- cbind(rep(1:3, each = 2), rep(1:2, times = 3))
out <- smoothclust(input, spatial_coords = spatial_coords, bandwidth = 0.6)


test_that("smoothclust example output has correct class", {
  expect_s4_class(out, "dgCMatrix")
  expect_s4_class(out, "Matrix")
})

test_that("smoothclust example output has correct dimensions", {
  expect_equal(dim(out), c(10, 6))
})

test_that("smoothclust example output values are stable", {
  size_factors <- colSums(input) / mean(colSums(input))
  input_norm <- sweep(input, 2, size_factors, "/")
  out_expected <- smoothclust(
    input_norm,
    spatial_coords = spatial_coords,
    bandwidth = 0.6,
    normalize = FALSE,
    log = FALSE
  )
  out_expected <- log2(out_expected + 1)

  expect_equal(
    as.matrix(out),
    as.matrix(out_expected)
  )
})

test_that("smoothclust retains original smoothing behavior when requested", {
  out_legacy <- smoothclust(
    input,
    spatial_coords = spatial_coords,
    bandwidth = 0.6,
    normalize = FALSE,
    log = FALSE
  )

  expect_equal(
    unname(signif(as.numeric(out_legacy[1:6, 1]), 6)),
    c(8.33333, 9.33333, 9.33333, 8, 10.3333, 11.3333)
  )
})
