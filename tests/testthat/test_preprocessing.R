manual_uniform_smooth <- function(x) {
  first_pair <- rowMeans(x[, 1:2, drop = FALSE])
  out <- cbind(first_pair, first_pair,
               x[, 3, drop = FALSE], x[, 4, drop = FALSE])
  dimnames(out) <- dimnames(x)
  out
}

test_that("smoothclust implements every preprocessing flag combination", {
  input <- matrix(
    c(2, 6, 4, 2, 8, 4, 12, 8),
    nrow = 2,
    dimnames = list(paste0("gene", 1:2), paste0("spot", 1:4))
  )
  spatial_coords <- cbind(c(0, 1, 3, 6), 0)
  size_factors <- colSums(input) / mean(colSums(input))
  normalized <- sweep(input, 2, size_factors, "/")
  raw_smooth <- manual_uniform_smooth(input)
  normalized_smooth <- manual_uniform_smooth(normalized)

  out_raw <- smoothclust(
    input, spatial_coords = spatial_coords, bandwidth = 1.1 / 6,
    normalize = FALSE, log = FALSE
  )
  out_normalized <- smoothclust(
    input, spatial_coords = spatial_coords, bandwidth = 1.1 / 6,
    normalize = TRUE, log = FALSE
  )
  out_logged <- smoothclust(
    input, spatial_coords = spatial_coords, bandwidth = 1.1 / 6,
    normalize = FALSE, log = TRUE
  )
  out_default <- smoothclust(
    input, spatial_coords = spatial_coords, bandwidth = 1.1 / 6,
    normalize = TRUE, log = TRUE
  )

  expect_equal(as.matrix(out_raw), raw_smooth)
  expect_equal(as.matrix(out_normalized), normalized_smooth)
  expect_equal(as.matrix(out_logged), log2(raw_smooth + 1))
  expect_equal(as.matrix(out_default), log2(normalized_smooth + 1))
  expect_identical(dimnames(out_default), dimnames(input))
})

test_that("automatic library-size factors are centered to mean one", {
  input <- matrix(
    c(2, 6, 4, 2, 8, 4, 12, 8),
    nrow = 2,
    dimnames = list(NULL, paste0("spot", 1:4))
  )
  spatial_coords <- cbind(c(0, 1, 3, 6), 0)
  expected_factors <- colSums(input) / mean(colSums(input))

  out <- smoothclust(
    input,
    spatial_coords = spatial_coords,
    bandwidth = 0.01,
    normalize = TRUE,
    log = FALSE
  )

  expect_equal(as.matrix(out), sweep(input, 2, expected_factors, "/"))
})

test_that("supplied size factors are aligned by name and centered", {
  input <- matrix(
    c(2, 6, 4, 2, 8, 4, 12, 8),
    nrow = 2,
    dimnames = list(NULL, paste0("spot", 1:4))
  )
  spatial_coords <- cbind(c(0, 1, 3, 6), 0)
  supplied <- c(spot4 = 4, spot2 = 2, spot1 = 1, spot3 = 3)
  expected_factors <- c(spot1 = 1, spot2 = 2, spot3 = 3, spot4 = 4) / 2.5

  out <- smoothclust(
    input,
    spatial_coords = spatial_coords,
    bandwidth = 0.01,
    normalize = TRUE,
    log = FALSE,
    size_factors = supplied
  )

  expect_equal(as.matrix(out), sweep(input, 2, expected_factors, "/"))
})

test_that("supplied size factors permit zero-sum filtered columns", {
  input <- matrix(
    c(2, 6, 0, 0, 4, 8),
    nrow = 2,
    dimnames = list(NULL, paste0("spot", 1:3))
  )
  spatial_coords <- cbind(c(0, 5, 10), 0)
  supplied <- c(2, 4, 6)
  centered <- supplied / mean(supplied)
  expected <- sweep(input, 2, centered, "/")

  out <- smoothclust(
    input,
    spatial_coords = spatial_coords,
    bandwidth = 0.01,
    normalize = TRUE,
    log = FALSE,
    size_factors = supplied
  )
  out_sparse <- smoothclust(
    Matrix::Matrix(input, sparse = TRUE),
    spatial_coords = spatial_coords,
    bandwidth = 0.01,
    normalize = TRUE,
    log = FALSE,
    size_factors = supplied
  )
  rownames(spatial_coords) <- colnames(input)
  spe <- SpatialExperiment::SpatialExperiment(
    assays = list(counts = input), spatialCoords = spatial_coords
  )
  out_spe <- smoothclust(
    spe,
    bandwidth = 0.01,
    normalize = TRUE,
    log = FALSE,
    size_factors = supplied
  )

  expect_equal(as.matrix(out), expected)
  expect_equal(as.matrix(out_sparse), expected)
  expect_equal(
    as.matrix(SummarizedExperiment::assays(out_spe)$normcounts_smooth),
    expected
  )
})

test_that("normalization occurs before spatial smoothing", {
  input <- matrix(c(10, 40, 30), nrow = 1)
  spatial_coords <- cbind(c(0, 1, 10), 0)
  size_factors <- c(0.75, 1.5, 0.75)

  out <- smoothclust(
    input,
    spatial_coords = spatial_coords,
    bandwidth = 0.11,
    normalize = TRUE,
    log = FALSE,
    size_factors = size_factors
  )

  normalized <- sweep(input, 2, size_factors, "/")
  expected <- matrix(c(mean(normalized[1, 1:2]),
                       mean(normalized[1, 1:2]), normalized[1, 3]), nrow = 1)
  raw_smoothed <- matrix(c(25, 25, 30), nrow = 1)
  smooth_then_normalize <- sweep(raw_smoothed, 2, size_factors, "/")

  expect_equal(as.matrix(out), expected)
  expect_false(isTRUE(all.equal(as.matrix(out), smooth_then_normalize)))
})

test_that("log transformation occurs after smoothing and honors pseudocount", {
  input <- matrix(c(0, 8, 16), nrow = 1)
  spatial_coords <- cbind(c(0, 1, 10), 0)

  out <- smoothclust(
    input,
    spatial_coords = spatial_coords,
    bandwidth = 0.11,
    normalize = FALSE,
    log = TRUE,
    pseudo_count = 2
  )

  smoothed <- matrix(c(4, 4, 16), nrow = 1)
  expected <- log2(smoothed + 2)
  log_then_smooth <- matrix(
    c(mean(log2(input[1, 1:2] + 2)),
      mean(log2(input[1, 1:2] + 2)), log2(input[1, 3] + 2)),
    nrow = 1
  )

  expect_equal(as.matrix(out), expected)
  expect_false(isTRUE(all.equal(as.matrix(out), log_then_smooth)))
})

test_that("dense and sparse inputs give equivalent sparse default results", {
  input <- matrix(
    c(0, 4, 2, 0, 3, 6, 8, 1),
    nrow = 2,
    dimnames = list(paste0("gene", 1:2), paste0("spot", 1:4))
  )
  input_sparse <- Matrix::Matrix(input, sparse = TRUE)
  spatial_coords <- cbind(c(0, 1, 3, 6), 0)

  out_dense <- smoothclust(
    input, spatial_coords = spatial_coords, bandwidth = 1.1 / 6
  )
  out_sparse <- smoothclust(
    input_sparse, spatial_coords = spatial_coords, bandwidth = 1.1 / 6
  )

  expect_s4_class(out_dense, "dgCMatrix")
  expect_s4_class(out_sparse, "dgCMatrix")
  expect_equal(as.matrix(out_sparse), as.matrix(out_dense))
  expect_identical(dimnames(out_sparse), dimnames(input_sparse))
})

test_that("SpatialExperiment output uses automatic assay names", {
  counts <- matrix(
    c(2, 6, 4, 2, 8, 4, 12, 8),
    nrow = 2,
    dimnames = list(paste0("gene", 1:2), paste0("spot", 1:4))
  )
  spatial_coords <- cbind(x = c(0, 1, 3, 6), y = 0)
  rownames(spatial_coords) <- colnames(counts)
  spe <- SpatialExperiment::SpatialExperiment(
    assays = list(counts = counts), spatialCoords = spatial_coords
  )

  cases <- list(
    list(normalize = FALSE, log = FALSE, assay = "counts_smooth"),
    list(normalize = TRUE, log = FALSE, assay = "normcounts_smooth"),
    list(normalize = FALSE, log = TRUE, assay = "logcounts_smooth"),
    list(normalize = TRUE, log = TRUE, assay = "logcounts_smooth")
  )

  for (case in cases) {
    out <- smoothclust(
      spe,
      bandwidth = 1.1 / 6,
      normalize = case$normalize,
      log = case$log
    )
    matrix_out <- smoothclust(
      counts,
      spatial_coords = spatial_coords,
      bandwidth = 1.1 / 6,
      normalize = case$normalize,
      log = case$log
    )

    expect_true(case$assay %in% SummarizedExperiment::assayNames(out))
    expect_equal(
      as.matrix(SummarizedExperiment::assays(out)[[case$assay]]),
      as.matrix(matrix_out)
    )
    expect_identical(SummarizedExperiment::assays(out)$counts, counts)

    provenance <- S4Vectors::metadata(out)$smoothclust[[case$assay]]
    expect_identical(provenance$normalize, case$normalize)
    expect_identical(provenance$log, case$log)
    expect_equal(provenance$size_factor_source,
                 if (case$normalize) "automatic" else "none")
    if (case$normalize) {
      expect_equal(mean(provenance$size_factors), 1)
    } else {
      expect_null(provenance$size_factors)
    }
    if (case$log) {
      expect_equal(provenance$pseudo_count, 1)
    } else {
      expect_true(is.na(provenance$pseudo_count))
    }
  }
})

test_that("SpatialExperiment output supports custom names and provenance", {
  counts <- matrix(
    c(2, 6, 4, 2, 8, 4, 12, 8),
    nrow = 2,
    dimnames = list(paste0("gene", 1:2), paste0("spot", 1:4))
  )
  spatial_coords <- cbind(x = c(0, 1, 3, 6), y = 0)
  rownames(spatial_coords) <- colnames(counts)
  spe <- SpatialExperiment::SpatialExperiment(
    assays = list(counts = counts), spatialCoords = spatial_coords
  )
  supplied <- setNames(c(1, 2, 3, 4), colnames(counts))

  out <- smoothclust(
    spe,
    bandwidth = 1.1 / 6,
    size_factors = supplied,
    output_assay_name = "custom_smooth"
  )

  expect_true("custom_smooth" %in% SummarizedExperiment::assayNames(out))
  provenance <- S4Vectors::metadata(out)$smoothclust$custom_smooth
  expect_equal(provenance$input_assay, "counts")
  expect_equal(provenance$output_assay, "custom_smooth")
  expect_true(provenance$normalize)
  expect_true(provenance$log)
  expect_equal(provenance$size_factor_source, "supplied")
  expect_equal(unname(provenance$size_factors),
               unname(supplied / mean(supplied)))
  expect_identical(names(provenance$size_factors), colnames(counts))
  expect_equal(provenance$pseudo_count, 1)
})

test_that("SpatialExperiment output refuses to overwrite an existing assay", {
  counts <- matrix(1:8, nrow = 2,
                   dimnames = list(NULL, paste0("spot", 1:4)))
  spatial_coords <- cbind(x = 0:3, y = 0)
  rownames(spatial_coords) <- colnames(counts)
  spe <- SpatialExperiment::SpatialExperiment(
    assays = list(counts = counts, logcounts_smooth = counts),
    spatialCoords = spatial_coords
  )

  expect_error(smoothclust(spe, bandwidth = 0.4))
})
