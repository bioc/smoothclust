#' smoothclust
#' 
#' Method for identification of spatial domains and spatially-aware clustering.
#' 
#' Method for identification of spatial domains and spatially-aware clustering
#' in spatial transcriptomics data.
#' 
#' Method for identification of spatial domains and spatially-aware clustering
#' in spatial transcriptomics data. The method generates spatial domains with
#' smooth boundaries by smoothing gene expression profiles across neighboring
#' spatial locations, followed by unsupervised clustering. Spatial domains
#' consisting of consistent mixtures of cell types may then be further
#' investigated by applying cell type compositional analyses or differential
#' analyses.
#'
#' When \code{normalize = TRUE}, library size normalization is performed before
#' spatial smoothing. This ensures that neighboring locations contribute
#' according to the spatial weights rather than their sequencing depths. When
#' \code{log = TRUE}, log2 transformation is also performed after smoothing.
#' Size factors should be calculated from the desired reference set of genes
#' before filtering and supplied using \code{size_factors}.
#' 
#' 
#' @param input Input data, which can be provided as either a
#'   \code{SpatialExperiment} object or a numeric matrix. With the default
#'   preprocessing settings, the selected assay or matrix is assumed to contain
#'   non-negative raw expression counts. Pre-normalized or log-transformed
#'   values may also be supplied if \code{normalize} and \code{log} are set
#'   appropriately. Spatial coordinates are obtained from the
#'   \code{spatialCoords} slot for a \code{SpatialExperiment}, or must be
#'   supplied separately using \code{spatial_coords} for a numeric matrix.
#' 
#' @param assay_name For a \code{SpatialExperiment} input object, this argument
#'   specifies the name of the \code{assay} containing the expression values to
#'   be processed and smoothed. With the default preprocessing settings, this
#'   should contain raw expression counts. If an assay containing normalized or
#'   log-transformed values is used, the \code{normalize} and \code{log}
#'   arguments should be set to avoid applying either step twice. This argument
#'   is only used if the input is a \code{SpatialExperiment} object. Default =
#'   \code{counts}.
#' 
#' @param spatial_coords Numeric matrix of spatial coordinates, assumed to
#'   contain x coordinates in first column and y coordinates in second column.
#'   This argument is only used if the input is a numeric matrix.
#' 
#' @param method Method used for smoothing. Options are \code{uniform},
#'   \code{kernel}, and \code{knn}. The \code{uniform} method calculates
#'   unweighted averages across spatial locations within a circular window with
#'   radius \code{bandwidth} at each spatial location, which smooths out spatial
#'   variability as well as sparsity due to sampling variability. The
#'   \code{kernel} method calculates a weighted average using a truncated
#'   exponential kernel applied to Euclidean distances with a length scale
#'   parameter equal to \code{bandwidth}, which provides a more sophisticated
#'   approach to smoothing out spatial variability but may be affected by
#'   sparsity due to sampling variability (especially sparsity at the index
#'   point), and is computationally slower. The \code{knn} method calculates an
#'   unweighted average across the index point and its k nearest neighbors, and
#'   is the fastest method. Default = \code{uniform}.
#' 
#' @param bandwidth Bandwidth parameter for smoothing, expressed as proportion
#'   of width or height (whichever is greater) of tissue area. Only used for
#'   \code{method = "uniform"} or \code{method = "kernel"}. For \code{method =
#'   "uniform"}, the bandwidth represents the radius of a circle, and unweighted
#'   averages are calculated across neighboring points within this circle. For
#'   \code{method = "kernel"}, the averaging is weighted by distances scaled
#'   using a truncated exponential kernel applied to Euclidean distances. For
#'   example, a bandwidth of 0.05 will smooth values across neighbors weighted
#'   by distances scaled using a truncated exponential kernel with length scale
#'   equal to 5% of the width or height (whichever is greater) of the tissue
#'   area. Weights for \code{method = "kernel"} are truncated at small values
#'   for computational efficiency. Default = 0.05.
#' 
#' @param k Number of nearest neighbors parameter for \code{method = "knn"}.
#'   Only used for \code{method == "knn"}. Unweighted averages are calculated
#'   across the index point and its k nearest neighbors. Default = 18 (based on
#'   two layers in honeycomb pattern for 10x Genomics Visium platform).
#' 
#' @param truncate Truncation threshold parameter if \code{method = "kernel"}.
#'   Kernel weights below this value are set to zero for computational
#'   efficiency. Only used for \code{method = "kernel"}. Default = 0.05.
#' 
#' @param n_threads Number of threads to use for nearest-neighbor searches.
#'   Default = 1.
#'
#' @param normalize Whether to perform library size normalization before spatial
#'   smoothing. Expression values in each column are divided by a size factor
#'   centered to have mean one. Default = \code{TRUE}.
#'
#' @param log Whether to apply a log2 transformation after spatial smoothing.
#'   The transformation is \code{log2(x + pseudo_count)}. Default =
#'   \code{TRUE}.
#'
#' @param size_factors Optional positive numeric vector with one value per
#'   spatial location, used for normalization when \code{normalize = TRUE}.
#'   Supplied values are centered to have mean one. If the vector is named, its
#'   names must match the column names of \code{input}, and values are reordered
#'   to match the input columns. If \code{NULL}, size factors are calculated
#'   from column sums of the selected, unsmoothed input assay or matrix. To
#'   calculate size factors from a broader set of genes than those supplied to
#'   \code{smoothclust}, calculate them before gene filtering and provide them
#'   here. When valid size factors are supplied, columns with zero sums in the
#'   selected input are permitted; when \code{NULL}, every input column must
#'   have a positive sum. Default = \code{NULL}.
#'
#' @param pseudo_count Positive numeric scalar added before log2 transformation
#'   when \code{log = TRUE}. The default value of 1 preserves sparse zeros.
#'   Non-default values may produce a dense output matrix. Default = 1.
#'
#' @param output_assay_name Optional name for the output assay when \code{input}
#'   is a \code{SpatialExperiment}. If \code{NULL}, the name is
#'   \code{logcounts_smooth} when \code{log = TRUE},
#'   \code{normcounts_smooth} when \code{normalize = TRUE} and
#'   \code{log = FALSE}, and \code{<assay_name>_smooth} otherwise. An error is
#'   raised if the output assay already exists. This argument is ignored for a
#'   numeric matrix input. Default = \code{NULL}.
#'
#'
#' @return Returns processed and spatially smoothed expression values for use in
#'   downstream analyses. Results are returned either as a
#'   \code{SpatialExperiment} object containing a new assay, or as a numeric
#'   matrix-like object, depending on the input type.
#' 
#' @importFrom SpatialExperiment spatialCoords
#' @importFrom SummarizedExperiment assays 'assays<-' assayNames
#' @importFrom BiocNeighbors findNeighbors findKNN
#' @importFrom Matrix sparseMatrix
#' @importFrom methods is as
#' @importFrom S4Vectors metadata 'metadata<-'
#' @importFrom utils txtProgressBar setTxtProgressBar
#' 
#' @export
#' 
#' @examples
#' set.seed(123)
#' input <- matrix(rpois(60, lambda = 10), nrow = 10)
#' spatial_coords <- cbind(rep(1:3, each = 2), rep(1:2, times = 3))
#' 
#' # normalize, smooth, and log-transform a numeric expression matrix
#' out <- smoothclust(input, spatial_coords = spatial_coords, bandwidth = 0.6)
#' dim(out)
#'
#' # supply size factors calculated before any gene filtering
#' size_factors <- colSums(input)
#' out_precomputed <- smoothclust(input, spatial_coords = spatial_coords,
#'                                bandwidth = 0.6,
#'                                size_factors = size_factors)
#' dim(out_precomputed)
#' 
smoothclust <- function(input, assay_name = "counts", spatial_coords = NULL, 
                        method = c("uniform", "kernel", "knn"), 
                        bandwidth = 0.05, k = 18, truncate = 0.05, 
                        n_threads = 1, normalize = TRUE, log = TRUE,
                        size_factors = NULL, pseudo_count = 1,
                        output_assay_name = NULL) {
  
  method <- match.arg(method, c("uniform", "kernel", "knn"))
  
  stopifnot(is.character(assay_name) && length(assay_name) == 1 &&
              !is.na(assay_name) && nzchar(assay_name))
  stopifnot(is.numeric(bandwidth) && length(bandwidth) == 1 && 
              is.finite(bandwidth) && bandwidth > 0)
  stopifnot(is.numeric(k) && length(k) == 1 && is.finite(k) && 
              k > 0 && k == floor(k))
  stopifnot(is.numeric(truncate) && length(truncate) == 1 && 
              is.finite(truncate) && truncate > 0 && truncate < 1)
  stopifnot(is.numeric(n_threads) && length(n_threads) == 1 && 
              is.finite(n_threads) && n_threads > 0 && 
              n_threads == floor(n_threads))
  stopifnot(is.logical(normalize) && length(normalize) == 1 && !is.na(normalize))
  stopifnot(is.logical(log) && length(log) == 1 && !is.na(log))
  stopifnot(is.numeric(pseudo_count) && length(pseudo_count) == 1 &&
              is.finite(pseudo_count) && pseudo_count > 0)
  stopifnot(is.null(output_assay_name) ||
              (is.character(output_assay_name) &&
                 length(output_assay_name) == 1 &&
                 !is.na(output_assay_name) && nzchar(output_assay_name)))

  if (!normalize && !is.null(size_factors)) {
    stop("'size_factors' can only be supplied when 'normalize = TRUE'.")
  }
  
  if (is(input, "SpatialExperiment")) {
    spe <- input
    stopifnot(assay_name %in% assayNames(spe))
    if (assay_name == "logcounts" && (normalize || log)) {
      stop("Set 'normalize = FALSE' and 'log = FALSE' when smoothing a ",
           "'logcounts' assay.")
    }
    if (assay_name == "normcounts" && normalize) {
      stop("Set 'normalize = FALSE' when using a 'normcounts' assay.")
    }
    vals <- assays(spe)[[assay_name]]
    spatial_coords <- spatialCoords(spe)
  } else {
    stopifnot((is.numeric(input) || is(input, "Matrix")),
              length(dim(input)) == 2)
    vals <- input
  }
  
  stopifnot(!is.null(spatial_coords), 
            is.numeric(spatial_coords), 
            is.matrix(spatial_coords), 
            ncol(spatial_coords) == 2, 
            all(is.finite(spatial_coords)))
  
  # convert vals to CsparseMatrix for efficient multiplication
  vals <- as(vals, "CsparseMatrix")
  N <- ncol(vals)
  stopifnot(nrow(spatial_coords) == N)

  # validate and normalize expression values before spatial smoothing
  size_factor_source <- "none"
  size_factors_used <- NULL
  if (normalize || log) {
    if (!all(is.finite(vals))) {
      stop("Expression values must be finite when normalization or logging is requested.")
    }
    if (any(vals < 0)) {
      stop("Expression values must be non-negative when normalization or logging is requested.")
    }
  }

  if (normalize) {
    if (is.null(size_factors)) {
      input_library_sizes <- as.numeric(Matrix::colSums(vals))
      if (any(input_library_sizes <= 0)) {
        stop("All input columns must have positive sums when normalization is requested.")
      }
      size_factors_used <- input_library_sizes
      size_factor_source <- "automatic"
    } else {
      if (!is.numeric(size_factors) || !is.null(dim(size_factors)) ||
          length(size_factors) != N || any(!is.finite(size_factors)) ||
          any(size_factors <= 0)) {
        stop("'size_factors' must be a positive finite numeric vector with ",
             "one value per input column.")
      }

      size_factors_used <- as.numeric(size_factors)
      factor_names <- names(size_factors)
      input_names <- colnames(input)
      if (!is.null(factor_names)) {
        if (is.null(input_names)) {
          stop("Named 'size_factors' require non-NULL input column names.")
        }
        if (anyNA(factor_names) || any(!nzchar(factor_names)) ||
            anyDuplicated(factor_names) || anyDuplicated(input_names) ||
            !setequal(factor_names, input_names)) {
          stop("Names of 'size_factors' must uniquely match input column names.")
        }
        size_factors_used <- size_factors_used[match(input_names, factor_names)]
      }
      size_factor_source <- "supplied"
    }

    size_factors_used <- size_factors_used / mean(size_factors_used)
    names(size_factors_used) <- colnames(input)
    vals <- vals %*% Matrix::Diagonal(x = 1 / size_factors_used)
  }

  if (method == "knn") {
    stopifnot(k < N)
  }
  
  if (method %in% c("uniform", "kernel")) {
    # convert bandwidth to same units as distances
    range_x <- abs(diff(range(spatial_coords[, 1])))
    range_y <- abs(diff(range(spatial_coords[, 2])))
    range_max <- max(range_x, range_y)
    bandwidth_scaled <- bandwidth * range_max
  }
  
  # --- method-specific steps to construct weights matrix (W) ---
  
  if (method == "uniform") {
    # 1. fast neighbor search
    nn_data <- findNeighbors(spatial_coords, threshold = bandwidth_scaled, 
                             get.index = TRUE, get.distance = FALSE, 
                             num.threads = n_threads)
    nn_list <- Map(c, seq_along(nn_data$index), nn_data$index)
    
    # 2. get number of neighbors for each spot
    n_neighbors <- lengths(nn_list)
    
    # 3. prepare indices and values for sparse weights matrix W
    # row indices: the neighbors themselves
    i_idx <- unlist(nn_list, use.names = FALSE)
    # column indices: the spot being considered (repeated for each of its neighbors)
    j_idx <- rep(seq_along(nn_list), n_neighbors)
    # values: 1 / (number of neighbors for that column)
    x_val <- rep(1 / n_neighbors, n_neighbors)
    
  } else if (method == "kernel") {
    # 1. determine max search distance based on truncation threshold
    # to avoid slow all-vs-all distance calculation
    max_dist <- -bandwidth_scaled * base::log(truncate)
    
    # 2. fast neighbor search within this radius
    nn_data <- findNeighbors(spatial_coords, threshold = max_dist, 
                             get.index = TRUE, get.distance = TRUE, 
                             num.threads = n_threads)
    
    # 3. calculate raw exponential kernel weights
    nn_list <- Map(c, seq_along(nn_data$index), nn_data$index)
    dist_list <- Map(c, 0, nn_data$distance)
    i_idx <- unlist(nn_list, use.names = FALSE)
    j_idx <- rep(seq_along(nn_list), lengths(nn_list))
    dists <- unlist(dist_list, use.names = FALSE)
    
    raw_weights <- if (bandwidth_scaled == 0) {
      rep(1, length(dists))
    } else {
      exp(-dists / bandwidth_scaled)
    }
    
    # 4. normalize weights so each column in W sums to 1
    col_sums <- as.numeric(sparseMatrix(i = j_idx, 
                                        j = rep.int(1L, length(j_idx)), 
                                        x = raw_weights, dims = c(N, 1)))
    x_val <- raw_weights / col_sums[j_idx]
    
  } else if (method == "knn") {
    # 1. fast k-nearest neighbor search
    nn_data <- findKNN(spatial_coords, k = k, 
                       get.index = "transposed", get.distance = FALSE, 
                       num.threads = n_threads)
    
    # 2. prepare indices and values for W; include self with uniform weight
    nn_mat <- rbind(seq_len(N), nn_data$index)
    n_neighbors <- nrow(nn_mat)
    i_idx <- as.vector(nn_mat)
    j_idx <- rep(seq_len(N), each = n_neighbors)
    x_val <- rep(1 / n_neighbors, length(i_idx))
  }
  
  # --- construct weights matrix W and perform single matrix multiplication ---
  
  # 1. construct weights matrix
  W <- sparseMatrix(i = i_idx, j = j_idx, x = x_val, dims = c(N, N), repr = "C")
  
  # 2. perform smoothing operation in one matrix multiplication step
  vals_smooth <- vals %*% W

  # 3. apply log2 transformation after smoothing
  if (log) {
    if (pseudo_count == 1) {
      vals_smooth <- base::log1p(vals_smooth) / base::log(2)
    } else {
      vals_smooth <- base::log2(vals_smooth + pseudo_count)
    }
  }
  
  # --- return results ---
  
  stopifnot(nrow(vals_smooth) == nrow(input))
  stopifnot(ncol(vals_smooth) == ncol(input))
  
  rownames(vals_smooth) <- rownames(input)
  colnames(vals_smooth) <- colnames(input)
  
  # return results (smoothed values)
  if (is(input, "SpatialExperiment")) {
    if (is.null(output_assay_name)) {
      output_assay_name <- if (log) {
        "logcounts_smooth"
      } else if (normalize) {
        "normcounts_smooth"
      } else {
        paste0(assay_name, "_smooth")
      }
    }
    if (output_assay_name %in% assayNames(spe)) {
      stop("Output assay '", output_assay_name, "' already exists. Use ",
           "'output_assay_name' to choose a different name.")
    }

    assays(spe)[[output_assay_name]] <- vals_smooth

    object_metadata <- metadata(spe)
    if (is.null(object_metadata$smoothclust)) {
      object_metadata$smoothclust <- list()
    }
    object_metadata$smoothclust[[output_assay_name]] <- list(
      input_assay = assay_name,
      output_assay = output_assay_name,
      normalize = normalize,
      log = log,
      size_factor_source = size_factor_source,
      size_factors = size_factors_used,
      pseudo_count = if (log) pseudo_count else NA_real_,
      method = method,
      bandwidth = bandwidth,
      k = k,
      truncate = truncate,
      n_threads = n_threads
    )
    metadata(spe) <- object_metadata
    spe
  } else {
    vals_smooth
  }
}
