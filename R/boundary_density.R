#' Boundary density metric
#' 
#' Function for boundary density metric
#' 
#' Function to calculate the raw boundary density metric, defined as the average
#' fraction of nearest neighbors per point that are from a different cluster.
#' Raw BD summarizes observed label-transition frequency on the nearest-neighbor
#' graph.
#' 
#' An optional analytic composition adjustment provides the sample-level
#' relative BD. Relative BD is a composition-adjusted summary, calculated by
#' dividing raw BD by its analytic expectation under random relabeling with the
#' observed label counts. The expectation is calculated directly from label
#' counts. Values below, equal to, or above 1 indicate fewer, the same, or more
#' transitions than expected under this fixed-count random-label reference.
#' 
#' 
#' @param spatial_coords Numeric matrix containing spatial coordinates of
#'   points, formatted as nrow = number of points, ncol = 2 (assuming x and y
#'   dimensions). For example, `spatial_coords = spatialCoords(spe)` if using a
#'   \code{SpatialExperiment} object.
#' 
#' @param labels Atomic vector or factor containing cluster labels for each
#'   point. Missing values are not allowed. For example, `labels <-
#'   colData(spe)$label` if using a \code{SpatialExperiment} object.
#' 
#' @param k Number of k nearest neighbors to use in calculation. Default = 6
#'   (e.g. for hexagonal arrangement in 10x Genomics Visium platform).
#' 
#' @param n_threads Number of threads to use for nearest-neighbor searches.
#'   Default = 1.
#' 
#' @param adjust Whether to calculate the composition adjustment. With
#'   \code{FALSE}, only the raw boundary density is returned. With \code{TRUE},
#'   the analytic composition adjustment returns the expectation under random
#'   relabeling preserving observed label counts and the ratio of raw BD to this
#'   expectation (relative BD). Default = \code{FALSE}.
#' 
#' 
#' @return Returns a list containing values at each point (i.e. the number of
#'   nearest neighbors that are from a different cluster, the number of nearest
#'   neighbors, and the local raw boundary density) as well as the mean
#'   discordant neighbor count and the sample-level raw boundary density. If
#'   \code{adjust = TRUE}, the list also contains the expected boundary
#'   density and sample-level relative boundary density.
#' 
#' 
#' @importFrom BiocNeighbors findKNN
#' 
#' @export
#' 
#' @examples
#' spatial_coords <- cbind(c(0, 1, 3, 6), 0)
#' labels <- c("A", "B", "A", "A")
#' 
#' # calculate raw boundary density metric
#' res <- boundary_density(spatial_coords, labels, k = 2)
#' str(res)
#' res$n_discordant
#' res$mean_discordant
#' res$boundary_density
#' 
#' # calculate analytic composition adjustment
#' res_adj <- boundary_density(spatial_coords, labels, k = 2, 
#'                             adjust = TRUE)
#' res_adj$expected_boundary_density
#' res_adj$relative_boundary_density
#' 
boundary_density <- function(spatial_coords, labels, k = 6, n_threads = 1, 
                             adjust = FALSE) {
  
  stopifnot(is.logical(adjust) && length(adjust) == 1 && !is.na(adjust))
  
  stopifnot(!is.null(spatial_coords), 
            is.numeric(spatial_coords), 
            is.matrix(spatial_coords), 
            ncol(spatial_coords) == 2, 
            all(is.finite(spatial_coords)))
  stopifnot((is.atomic(labels) || is.factor(labels)) && 
              is.null(dim(labels)))
  stopifnot(length(labels) == nrow(spatial_coords))
  stopifnot(!anyNA(labels))
  stopifnot(is.numeric(k) && length(k) == 1 && is.finite(k) && 
              k > 0 && k == floor(k) && k < nrow(spatial_coords))
  stopifnot(is.numeric(n_threads) && length(n_threads) == 1 && 
              is.finite(n_threads) && n_threads > 0 && 
              n_threads == floor(n_threads))
  # --- fast k-nearest neighbor search ---
  
  # search for k nearest neighbors
  nn_data <- findKNN(spatial_coords, k = k, 
                     get.index = TRUE, get.distance = FALSE, 
                     num.threads = n_threads)
  neigh <- nn_data$index
  
  # --- vectorized label lookup and calculation ---
  
  # create matrix of neighbor labels using a single matrix-indexing operation
  neigh_labels <- matrix(labels[neigh], ncol = ncol(neigh))
  
  # compare 'labels' vector against each column of 'neigh_labels'
  n_discordant <- rowSums(labels != neigh_labels)
  n_neighbors <- rep.int(ncol(neigh), nrow(neigh))
  local_boundary_density <- n_discordant / n_neighbors
  boundary_density <- mean(local_boundary_density)
  
  # --- return results ---
  
  # return local and sample-level values
  out <- list(n_discordant = n_discordant, 
              n_neighbors = n_neighbors, 
              local_boundary_density = local_boundary_density, 
              mean_discordant = mean(n_discordant), 
              boundary_density = boundary_density)
  
  if (adjust) {
    label_counts <- as.numeric(table(labels))
    n <- length(labels)
    expected_boundary_density <- 
      1 - sum(label_counts * (label_counts - 1)) / (n * (n - 1))
    relative_boundary_density <- 
      if (expected_boundary_density > 0) {
        boundary_density / expected_boundary_density
      } else {
        NA_real_
      }
    out <- c(out, 
             list(expected_boundary_density = expected_boundary_density, 
                  relative_boundary_density = relative_boundary_density))
  }
  
  out
}
