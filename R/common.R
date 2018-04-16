#' @importFrom dplyr %>%
#' @export
dplyr::`%>%`

#' @importFrom rlang :=
#' @export
rlang::`:=`

#' @importFrom rlang .data
#' @export
rlang::.data

#' @importFrom roahd MBD
#' @export
roahd::MBD

get_hausdorff_distance_impl <- function(point, streamline) {
  streamline %>%
    dplyr::mutate(
      distance = purrr::pmap(list(x, y, z), c) %>%
        purrr::map_dbl(~ sqrt(sum((. - point)^2)))
    ) %>%
    dplyr::summarise(dist = min(dist)) %>%
    dplyr::pull(distance)
}

#' Tract Clustering
#'
#' @param distance_matrix A matrix of size \code{n x n} where \code{n} is the
#'   number of streamlines in the tract and each entry of the matrix evaluates
#'   the distance between two streamlines.
#' @param k An integer specifying the number of clusters (default: \code{1L}).
#' @param max_iter An integer specifying the maximum number of iterations
#'   (default: \code{100L}).
#' @param nstart An integer specifying the number of random initializations
#'   (default: \code{100L}).
#' @param ncores An integer specifying the number of cores available for the computations (default: \code{1L}).
#'
#' @return A list with two components: \code{groups} gives membership of each
#'   individual streamline while \code{centroids} gives the labels of the
#'   groups' centroids.
#' @export
#'
#' @examples
#' D <- matrix(runif(100), 10L, 10L)
#' D <- (D + t(D)) / 2
#' diag(D) <- 0
#' find_clusters(D, k = 4L)
find_clusters <- function(distance_vector, k = 1L, max_iter = 100L, nstart = 100L, ncores = 1L) {
  distance_matrix <- get_distance_matrix(distance_vector)
  ncores <- min(ncores, nstart)
  parallel <- (ncores > 1L && requireNamespace("multidplyr", quietly = TRUE))

  if (parallel) {
    cl <- multidplyr::create_cluster(cores = ncores) %>%
      multidplyr::cluster_copy(find_clusters_internal) %>%
      multidplyr::cluster_copy(distance_matrix) %>%
      multidplyr::cluster_copy(k) %>%
      multidplyr::cluster_copy(max_iter)
    tmp <- tibble::tibble(id = seq_len(nstart)) %>%
      multidplyr::partition(cluster = cl) %>%
      dplyr::mutate(
        res = purrr::map(id, ~ find_clusters_internal(distance_matrix, k, max_iter))
      ) %>%
      dplyr::collect() %>%
      dplyr::ungroup() %>%
      dplyr::arrange(id)
    res <- tmp$res
  } else {
    res <- list()
    for (i in seq_len(nstart))
      res[[i]] <- find_clusters_internal(distance_matrix, k, max_iter)
  }
  idx <- which.min(purrr::map_dbl(res, "wmse"))
  res[[idx]]
}

find_clusters_internal <- function(distance_matrix, k = 1L, max_iter = 100L) {
  n <- nrow(distance_matrix)
  labels <- seq_len(n)
  init_clusters <- sample.int(n, k)
  within <- rep(0, k)
  final_wmse <- 0
  for (iter in seq_len(max_iter)) {
    # Assign curves to clusters
    groups <- labels %>%
      purrr::map_int(~ which.min(distance_matrix[., init_clusters, drop = TRUE]))
    # Find new cluster centers
    for (i in seq_len(k)) {
      group <- (groups == i)
      group_labels <- labels[group]
      distances <- distance_matrix[group_labels, group_labels, drop = FALSE]
      min_norm <- 0
      for (j in seq_along(group_labels)) {
        tmp_distances <- distances[-j, -j, drop = FALSE]
        tmp_norm <- norm(tmp_distances, type = "F")
        if (tmp_norm < min_norm || j == 1) {
          init_clusters[i] <- group_labels[j]
          min_norm <- tmp_norm
        }
      }
      within[i] <- min_norm
    }
    wmse <- mean(within)
    if (wmse < final_wmse || iter == 1L) {
      final_groups <- groups
      final_centroids <- init_clusters
      final_within <- within
      final_wmse <- wmse
    } else break
  }
  list(groups = final_groups, cendroids = final_centroids,
       within = final_within, wmse = final_wmse)
}

get_distance_matrix <- function(distance_vector) {
  p <- length(distance_vector)
  n <- (1 + sqrt(1 + 8 * p)) / 2
  res <- matrix(0, n, n)
  offset <- 0
  for (i in 1:(n-1)) {
    res[i, (i+1):n] <- distance_vector[offset + 1:(n-i)]
    offset <- offset + (n - i)
  }
  res[lower.tri(res)] <- t(res)[lower.tri(res)]
  res
}

depth_cost <- function(x, depth_data) {
  (mean(depth_data >= x * max(depth_data)) - 0.5)^2
}

alignment_cost <- function(param, str, xfun, yfun, zfun) {
  a <- param[2]
  b <- param[1]
  mean(
    abs(str$x - xfun(a * str$s + b))  +
      abs(str$y - yfun(a * str$s + b)) +
      abs(str$z - zfun(a * str$s + b)),
    na.rm = TRUE
  )
}

is_mfData <- function(x) {
  "mfData" %in% class(x)
}
