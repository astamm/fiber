#' @importFrom dplyr %>%
#' @export
dplyr::`%>%`

#' @importFrom roahd MBD
#' @export
roahd::MBD

#' Vector-To-Tensor Representation
#'
#' @param vector A vector of size 6 storing the 6 unique components of a tensor.
#' @param twice A boolean that says whether off-diagonal tensor elements were
#'   doubled in vector representation (default: \code{FALSE}).
#'
#' @return A 3x3 symmetric postive definite matrix as output tensor.
#' @export
#'
#' @examples
#' v <- seq_len(6L)
#' as_tensor(v)
as_tensor <- function(vector, twice = FALSE) {
  if (!is.numeric(vector))
    stop("Input should be a numeric vector")
  else {
    if (length(vector) != 6L)
      stop("Input should be of dimension 6")
  }

  xx <- vector[1]
  yx <- vector[2]
  yy <- vector[3]
  zx <- vector[4]
  zy <- vector[5]
  zz <- vector[6]

  if (twice) {
    yx <- yx / 2
    zx <- zx / 2
    zy <- zy / 2
  }

  cbind(c(xx, yx, zx), c(yx, yy, zy), c(zx, zy, zz))
}

#' Diffusion Biomarkers
#'
#' \code{get_fractional_anisotropy} computes the fractional anisotropy of a
#' diffusion tensor.
#' \code{get_mean_diffusivity} computes the mean diffusivity of a diffusion tensor.
#'
#' @param tensor A \code{3x3} symmetric definite positive matrix.
#' @param validate A boolean which is \code{TRUE} if the input type should be
#'   checked or \code{FALSE} otherwise (default: \code{TRUE}).
#'
#' @return A scalar giving the desired diffusion biomarker extracted from the input diffusion tensor.
#' @name get-biomarker
#' @examples
#' DT <- diag(c(1.71e-3, 3e-4, 1e-4))
#' get_fractional_anisotropy(DT)
#' get_mean_diffusivity(DT)
NULL

#' @rdname get-biomarker
#' @export
get_fractional_anisotropy <- function(tensor, validate = TRUE) {
  if (validate) {
    print("Do some checking on input tensor...")
  }

  m <- get_mean_diffusivity(tensor, FALSE)
  vals <- eigen(tensor, symmetric = TRUE, only.values = TRUE)$values
  sqrt(1.5 * sum((vals - m)^2) / sum(vals^2))
}

#' @rdname get-biomarker
#' @export
get_mean_diffusivity <- function(tensor, validate = TRUE) {
  if (validate) {
    print("Do some checking on input tensor...")
  }

  mean(diag(tensor))
}

add_diffusion_information <- function(data, model = "None") {
  switch(
    model,
    None = data,
    DTI = data %>%
      dplyr::mutate(
        Tensors = purrr::pmap(list(`Tensors#0`, `Tensors#1`, `Tensors#2`,
                                   `Tensors#3`,`Tensors#4`, `Tensors#5`), c) %>%
          purrr::map(as_tensor)
      ) %>%
      dplyr::select(-dplyr::starts_with("Tensors#"))
  )
}

get_hausdorff_distance_internal <- function(point, streamline) {
  tmp <- streamline %>%
    dplyr::mutate(
      Points = purrr::pmap(list(x, y, z), c),
      dist = purrr::map_dbl(Points, ~ sqrt(sum((. - point)^2)))
    ) %>%
    dplyr::summarise(dist = min(dist))
  tmp$dist
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

cost_L2 <- function(param, str, xfun, yfun, zfun) {
  sqrt(mean((str$x - xfun(param[1] * str$s))^2  +
              (str$y - yfun(param[1] * str$s))^2 +
              (str$z - zfun(param[1] * str$s))^2, na.rm = TRUE))
}

cost_L1 <- function(param, str, xfun, yfun, zfun) {
  mean(abs(str$x - xfun(param[1] * str$s))  +
         abs(str$y - yfun(param[1] * str$s)) +
         abs(str$z - zfun(param[1] * str$s)), na.rm = TRUE)
}
