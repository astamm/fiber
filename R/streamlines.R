#' Streamline Constructor
#'
#' \code{streamline} is the constructor for objects of class \code{streamline}.
#'
#' @param ... A set of name-value pairs. Arguments are evaluated sequentially,
#'   so you can refer to previously created variables. To be a valid tract, the
#'   set should contain at least the fields x, y, z.
#' @param validate A boolean specifying whether the class of the input object
#'   should be checked (default: \code{TRUE}).
#'
#' @return A \code{\link{streamline}}.
#' @export
#' @examples
#' st <- streamline(X = 1:10, Y = 1:10, Z = 1:10, Pie = 1:10)
streamline <- function(..., validate = TRUE) {
  as_streamline(tibble::lst(...), validate = validate)
}

#' Streamline Coercion
#'
#' Coerces a dataset into a \code{\link{streamline}}.
#'
#' @param streamline Input dataset.
#' @param validate A boolean specifying whether the class of the input object
#'   should be checked (default: \code{TRUE}).
#'
#' @return A \code{\link{streamline}}.
#' @export
#'
#' @examples
#' file <- system.file("extdata", "CC_SingleTensor.csv", package = "fiber")
#' cc <- read_tract(file)
#' as_streamline(cc$data[[1]])
as_streamline <- function(input, validate = TRUE, ...) {
  names(input) <- stringr::str_to_lower(names(input))
  vars <- names(input)

  if (validate) {
    if (!all(c("x", "y", "z") %in% vars))
      stop("The input list should have at least x, y and z entries.")
  }

  input <- tibble::as_tibble(input, validate = validate, ...)

  if (!("s" %in% vars)) {
    input <- input %>%
      dplyr::mutate_at(
        .vars = dplyr::vars(dx = x, dy = y, dz = z),
        .funs = dplyr::funs(. - lag(.))
      ) %>%
      tidyr::replace_na(list(dx = 0, dy = 0, dz = 0)) %>%
      dplyr::mutate(norm = sqrt(dx^2 + dy^2 + dz^2), s = cumsum(norm)) %>%
      dplyr::select(-dx, -dy, -dz, -norm)
  }

  input <- input %>% dplyr::select(s, x, y, z, dplyr::everything())

  class(input) <- c("streamline", class(input))

  input
}

#' Streamline Format Verification
#'
#' Check whether the input object is of class \code{streamline}.
#'
#' @param streamline An input R object.
#'
#' @return A boolean which is \code{TRUE} if the input object is of class
#'  \code{streamline} and \code{FALSE} otherwise.
#' @export
#'
#' @examples
#' file <- system.file("extdata", "CC_SingleTensor.csv", package = "fiber")
#' cc <- read_tract(file)
#' is_streamline(cc$data[[1]])
is_streamline <- function(x) {
  "streamline" %in% class(x)
}

#' Shape Characteristics of a Streamline
#'
#' \code{get_euclidean_length} computes the Euclidean length of a streamline,
#' which is the Euclidean distance between its endpoints.
#' \code{get_curvilinear_length} computes the curvilinear length of a
#' streamline, which is the range of its curvilinear abscissa.
#' \code{get_sinuosity} computes the sinuosity of a streamline, which is the
#' ratio between its curvilinear and Euclidean lengths. \code{get_curvature}
#' computes the curvature of a streamline along its curvilinear abscissa.
#' \code{get_curvature_max} computes the maximal curvature of a streamline along
#' its curvilinear abscissa. \code{get_curvature_mean} computes the mean
#' curvature of a streamline along its curvilinear abscissa.
#' \code{get_curvature_sd} computes the standard deviation of the curvature of a
#' streamline along its curvilinear abscissa.
#'
#' @param streamline A \code{\link{streamline}}.
#' @param validate A boolean which is \code{TRUE} if the input type should be
#'   checked or \code{FALSE} otherwise (default: \code{TRUE}).
#' @param direction String for choosing along which axis to compute the
#'   curvature (default: \code{NULL} which computes the 3D curvature)
#'
#' @return A scalar giving the desired shape characteristic of the input
#'   streamline, except for \code{get_curvature} which returns a
#'   \code{\link[tibble]{tibble}}.
#' @name get-shape
#' @examples
#' file <- system.file("extdata", "CC_SingleTensor.csv", package = "fiber")
#' cc <- read_tract(file)
#' st <- cc$data[[1]]
#' get_euclidean_length(st)
#' get_curvilinear_length(st)
#' get_sinuosity(st)
#' get_curvature(st)
#' get_curvature_max(st)
#' get_curvature_mean(st)
#' get_curvature_sd(st)
NULL

#' @rdname get-shape
#' @export
get_euclidean_length <- function(streamline, validate = TRUE) {
  if (validate) {
    if (!is_streamline(streamline))
      stop("The input dataset is not of class streamline.")
  }

  n <- nrow(streamline)
  dx <- streamline$x[n] - streamline$x[1]
  dy <- streamline$x[n] - streamline$y[1]
  dz <- streamline$x[n] - streamline$z[1]

  sqrt(dx^2 + dy^2 + dz^2)
}

#' @rdname get-shape
#' @export
get_curvilinear_length <- function(streamline, validate = TRUE) {
  if (validate) {
    if (!is_streamline(streamline))
      stop("The input dataset is not of class streamline.")
  }

  max(streamline$s) - min(streamline$s)
}

#' @rdname get-shape
#' @export
get_sinuosity <- function(streamline, validate = TRUE) {
  if (validate) {
    if (!is_streamline(streamline))
      stop("The input dataset is not of class streamline.")
  }

  cl <- get_curvilinear_length(streamline, validate = FALSE)
  el <- get_euclidean_length(streamline, validate = FALSE)

  cl / el
}

myeval <- function(model, test) {
  test <- test %>%
    tibble::as_tibble() %>%
    dplyr::filter(s >= model$basis$rangeval[1], s <= model$basis$rangeval[2])
  obs <- as.matrix(test[, 2:4])
  pred <- fda::eval.fd(test$s, model)
  norm(obs - pred, "F")
}

cost_gcv <- function(x, st) {
  cv <- modelr::crossv_mc(st, 100)
  models <- purrr::map(cv$train, cost_fit, x = x)
  errors <- purrr::map2_dbl(models, cv$test, myeval)
  mean(errors)
}

cost_fit <- function(st, x) {
  st <- tibble::as_tibble(st)
  grid <- st$s
  data <- array(dim = c(length(grid), 1L, 3L))
  data <- matrix(0, nrow = length(grid), 3L)
  data[, 1] <- st$x
  data[, 2] <- st$y
  data[, 3] <- st$z
  fda::Data2fd(grid, data, basisobj = 5L, lambda = x)
}

#' @rdname get-shape
#' @export
get_curvature <- function(st, validate = TRUE) {
  if (validate)
    if (!is_streamline(st))
      stop("The input dataset is not of class streamline.")

  grid <- st$s
  data <- array(dim = c(length(grid), 1L, 3L))
  data <- matrix(0, nrow = length(grid), 3L)
  data[, 1] <- st$x
  data[, 2] <- st$y
  data[, 3] <- st$z

  lambda_opt <- optimise(cost_gcv, c(1e-6, 100), st = st)$minimum
  writeLines(paste("Optimal roughness penalization (minimal GCV):", lambda_opt))
  fd <- fda::Data2fd(grid, data, basisobj = 5L, lambda = lambda_opt)
  d1 <- fda::eval.fd(grid, fd, Lfdobj = 1)
  d2 <- fda::eval.fd(grid, fd, Lfdobj = 2)
  dx <- d1[, 1]
  dy <- d1[, 2]
  dz <- d1[, 3]
  d2x <- d2[, 1]
  d2y <- d2[, 2]
  d2z <- d2[, 3]
  k <- sqrt( ( (d2z * dy - d2y * dz)^2 + (d2x * dz - d2z * dx)^2 + (d2y * dx - d2x * dy)^2 )
             / (dx^2 + dy^2 + dz^2)^3 )
  tibble::tibble(s = grid, curv = k)
}

#' @rdname get-shape
#' @export
get_curvature_max <- function(streamline, validate = TRUE) {
  if (validate)
    if (!is_streamline(streamline))
      stop("The input dataset is not of class streamline.")

  curvature <- get_curvature(streamline, FALSE)
  max(curvature$k)
}

#' @rdname get-shape
#' @export
get_curvature_mean <- function(streamline, validate = TRUE) {
  if (validate)
    if (!is_streamline(streamline))
      stop("The input dataset is not of class streamline.")

  curvature <- get_curvature(streamline, FALSE, direction)
  mean(curvature$k)
}

#' @rdname get-shape
#' @export
get_curvature_sd <- function(streamline, validate = TRUE) {
  if (validate)
    if (!is_streamline(streamline))
      stop("The input dataset is not of class streamline.")

  curvature <- get_curvature(streamline, FALSE, direction)
  sd(curvature$k)
}

#' Distances Between Streamlines
#'
#' \code{get_hausdorff_distance} computes the Hausdorff distance between the two input
#'   streamlines.
#' \code{get_L2_distance} computes the L2 distance between the two input
#'   streamlines.
#'
#' @param streamline1 A \code{\link{streamline}}.
#' @param streamline2 A \code{\link{streamline}}.
#'
#' @return A non-negative scalar representing the chosen distance between the two input
#'   streamlines.
#' @name get-distance
#' @examples
#' file <- system.file("extdata", "CC_SingleTensor.csv", package = "fiber")
#' cc <- read_tract(file)
#' tmp <- reparametrize_tract(cc, grid_length = 50L)
#' str1 <- tmp$data[[1]]
#' str2 <- tmp$data[[2]]
#' get_hausdorff_distance(str1, str2)
#' get_L2_distance(str1, str2)
#' get_L1_distance(str1, str2)
NULL

#' @rdname get-distance
#' @export
get_hausdorff_distance <- function(streamline1, streamline2) {
  dist1 <- streamline1 %>%
    dplyr::mutate(
      dist = purrr::pmap(list(x, y, z), c) %>%
        purrr::map_dbl(get_hausdorff_distance_impl, streamline = streamline2)
    ) %>%
    dplyr::summarise(dist = max(dist)) %>%
    dplyr::pull(dist)

  dist2 <- streamline2 %>%
    dplyr::mutate(
      dist = purrr::pmap(list(x, y, z), c) %>%
        purrr::map_dbl(get_hausdorff_distance_impl, streamline = streamline1)
    ) %>%
    dplyr::summarise(dist = max(dist)) %>%
    dplyr::pull(dist)

  max(dist1, dist2)
}

#' @rdname get-distance
#' @export
get_L2_distance <- function(streamline1, streamline2) {
  if (nrow(streamline1) > nrow(streamline2)) {
    str <- streamline1
    streamline1 <- streamline2
    streamline2 <- str
  }
  aln <- align(
    fixed_streamline = streamline1,
    moving_streamline = streamline2,
    cost_function = cost_L2
  )
  aln$objective
}

#' @rdname get-distance
#' @export
get_L1_distance <- function(streamline1, streamline2) {
  if (nrow(streamline1) > nrow(streamline2)) {
    str <- streamline1
    streamline1 <- streamline2
    streamline2 <- str
  }
  aln <- align(
    fixed_streamline = streamline1,
    moving_streamline = streamline2,
    cost_function = cost_L1
  )
  aln$objective
}

align <- function(fixed_streamline, xfun, yfun, zfun, cost_function) {
  optim(
    par = c(0, 1),
    fn = cost_function,
    str = fixed_streamline,
    xfun = xfun,
    yfun = yfun,
    zfun = zfun,
    method = "Nelder-Mead"
  )
}

align_streamline <- function(fixed_streamline, moving_streamline, cost_function) {
  xfun <- approxfun(moving_streamline$s, moving_streamline$x)
  yfun <- approxfun(moving_streamline$s, moving_streamline$y)
  zfun <- approxfun(moving_streamline$s, moving_streamline$z)
  fafun <- approxfun(
    moving_streamline$s,
    purrr::map_dbl(moving_streamline$diffusion, get_fractional_anisotropy, validate = FALSE)
  )
  opt <- align(fixed_streamline, xfun, yfun, zfun, cost_function)
  a <- opt$par[2]
  b <- opt$par[1]
  abscissa <- a * fixed_streamline$s + b
  streamline(
    s = fixed_streamline$s,
    x = xfun(abscissa),
    y = yfun(abscissa),
    z = zfun(abscissa),
    fa = fafun(abscissa),
    validate = FALSE
  )
}

align_streamline2 <- function(fixed_streamline, moving_streamline, cost_function) {
  xfun <- approxfun(moving_streamline$s, moving_streamline$x)
  yfun <- approxfun(moving_streamline$s, moving_streamline$y)
  zfun <- approxfun(moving_streamline$s, moving_streamline$z)
  fafun <- approxfun(moving_streamline$s, moving_streamline$fa)
  opt <- align(fixed_streamline, xfun, yfun, zfun, cost_function)
  a <- opt$par[2]
  b <- opt$par[1]
  print(a)
  print(b)
  abscissa <- a * fixed_streamline$s + b
  streamline(
    s = fixed_streamline$s,
    x = xfun(abscissa),
    y = yfun(abscissa),
    z = zfun(abscissa),
    fa = fafun(abscissa),
    validate = FALSE
  )
}

#' Streamline Plotting
#'
#' This function plots a strewamline in 3D using the \code{rgl} package. If
#' diffusion tensor information is available, it is also plotted as ellipsoids
#' at each point of the streamline.
#'
#' @param streamline A \code{\link{streamline}} object.
#' @param validate A boolean specifying whether the class of the input object
#'   should be checked (default: \code{TRUE}).
#' @param plot_microstructure A boolean specifying if microstruture should be
#'   superimposed on the displayed \code{\link{streamline}}s.
#' @param new_window A boolean specifying whether a new graphical window should
#'   be used (default: \code{TRUE}).
#' @param scale A scalar that handles the ellipsoid scale (default: 4).
#' @param ... Additional plotting parameters to be passed to \code{rgl} basic
#'   plotting functions.
#'
#' @return The function is used for its side effect of plotting.
#' @export
#'
#' @examples
plot_streamline <- function(
  streamline,
  validate = TRUE,
  plot_microstructure = FALSE,
  new_window = TRUE,
  scale = 4,
  ...) {
  if (validate)
    if (!is_streamline(streamline))
      stop("The input object to be plotted should be of class streamline.")

  if (new_window) rgl::open3d()
  rgl::lines3d(x = streamline$x, y = streamline$y, z = streamline$z, ...)

  if (plot_microstructure) {
    if ("t" %in% names(streamline)) {
      size <- get_curvilinear_length(streamline)
      npts <- nrow(streamline)
      scl <- scale * size / npts
      for (i in 1:npts)
        rgl::plot3d(
          rgl::ellipse3d(
            streamline$t[[i]],
            centre = c(streamline$x[[i]], streamline$y[[i]], streamline$z[[i]]),
            scale = rep(scl, 3L)
          ),
          add = TRUE,
          ...
        )
    }
  }
}