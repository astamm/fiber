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

#' @export
#' @rdname streamline
streamline_ <- function(xs, validate = TRUE) {
  as_streamline(tibble::lst_(xs), validate = validate)
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
#' file <- system.file("extdata", "Case001_CST_Left.csv", package = "fdatractography")
#' cst_left <- read_tract(file)
#' as_streamline(cst_left$data[[1]])
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
#' file <- system.file("extdata", "Case001_CST_Left.csv", package = "fdatractography")
#' cst_left <- read_tract(file)
#' is_streamline(cst_left$data[[1]])
is.streamline <- function(x) {
  "streamline" %in% class(x)
}

#' @export
#' @rdname is.streamline
is_streamline <- is.streamline

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
#' file <- system.file("extdata", "Case001_CST_Left.csv", package = "fdatractography")
#' cst_left <- read_tract(file)
#' st <- cst_left$data[[1]]
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

  max(streamline$s)
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

#' @rdname get-shape
#' @export
get_curvature <- function(streamline, validate = TRUE, direction = NULL) {
  if (validate) {
    if (!is_streamline(streamline))
      stop("The input dataset is not of class streamline.")
    if (!(is.null(direction) ||
          stringr::str_to_lower(direction) %in% c("x", "y", "z")))
      stop("Possible choices for direction argument are NULL (3D curvature), x, y or z.")
  }

  n <- nrow(streamline)

  if (is.null(direction)) {
    streamline <- streamline %>%
      dplyr::do(dplyr::tibble(
        s = .$s,
        x = approx(x = .$s, y = .$x, n = n)$y,
        y = approx(x = .$s, y = .$y, n = n)$y,
        z = approx(x = .$s, y = .$z, n = n)$y
      )) %>%
      dplyr::mutate_at(
        .vars = dplyr::vars(d2x = x, d2y = y, d2z = z),
        .funs = dplyr::funs(dplyr::lead(.) + dplyr::lag(.) - 2 * .)
      ) %>%
      dplyr::mutate(k = (n - 1) ^ 2 * sqrt(d2x ^ 2 + d2y ^ 2 + d2z ^ 2))
  } else {
    streamline <- streamline %>%
      dplyr::do(dplyr::tibble(s = .$s,
                              direction = approx(
                                x = .$s, y = .[[direction]], n = n
                              )$y)) %>%
      dplyr::mutate(
        d2 = dplyr::lead(direction) + dplyr::lag(direction) - 2 * direction,
        k = (n - 1) ^ 2 * abs(d2)
      )
  }

  streamline %>%
    dplyr::select(s, k) %>%
    tidyr::drop_na()
}

#' @rdname get-shape
#' @export
get_curvature_max <- function(streamline, validate = TRUE, direction = NULL) {
  if (validate) {
    if (!is_streamline(streamline))
      stop("The input dataset is not of class streamline.")
    if (!(is.null(direction) ||
          stringr::str_to_lower(direction) %in% c("x", "y", "z")))
      stop("Possible choices for direction argument are NULL (3D curvature), x, y or z.")
  }

  curvature <- get_curvature(streamline, FALSE, direction)
  max(curvature$k)
}

#' @rdname get-shape
#' @export
get_curvature_mean <- function(streamline, validate = TRUE, direction = NULL) {
  if (validate) {
    if (!is_streamline(streamline))
      stop("The input dataset is not of class streamline.")
    if (!(is.null(direction) ||
          stringr::str_to_lower(direction) %in% c("x", "y", "z")))
      stop("Possible choices for direction argument are NULL (3D curvature), x, y or z.")
  }

  curvature <- get_curvature(streamline, FALSE, direction)
  mean(curvature$k)
}

#' @rdname get-shape
#' @export
get_curvature_sd <- function(streamline, validate = TRUE, direction = NULL) {
  if (validate) {
    if (!is_streamline(streamline))
      stop("The input dataset is not of class streamline.")
    if (!(is.null(direction) ||
          stringr::str_to_lower(direction) %in% c("x", "y", "z")))
      stop("Possible choices for direction argument are NULL (3D curvature), x, y or z.")
  }

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
#' file <- system.file("extdata", "Case001_CST_Left.csv", package = "fdatractography")
#' cst_left <- read_tract(file)
#' tmp <- reparametrize_tract(cst_left, grid_length = 50L)
#' str1 <- tmp$data[[1]]
#' str2 <- tmp$data[[2]]
#' get_hausdorff_distance(str1, str2)
#' get_L2_distance(str1, str2)
#' get_L1_distance(str1, str2)
NULL

#' @rdname get-distance
#' @export
get_hausdorff_distance <- function(streamline1, streamline2) {
  tmp <- streamline1 %>%
    dplyr::mutate(
      Points = purrr::pmap(list(x, y, z), c),
      dist = purrr::map_dbl(Points, get_hausdorff_distance_internal, streamline2)
    ) %>%
    dplyr::summarise(dist = max(dist))
  max_dist1 <- tmp$dist

  tmp <- streamline2 %>%
    dplyr::mutate(
      Points = purrr::pmap(list(x, y, z), c),
      dist = purrr::map_dbl(Points, get_hausdorff_distance_internal, streamline1)
    ) %>%
    dplyr::summarise(dist = max(dist))
  max_dist2 <- tmp$dist

  max(max_dist1, max_dist2)
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

align <- function(fixed_streamline, moving_streamline, cost_function) {
  xfun <- approxfun(moving_streamline$s, moving_streamline$x)
  yfun <- approxfun(moving_streamline$s, moving_streamline$y)
  zfun <- approxfun(moving_streamline$s, moving_streamline$z)
  optimize(f = cost_function, interval = c(0, 2), str = fixed_streamline,
           xfun = xfun, yfun = yfun, zfun = zfun)
}

align_streamline <- function(fixed_streamline, moving_streamline, cost_function) {
  opt <- align(fixed_streamline, moving_streamline, cost_function)
  streamline(
    s = moving_streamline$s / opt$minimum,
    x = moving_streamline$x,
    y = moving_streamline$y,
    z = moving_streamline$z,
    validate = FALSE
  )
}
