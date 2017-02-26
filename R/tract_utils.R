#' Tract Constructor
#'
#' \code{tract} is the constructor for objects of class \code{tract}.
#'
#' @param ... A set of name-value pairs. Arguments are evaluated sequentially,
#'   so you can refer to previously created variables. To be a valid tract, the
#'   set should contain at least the fields \code{name} with the name of the
#'   tract, \code{pid} with the identifier of the subject to which the tract
#'   belongs, \code{side} with the hemisphere to which the tract belongs and
#'   \code{data} with the list of \code{streamline}s composing the tract.
#'
#' @return A \code{\link{tract}}.
#' @export
#'
#' @examples
#' tr <- tract(name = "CST", pid = "H018947", side = "Left", data = list())
tract <- function(...) {
  as_tract(tibble::lst(...))
}

#' @export
#' @rdname tract
tract_ <- function(xs) {
  as_tract(tibble::lst_(xs))
}

#' Tract Coercion
#'
#' \code{as_tract} coerces an input \code{list} into a \code{tract}.
#'
#' @param list Input list.
#' @param validate A boolean specifying whether the class of the input object
#'   should be checked (default: \code{TRUE}).
#'
#' @return A \code{\link{tract}}.
#' @export
#'
#' @examples
#' file <- system.file("extdata", "Case001_CST_Left.csv", package = "fdatractography")
#' cst_left <- read_tract(file)
#' tr <- as_tract(cst_left)
as_tract <- function(list, validate = TRUE) {
  vars <- stringr::str_to_lower(names(list))
  names(list) <- vars

  if (validate) {
    if (!is.list(list))
      stop("The input object should be a list.")

    if (!all(c("name", "pid", "side", "data") %in% vars))
      stop("Input list should contain fields name, pid, side and data.")

    if (!is.character(list$name) | length(list$name) != 1L)
      stop("The name fieled should be a character vector of length 1.")

    if (!is.character(list$pid) | length(list$pid) != 1L)
      stop("The pid field should be a character vector of length 1.")

    if (!is.character(list$side) | length(list$side) != 1L)
      stop("The side field should be a character vector of length 1.")

    data_ok <- is.list(list$data)

    if (!data_ok)
      stop("The data field should be a list.")

    for (i in seq_along(list$data)) {
      if (!is_streamline(list$data[[i]])) {
        data_ok <- FALSE
        break
      }
    }

    if (!data_ok)
      stop("The data field should only contain streamline objects.")
  }

  res <- list(name = list$name, pid = list$pid, side = list$side, data = list$data)

  class(res) <- c("tract", class(res))

  res
}

#' Tract Format Verification
#'
#' \code{is_tract} check whether an input R object is of class \code{tract}.
#'
#' @param object An input R object.
#'
#' @return A boolean which is \code{TRUE} if the input object is of class
#'  \code{tract} and \code{FALSE} otherwise.
#' @export
#'
#' @examples
#' file <- system.file("extdata", "Case001_CST_Left.csv", package = "fdatractography")
#' cst_left <- read_tract(file)
#' is_tract(cst_left)
is.tract <- function(object) {
  "tract" %in% class(object)
}

#' @export
#' @rdname is.tract
is_tract <- is.tract

#' Tract Reader
#'
#' Read tractography data stored in a CSV file as output by VTPtoCSV.
#'
#' @param path A string indicating the path to the CSV file to be read. The CSV
#'   filename should be formatted as \code{SubjectID_TractName.csv}.
#'
#' @return A \code{tract} object.
#' @export
#'
#' @examples
#' file <- system.file("extdata", "Case001_CST_Left.csv", package = "fdatractography")
#' cst_left <- read_tract(file)
read_tract <- function(path) {
  filename <- path %>%
    stringr::str_split("/") %>%
    purrr::map_chr(dplyr::last) %>%
    stringr::str_replace(".csv", "") %>%
    stringr::str_split("_")

  name <- filename %>% purrr::map_chr(2L)
  pid <- filename %>% purrr::map_chr(1L)
  side <- filename %>% purrr::map_chr(3L)

  data <- readr::read_csv(path)

  if (nrow(data) == 0)
    return(tract(name = name, pid = pid, side = side, data = list()))

  if (!all(c("LineID", "LPointID", "X", "Y", "Z") %in% names(data)))
    stop("The data tibble should contain at least the variables Side, LineID, LPointID, X, Y and Z.")

  data <- data %>%
    dplyr::arrange(LineID, LPointID) %>%
    dplyr::group_by(LineID) %>%
    dplyr::do(streamlines = streamline(x = .$X, y = .$Y, z = .$Z)) %>%
    dplyr::ungroup()

  tract(name = name, pid = pid, side = side, data = data$streamlines)
}

#' Tract Writer
#'
#' Write tractography data in a CSV file formatted as required by CSVtoVTP.
#'
#' @param tract A \code{\link{tract}}.
#' @param folder A string giving the path to the folder where the CSV file
#'   should be saved.
#'
#' @return TO DO.
#' @export
write_tract <- function(tract, folder) {
  if (!is_tract(tract))
    stop("The input object is not a tract.")

  path <- stringr::str_c(tract$pid, tract$name, tract$side, sep = "_") %>%
    stringr::str_c(folder, ., ".csv")

  data <- tract %>%
    tibble::as_tibble() %>%
    dplyr::mutate(
      LPointID = purrr::map(data, ~ seq_len(nrow(.)) - 1L),
      LineID = seq_len(nrow(.)) - 1L
    ) %>%
    tidyr::unnest() %>%
    dplyr::select(x, y, z, LPointID, LineID)

  readr::write_csv(data, path)
}

#' Tract Plotting
#'
#' \code{plot_tract} plots a \code{tract} object.
#'
#' @param tract A \code{\link{tract}} to be plotted.
#'
#' @return Invisibly returns the input \code{\link{tract}}.
#' @export
#'
#' @examples
#' file <- system.file("extdata", "Case001_CST_Left.csv", package = "fdatractography")
#' cst_left <- read_tract(file)
#' file <- system.file("extdata", "Case001_CST_Right.csv", package = "fdatractography")
#' cst_right <- read_tract(file)
#' cst <- bind_tracts(cst_left, cst_right)
#' plot_tract(cst)
plot_tract <- function(tract) {
  if (is_tract(tract))
    tract <- tibble::as_tibble(tract)
  else if (tibble::is_tibble(tract)) {
    if (!all(c("name", "pid", "side", "data") %in% names(tract)))
      stop("The input tibble does not contain all the required fields.")
  } else
    stop("Input object is neither of class tract nor of class tibble.")

  name <- unique(tract$name)
  if (length(name) > 1L)
    stop("It is not possible to simulaneously plot different tracts.")

  pid <- unique(tract$pid)
  if (length(name) > 1L)
    stop("It is not possible to simulaneously plot a tract for different subjects.")

  print(
    tract %>%
      tidyr::unnest(.id = "LineID") %>%
      tidyr::gather(CoordName, CoordVal, x, y, z) %>%
      ggplot2::ggplot(ggplot2::aes(
        x = s,
        y = CoordVal,
        col = side,
        fill = LineID
      )) +
      ggplot2::geom_line(size = 0.5, alpha = 0.3) +
      ggplot2::facet_wrap(~CoordName, scales = "free") +
      ggplot2::theme_bw() +
      ggplot2::xlab("Arc Length (mm)") +
      ggplot2::ylab("Coordinate (mm)") +
      ggplot2::guides(fill = FALSE) +
      ggplot2::theme(legend.position = "top") +
      ggplot2::ggtitle(paste(name, "of Subject", pid))
  )

  invisible(tract)
}

#' Tract Combination
#'
#' @param ... A set of \code{\link{tract}}s.
#'
#' @return A \code{\link[tibble]{tibble}} with the following 4 variables:
#'   \code{name} with the name of the tract, \code{pid} with the identifier of
#'   the subject to which the tract belongs, \code{side} with the hemisphere to
#'   which the tract belongs and \code{data} with the list of \code{streamline}s
#'   composing the tract.
#' @export
#'
#' @examples
#' file <- system.file("extdata", "Case001_CST_Left.csv", package = "fdatractography")
#' cst_left <- read_tract(file)
#' file <- system.file("extdata", "Case001_CST_Right.csv", package = "fdatractography")
#' cst_right <- read_tract(file)
#' bind_tracts(cst_left, cst_right)
bind_tracts <- function(...) {
  tract_list <- list(...)

  res <- NULL

  for (i in seq_along(tract_list)) {
    tract <- tract_list[[i]]
    if (!is_tract(tract))
      stop("One of the inputs is not a tract.")
    res <- res %>% dplyr::bind_rows(tibble::as_tibble(tract))
  }

  res
}

#' Tract Reparametrization
#'
#' @param tract A \code{\link{tract}}.
#' @param grid_length Size of the grid for uniformly distributed curvilinear
#'   abscissa (default: uses the mean grid size across streamlines composing the
#'   tract).
#'
#' @return A \code{\link{tract}} with streamlines evaluated on the same uniform
#'   grid.
#' @export
#'
#' @examples
#' file <- system.file("extdata", "Case001_CST_Left.csv", package = "fdatractography")
#' cst_left <- read_tract(file)
#' file <- system.file("extdata", "Case001_CST_Right.csv", package = "fdatractography")
#' cst_right <- read_tract(file)
#' reparametrize_tract(cst_left)
reparametrize_tract <- function(tract, grid_length = NULL) {
  if (is.null(grid_length)) {
    grid_length <- tract$data %>%
      purrr::map_int(nrow) %>%
      mean() %>%
      round()
  }

  tract$data <- tract$data %>%
    purrr::map(dplyr::do, tibble::tibble(
      s = modelr::seq_range(.$s, n = grid_length),
      x = approx(.$s, .$x, xout = s)$y,
      y = approx(.$s, .$y, xout = s)$y,
      z = approx(.$s, .$z, xout = s)$y
    ))

  tract
}
