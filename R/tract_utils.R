#' Tract Constructor
#'
#' \code{tract} is the constructor for objects of class \code{tract}.
#'
#' @param ... A set of name-value pairs. Arguments are evaluated sequentially,
#'   so you can refer to previously created variables. To be a valid tract, the
#'   set should contain at least the fields \code{name} with the name of the
#'   tract, \code{case} with the identifier of the subject to which the tract
#'   belongs, \code{scan} to identify repeated measures for the same subject,
#'   \code{side} with the hemisphere to which the tract belongs and
#'   \code{data} with the list of \code{streamline}s composing the tract.
#'
#' @return A \code{\link{tract}}.
#' @export
#'
#' @examples
#' tr <- tract(name = "CST", case = "H018947", scan = "01", side = "L",
#'   data = list())
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
as_tract <- function(input, ...) {
  UseMethod("as_tract", input)
}

#' @export
#' @rdname as_tract
as_tract.list <- function(input, validate = TRUE, ...) {
  vars <- stringr::str_to_lower(names(input))
  names(input) <- vars

  if (validate) {
    if (!is.list(input))
      stop("The input object should be a list.")

    if (!all(c("name", "case", "scan", "side", "data") %in% vars))
      stop("Input list should contain fields name, case, scan, side and data.")

    if (!is.character(input$name) | length(input$name) != 1L)
      stop("The name field should be a character vector of length 1.")

    if (!is.character(input$case) | length(input$case) != 1L)
      stop("The case field should be a character vector of length 1.")

    if (!is.character(input$scan) | length(input$scan) != 1L)
      stop("The scan field should be a character vector of length 1.")

    if (!is.character(input$side) | length(input$side) != 1L)
      stop("The side field should be a character vector of length 1.")

    data_ok <- is.list(input$data)

    if (!data_ok)
      stop("The data field should be a list.")

    for (i in seq_along(input$data)) {
      if (!is_streamline(input$data[[i]])) {
        data_ok <- FALSE
        break
      }
    }

    if (!data_ok)
      stop("The data field should only contain streamline objects.")
  }

  res <- list(
    name = input$name,
    case = input$case,
    scan = input$scan,
    side = input$side,
    data = input$data
  )

  class(res) <- c("tract", class(res))

  res
}

#' @export
#' @rdname as_tract
as_tract.mfData <- function(input, metadata, ...) {
  s <- seq(input$t0, input$tP, length.out = input$P)
  metadata$data <- seq_len(input$N) %>%
    purrr::map(~ streamline(
      s = s,
      x = as.numeric(input$fDList[[1]]$values[., ]),
      y = as.numeric(input$fDList[[2]]$values[., ]),
      z = as.numeric(input$fDList[[3]]$values[., ])
    ))
  as_tract(metadata)
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

  n <- purrr::map_int(filename, length)
  if (n < 3L || n > 4L)
    stop("The CSV file name should contain 3 or 4 fields separated by '_'")

  tract_name <- purrr::map_chr(filename, 3L)
  case <- purrr::map_chr(filename, 1L)
  scan <- purrr::map_chr(filename, 2L)
  side <- NA_character_
  if (n == 4L)
    side <- purrr::map_chr(filename, 4L)

  data <- readr::read_csv(path)

  if (nrow(data) == 0)
    return(tract(
      name = tract_name, case = case, scan = scan, side = side, data = list()
    ))

  required_vars <- c("LineID", "LPointID", "X", "Y", "Z")
  if (!all(required_vars %in% names(data)))
    stop("The data tibble should contain at least the variables LineID,
LPointID, X, Y, and Z.")

  model <- "None"
  check_dti <- data %>%
    names() %>%
    stringr::str_locate("Tensors#") %>%
    tibble::as_tibble() %>%
    tidyr::drop_na() %>%
    nrow()
  if (check_dti == 6L)
    model <- "DTI"

  data <- data %>%
    dplyr::arrange(LineID, LPointID) %>%
    add_diffusion_information(model) %>%
    dplyr::group_by(LineID) %>%
    dplyr::do(streamlines = streamline(
      x = .$X,
      y = .$Y,
      z = .$Z,
      diffusion = switch(model,
                         None = NA_real_,
                         DTI = .$Tensors)
    )) %>%
    dplyr::ungroup()

  tract(name = tract_name, case = case, scan = scan, side = side,
        data = data$streamlines)
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

  path <-
    stringr::str_c(tract$case, tract$scan, tract$name, tract$side, sep = "_") %>%
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
plot_tract <- function(tract, transparency = 1) {
  if (is_tract(tract))
    tract <- tibble::as_tibble(tract)
  else if (tibble::is_tibble(tract)) {
    if (!all(c("name", "case", "scan", "side", "data") %in% names(tract)))
      stop("The input tibble does not contain all the required fields.")
  } else
    stop("Input object is neither of class tract nor of class tibble.")

  name <- unique(tract$name)
  if (length(name) > 1L)
    stop("It is not possible to simulaneously plot different tracts.")

  case <- unique(tract$case)
  if (length(case) > 1L)
    stop("It is not possible to simulaneously plot a tract for different subjects.")

  print(
    tract %>%
      tidyr::unnest(.id = "LineID") %>%
      tidyr::gather(CoordName, CoordVal, x, y, z) %>%
      ggplot2::ggplot(ggplot2::aes(
        x = s,
        y = CoordVal,
        col = scan,
        fill = LineID
      )) +
      ggplot2::geom_line(size = 0.5, alpha = transparency) +
      ggplot2::facet_wrap(~ CoordName, scales = "free", nrow = 3L) +
      ggplot2::theme_bw() +
      ggplot2::xlab("Arc Length (mm)") +
      ggplot2::ylab("Coordinate (mm)") +
      ggplot2::guides(fill = FALSE) +
      ggplot2::theme(legend.position = "top") +
      ggplot2::ggtitle(paste(name, "of Subject", case))
  )

  invisible(tract)
}

#' Tract Combination
#'
#' @param ... A set of \code{\link{tract}}s.
#'
#' @return A \code{\link[tibble]{tibble}} with the following 5 variables:
#'   \code{name} with the name of the tract, \code{case} with the identifier of
#'   the subject to which the tract belongs, \code{scan} with the identifier of
#'   the repeated measure for a given subject, \code{side} with the hemisphere
#'   to which the tract belongs and \code{data} with the list of
#'   \code{streamline}s composing the tract.
#' @export
#'
#' @examples
#' file <- system.file("extdata", "Case001_CST_Left.csv", package = "fdatractography")
#' cst_left <- read_tract(file)
#' file <- system.file("extdata", "Case001_CST_Right.csv", package = "fdatractography")
#' cst_right <- read_tract(file)
#' bind_tracts(cst_left, cst_right)
bind_tracts <- function(...) {
  tract_list <- rlang::flatten(rlang::dots_values(...))
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
#' @param grid Uniform grid for the curvilinear abscissa (default: \code{0L}
#'   which uses the average number of points across streamlines defining the
#'   tract). Can be specificed either as an integer in which case the abscissa
#'   range of each streamline is used and abscissa is uniformly resampled within
#'   this range or as a numeric vector that will be taken as the new common
#'   abscissa for all streamlines.
#' @param validate A boolean specifying whether to check that first input is
#'   indeed a \code{\link{tract}}. should be checked (default: \code{TRUE}).
#'
#' @return A \code{\link{tract}} reparametrized according to the input grid.
#' @export
#'
#' @examples
#' file <- system.file("extdata", "Case001_CST_Left.csv", package = "fdatractography")
#' cst_left <- read_tract(file)
#' reparametrize_tract(cst_left)
reparametrize_tract <- function(tract, grid = 0L, validate = TRUE) {
  if (validate) {
    if (!is_tract(tract))
      stop("First input should be a tract object.")
  }

  UseMethod("reparametrize_tract", grid)
}

#' @export
#' @rdname reparametrize_tract
reparametrize_tract.integer <- function(tract, grid = 0L, validate = TRUE) {
  if (grid == 0L) {
    grid <- tract$data %>%
      purrr::map_int(nrow) %>%
      mean() %>%
      round()
  }

  tract$data <- tract$data %>%
    purrr::map(dplyr::do, streamline(
      s = modelr::seq_range(.$s, n = grid),
      x = approx(.$s, .$x, xout = s)$y,
      y = approx(.$s, .$y, xout = s)$y,
      z = approx(.$s, .$z, xout = s)$y,
      validate = FALSE
    ))

  tract
}

#' @export
#' @rdname reparametrize_tract
reparametrize_tract.numeric <- function(tract, grid = numeric(), validate = TRUE) {
  if (length(grid) == 0L) {
    grid_length <- tract$data %>%
      purrr::map_int(nrow) %>%
      mean() %>%
      round()
    max_abs <- tract$data %>%
      purrr::map_dbl(get_curvilinear_length) %>%
      min()
    grid <- seq(0, max_abs, length.out = grid_length)
  }

  tract$data <- tract$data %>%
    purrr::map(dplyr::do, streamline(
      s = grid,
      x = approx(.$s, .$x, xout = s)$y,
      y = approx(.$s, .$y, xout = s)$y,
      z = approx(.$s, .$z, xout = s)$y,
      validate = FALSE
    ))

  tract
}

#' Tract Hausdorff Distance
#'
#' @param tract A \code{\link{tract}}.
#' @param grid_length An integer specifying the grid length for uniform resampling of the 3D coordinates (default: \code{50L}).
#' @param ncores An integer specifying the number of cores available for the computations (default: \code{1L}).
#' @param nobs An integer specifying how many random \code{\link{streamline}}s should be kept from the input tract (default: \code{10L}).
#'
#' @return A vector of size \code{nobs * (nobs - 1) / 2} storing optimally the matrix of Hausdorff distances between \code{\link{streamline}}s that were kept as part of the input \code{\link{tract}}.
#' @export
#'
#' @examples
#' file <- system.file("extdata", "Case001_CST_Left.csv", package = "fdatractography")
#' cst_left <- read_tract(file)
#' get_distance_vector(cst_left)
get_distance_vector <- function(tract, grid_length = 50L, ncores = 1L, nobs = 10L, distance_fun = get_L2_distance) {
  tmp <- tract %>%
    reparametrize_tract(grid = grid_length, validate = FALSE) %>%
    tibble::as_tibble() %>%
    dplyr::slice(seq_len(nobs))

  parallel <- (ncores > 1L && requireNamespace("multidplyr", quietly = TRUE))
  if (parallel)
    cl <- multidplyr::create_cluster(cores = ncores) %>%
      multidplyr::cluster_copy(distance_fun)

  tmp <- tidyr::crossing(
    tmp %>% dplyr::select(str1 = data) %>% dplyr::mutate(id1 = seq_len(n())),
    tmp %>% dplyr::select(str2 = data) %>% dplyr::mutate(id2 = seq_len(n()))
  ) %>%
    dplyr::filter(id2 > id1)

  if (parallel) {
    tmp <- tmp %>%
      multidplyr::partition(cluster = cl) %>%
      dplyr::mutate(
        distances = purrr::map2_dbl(str1, str2, distance_fun)
      ) %>%
      dplyr::collect() %>%
      dplyr::ungroup() %>%
      dplyr::arrange(id1, id2)
  } else {
    tmp <- tmp %>%
      dplyr::mutate(
        distances = purrr::map2_dbl(str1, str2, distance_fun)
      ) %>%
      dplyr::arrange(id1, id2)
  }

  tmp$distances
}

align_tract <- function(tract) {
  min_index <- which.min(purrr::map_dbl(tract$data, get_curvilinear_length))
  min_streamline <- tract$data[[min_index]]
  tract$data <- tract$data %>%
    purrr::map(align_streamline,
               fixed_streamline = min_streamline,
               cost_function = cost_L1)
  tract
}

# align_tract <- function(tract) {
#   depths <- tract %>%
#     roahd::as.mfData() %>%
#     roahd::multiMBD()
#   median_index <- which.max(depths)
#   median_depth <- max(depths)
#
#   return(
#     list(
#       tract = tract,
#       median_idx = median_index,
#       max_depth = median_depth
#     )
#   )
#
#   current_index <- median_index + 1
#   pos <- 0
#   max_iter <- 10L
#   while (current_index != median_index && pos < max_iter) {
#     pos <- pos + 1
#     writeLines(paste0("    Obs. ", pos))
#     reference_streamline <- tract$data[[median_index]]
#     tract$data <- tract$data %>%
#       purrr::map(
#         align_streamline,
#         fixed_streamline = reference_streamline,
#         cost_function = cost_L1
#       )
#     tmp_depths <- tract %>%
#       roahd::as.mfData() %>%
#       roahd::multiMBD()
#     current_index <- median_index
#     median_index <- which.max(tmp_depths)
#     median_depth <- max(tmp_depths)
#   }
#
#   list(
#     tract = tract,
#     median_idx = median_index,
#     max_depth = median_depth
#   )
# }

#' Modified Band Depth for Tracts Objects
#'
#' This is a specialization of the \code{\link[roahd]{MBD}} function for
#' \code{\link{tract}} objects.
#'
#' @param tract A \code{\link{tract}} object.
#' @param validate A boolean that specifies whether the input format should be
#'   checked (default: \code{TRUE}).
#'
#' @return A numeric vector of depths for each \code{\link{streamline}} of the
#'   input \code{\link{tract}} object.
#' @export
#'
#' @examples
#' file <- system.file("extdata", "Case001_CST_Left.csv", package = "fdatractography")
#' cst_left <- read_tract(file)
#' MBD(cst_left)
MBD.tract <- function(tract) {
  tract %>%
    align_tract() %>%
    roahd::as.mfData() %>%
    roahd::multiMBD()
}

MBD_relative.tract <- function(tract_target, tract_reference, mfData_reference) {
  N_target <- length(tract_target$data)

  tract <- tract_target %>%
    bind_tracts(tract_reference) %>%
    align_tract() %>%
    roahd::as.mfData()

  relative_depths <- rep(0, tract$L)

  for (i in 1:tract$L) {
    fd_target <- tract$fDList[[i]][1:N_target, ]
    fd_reference <- tract$fDList[[i]][(N_target + 1):tract$N, ]
    relative_depths[i] <- roahd::MBD_relative(fd_target, fd_reference)
  }

  mean(relative_depths)
}

robust_clusterize <- function(tract, k = 1L, max_iter = 10L, nstart = 4L, ncores = 1L) {
  ncores <- min(ncores, nstart)
  parallel <- (ncores > 1L && requireNamespace("multidplyr", quietly = TRUE))

  if (parallel) {
    cl <- multidplyr::create_cluster(cores = ncores) %>%
      multidplyr::cluster_copy(clusterize) %>%
      multidplyr::cluster_copy(tract) %>%
      multidplyr::cluster_copy(k) %>%
      multidplyr::cluster_copy(max_iter)
    tmp <- tibble::tibble(id = seq_len(nstart)) %>%
      multidplyr::partition(cluster = cl) %>%
      dplyr::mutate(
        res = purrr::map(id, ~ clusterize(tract, k, max_iter))
      ) %>%
      dplyr::collect() %>%
      dplyr::ungroup() %>%
      dplyr::arrange(id)
    res <- tmp$res
  } else {
    res <- list()
    for (i in seq_len(nstart))
      res[[i]] <- clusterize(tract, k, max_iter)
  }
  idx <- which.min(purrr::map_dbl(res, "wmse"))
  res[[idx]]
}

clusterize <- function(tract, k = 1L, max_iter = 10L, validate = TRUE) {
  if (validate) {
    if (!is_tract(tract))
      stop("Input should be a tract object.")
  }

  n <- length(tract$data)
  labels <- seq_len(n)
  centroids <- sample.int(n, k)
  groups <- rep(0, n)
  reverse_groups <- list()
  for (i in seq_len(k)) {
    reverse_groups[[i]] <- centroids[i]
  }
  within <- rep(0, k)
  wmse <- rep(0, max_iter)

  clustered_tract <- list()

  for (iter in seq_len(max_iter)) {

    # Some verbose display
    str_iter <- formatC(
      iter,
      width = stringr::str_length(max_iter),
      format = "d",
      flag = "0"
    )
    writeLines(paste0("- Iteration ", str_iter, "/", max_iter))

    # # Assignment step: Phase 1 (corners)
    # # Handle corner cases to have 2 elements per cluster to start with,
    # # otherwise you cannot compute relative depths
    # corners <- centroids %>%
    #   purrr::map_int(function(idx) {
    #     distances <- rep(0, n)
    #     for (i in labels) {
    #       if (i %in% centroids)
    #         distances[i] <- 1e100
    #       else
    #         distances[i] <- get_L1_distance(tract$data[[i]], tract$data[[idx]])
    #     }
    #     which.min(distances)
    #   })
    # for (i in seq_len(k))
    #   reverse_groups[[i]] <- c(reverse_groups[[i]], corners[i])
    #
    # # Assignment step: Phase 2 (normal assignment)
    # for (i in labels) {
    #   writeLines(paste0("  * Assignment Step: Label ", i))
    #   skip_label <- FALSE
    #   for (j in seq_len(k)) {
    #     if (i %in% reverse_groups[[j]]) {
    #       skip_label <- TRUE
    #       break
    #     }
    #   }
    #   if (skip_label) next
    #   tract_target <- tract[i]
    #   relative_depths <- rep(0, k)
    #   for (j in seq_len(k)) {
    #     group_labels <- reverse_groups[[j]]
    #     tract_reference <- tract[group_labels]
    #     relative_depths[j] <- roahd::MBD_relative(tract_target, tract_reference)
    #   }
    #   optimal_cluster <- which.max(relative_depths)
    #   groups[i] <- optimal_cluster
    #   reverse_groups[[optimal_cluster]] <- unique(c(reverse_groups[[optimal_cluster]], i))
    # }

    # Assignment step
    # New attempt
    writeLines("  * Assign streamlines to clusters...")
    writeLines(paste0("    Centroids: ", centroids))
    if (iter == 1L) {
      # At first iteration, assign labels according to L1 distance to centroids
      groups <- labels %>%
        purrr::map_int(function(id) {
          distances <- rep(0, k)
          for (i in seq_len(k))
            distances[i] <- get_L1_distance(
              streamline1 = tract$data[[id]],
              streamline2 = tract$data[[centroids[i]]]
              )
          which.min(distances)
        })
    } else {
      # Assign labels to cluster with maximal relative depth
      groups <- labels %>%
        purrr::map_int(function(id) {
          relative_depths <- rep(0, k)
          for (i in seq_len(k)) {
            med <- roahd::median_mfData(clustered_tract[[i]])
            median_streamline <- streamline(
              s = seq(med$t0, med$tP, length.out = med$P),
              x = med$fDList[[1]]$values[1, ],
              y = med$fDList[[2]]$values[1, ],
              z = med$fDList[[3]]$values[1, ]
            )
            # current_streamline <- tract$data[[id]] %>%
            #   align_streamline(
            #     fixed_streamline = median_streamline,
            #     cost_function = cost_L1
            #   )
            current_streamline <- tract$data[[id]]
            current_tract <- tract[id]
            current_tract$data[[1]] <- current_streamline
            current_tract <- current_tract %>%
              reparametrize_tract(median_streamline$s, FALSE) %>%
              roahd::as.mfData()
            # TO DO: what happens if median is longer that current streamline?
            # Currently approx fill in with NA, how does roahd deal with it?
            relative_depths[i] <- 0
            for (j in seq_len(med$L)) {
              fd_target <- current_tract$fDList[[j]]
              fd_reference <- clustered_tract[[i]]$fDList[[j]]
              relative_depths[i] <- relative_depths[i] +
                roahd::MBD_relative(fd_target, fd_reference)
            }
          }
          which.max(relative_depths)
        })
    }

    for (i in seq_len(k))
      reverse_groups[[i]] <- labels[groups == i]

    # Find new cluster centers
    writeLines("  * Computing cluster medians...")
    for (i in seq_len(k)) {
      group <- (groups == i)
      group_labels <- labels[group]
      subtract <- tract[group_labels]
      aln <- align_tract(subtract)
      clustered_tract[[i]] <- roahd::as.mfData(aln$tract)
      centroids[i] <- group_labels[aln$median_idx]
      within[i] <- IQR(clustered_tract[[i]])# aln$max_depth
    }

    wmse[iter] <- mean(within)
    writeLines(paste0("    Sum of within-cluster maximum depths: ", wmse[iter]))
  }
  list(groups = groups, cendroids = centroids, within = within, wmse = wmse,
       clustered_tract = clustered_tract)
}

clusterize_test <- function(tract, k = 1L, max_iter = 10L) {
  n <- length(tract$data)
  labels <- seq_len(n)
  centroids <- sample.int(n, k)
  groups <- rep(0, n)
  wmse <- rep(0, max_iter)

  tract <- reparametrize_tract(tract, grid = numeric(), validate = FALSE)

  clustered_tract <- list()

  for (iter in seq_len(max_iter)) {

    # Some verbose display
    str_iter <- formatC(
      iter,
      width = stringr::str_length(max_iter),
      format = "d",
      flag = "0"
    )
    writeLines(paste0("- Iteration ", str_iter, "/", max_iter))

    reverse_groups <- list()
    for (i in seq_len(k)) {
      reverse_groups[[i]] <- centroids[i]
    }

    # Assignment step: Phase 1 (corners)
    # Handle corner cases to have 2 elements per cluster to start with,
    # otherwise you cannot compute relative depths
    corners <- centroids %>%
      purrr::map_int(function(idx) {
        distances <- rep(0, n)
        for (i in labels) {
          if (i %in% centroids)
            distances[i] <- 1e100
          else
            distances[i] <- get_L1_distance(tract$data[[i]], tract$data[[idx]])
        }
        which.min(distances)
      })
    for (i in seq_len(k))
      reverse_groups[[i]] <- c(reverse_groups[[i]], corners[i])

    # print(reverse_groups)

    # Assignment step: Phase 2 (normal assignment)
    for (i in labels) {
      writeLines(paste0("  * Assignment Step: Label ", i))
      skip_label <- FALSE
      for (j in seq_len(k)) {
        if (i %in% reverse_groups[[j]]) {
          skip_label <- TRUE
          break
        }
      }
      if (skip_label) next
      tract_target <- roahd::as.mfData(tract[i])
      # print(tract_target)
      relative_depths <- rep(0, k)
      for (j in seq_len(k)) {
        group_labels <- reverse_groups[[j]]
        tract_reference <- roahd::as.mfData(tract[group_labels])
        depth <- 0
        for (l in 1:tract_reference$L)
          depth <- depth + roahd::MBD_relative(tract_target$fDList[[l]], tract_reference$fDList[[l]])
        relative_depths[j] <- depth / tract_reference$L
      }
      optimal_cluster <- which.max(relative_depths)
      # depth_total <- depth_total + max(relative_depths)
      groups[i] <- optimal_cluster
      reverse_groups[[optimal_cluster]] <- unique(c(reverse_groups[[optimal_cluster]], i))
    }

    # Find new cluster centers
    writeLines("  * Computing cluster medians...")
    depth_total <- 0
    for (i in seq_len(k)) {
      group_labels <- reverse_groups[[i]]
      subtract <- tract[group_labels]
      clustered_tract[[i]] <- roahd::as.mfData(subtract)
      depths <- roahd::multiMBD(clustered_tract[[i]])
      centroids[i] <- group_labels[which.max(depths)]
      depth_total <- max(depths)
    }

    wmse[iter] <- depth_total
    writeLines(paste0("    Sum of within-cluster maximum depths: ", wmse[iter]))
  }
  list(groups = groups, cendroids = centroids, wmse = wmse,
       clustered_tract = clustered_tract)
}

simplified_clusterize <- function(tract, k = 1L, max_iter = 10L, validate = TRUE) {
  if (validate) {
    if (!is_tract(tract))
      stop("Input should be a tract object.")
  }

  n <- length(tract$data)
  labels <- seq_len(n)
  centroids <- sample.int(n, k)
  groups <- rep(0, n)
  reverse_groups <- list()
  for (i in seq_len(k)) {
    reverse_groups[[i]] <- centroids[i]
  }
  within <- rep(0, k)

  clustered_tract <- list()
  mf_tract <- roahd::as.mfData(tract)
  res <- list()

  for (iter in seq_len(max_iter)) {

    # Some verbose display
    str_iter <- formatC(
      iter,
      width = stringr::str_length(max_iter),
      format = "d",
      flag = "0"
    )
    writeLines(paste0("- Iteration ", str_iter, "/", max_iter))

    # Assignment step
    writeLines("  * Assign streamlines to clusters...")
    writeLines(paste0("    Centroids: ", centroids))
    total <- 0
    if (iter == 1L) {
      # At first iteration, assign labels according to L1 distance to centroids
      groups <- labels %>%
        purrr::map_int(function(id) {
          distances <- rep(0, k)
          for (i in seq_len(k))
            distances[i] <- get_L1_distance(
              streamline1 = tract$data[[id]],
              streamline2 = tract$data[[centroids[i]]]
            )
          which.min(distances)
        })
    } else {
      # Assign labels to cluster with maximal relative depth
      tmp <- labels %>%
        purrr::map(function(id) {
          relative_depths <- rep(0, k)
          for (i in seq_len(k)) {
            relative_depths[i] <- 0
            for (j in seq_len(mf_tract$L)) {
              fd_target <- mf_tract$fDList[[j]][id, ]
              fd_reference <- clustered_tract[[i]]$fDList[[j]]
              relative_depths[i] <- relative_depths[i] +
                roahd::MBD_relative(fd_target, fd_reference)
            }
          }
          list(maxi = max(relative_depths), argmaxi = which.max(relative_depths))
        })
      groups <- purrr::map_int(tmp, "argmaxi")
      total <- sum(purrr::map_dbl(tmp, "maxi"))
    }
    writeLines(paste0("Crit: ", total))

    for (i in seq_len(k))
      reverse_groups[[i]] <- labels[groups == i]

    # Find new cluster centers
    writeLines("  * Computing cluster medians...")
    for (i in seq_len(k)) {
      group <- (groups == i)
      group_labels <- labels[group]
      subtract <- mf_tract
      subtract$N <- length(group_labels)
      for (j in 1:mf_tract$L)
        subtract$fDList[[j]] <- mf_tract$fDList[[j]][group_labels, ]
      clustered_tract[[i]] <- subtract
      depths <- roahd::multiMBD(subtract)
      centroids[i] <- group_labels[which.max(depths)]
      within[i] <- IQR(subtract)# max(depths)
    }

    wmse <- mean(within)
    writeLines(paste0("    Sum of within-cluster maximum depths: ", wmse))

    res[[iter]] <- list(groups = groups, cendroids = centroids, within = within,
                        wmse = wmse, clustered_tract = clustered_tract)
  }
  res
}

find_best_cluster <- function(idx, k, groups, labels, depths, tract) {
  maximum_depths <- rep(0, k)
  for (j in seq_len(k)) {
    group <- (groups == j)
    group_labels <- labels[group]
    if (idx %in% group_labels) {
      maximum_depths[j] <- max(depths[group_labels])
      next
    }
    group_labels <- c(idx, group_labels)
    subtract <- tract
    subtract$data <- subtract$data[group_labels]
    maximum_depths[j] <- max(MBD(subtract))
  }
  which.max(maximum_depths)
}

as.mfData.tract <- function(x, grid_length = NULL, ...) {
  size <- min(purrr::map_dbl(x$data, get_curvilinear_length))

  if (is.null(grid_length))
    grid_length <- round(mean(purrr::map_int(x$data, nrow)))

  s <- seq(0, size, length.out = grid_length)

  tmp <- x %>%
    reparametrize_tract(grid = s, validate = FALSE) %>%
    tibble::as_tibble() %>%
    dplyr::select(data) %>%
    tidyr::unnest(.id = "LineID") %>%
    dplyr::group_by(LineID) %>%
    dplyr::mutate(PointID = seq_len(n())) %>%
    dplyr::ungroup(LineID)

  xmatrix <- tmp %>%
    dplyr::select(LineID, PointID, x) %>%
    tidyr::spread(PointID, x) %>%
    dplyr::select(-LineID) %>%
    as.matrix()

  ymatrix <- tmp %>%
    dplyr::select(LineID, PointID, y) %>%
    tidyr::spread(PointID, y) %>%
    dplyr::select(-LineID) %>%
    as.matrix()

  zmatrix <- tmp %>%
    dplyr::select(LineID, PointID, z) %>%
    tidyr::spread(PointID, z) %>%
    dplyr::select(-LineID) %>%
    as.matrix()

  roahd::mfData(s, list(xmatrix, ymatrix, zmatrix))
}

#' Inter-Quartile Range for arbitrary datasets
#'
#' This is a generic function to compute the inter-quartile range of a dataset.
#'
#' @param x Either a numeric vector or a \code{\link{tract}}.
#' @param na.rm A boolean for removing missing values (default: \code{FALSE}).
#'   For numeric vector inputs only.
#' @param type An integer selecting one of the many quantile algorithms, see
#'   \code{\link[stats]{quantile}}. For numeric vector inputs only.
#' @param validate A boolean that specifies whether the input format should be
#'   checked (default: \code{TRUE}). For \code{\link{tract}} inputs only.
#'
#' @return A positive scalar representing the inter-quartile range of the input
#'   dataset.
#' @export
#'
#' @examples
#' IQR(1:10)
#' file <- system.file("extdata", "Case001_CST_Left.csv", package = "fdatractography")
#' cst_left <- read_tract(file)
#' IQR(cst_left)
IQR <- function(x, ...) {
  UseMethod("IQR", x)
}

#' @export
#' @rdname IQR
IQR.numeric <- stats::IQR

#' @export
#' @rdname IQR
IQR.mfData <- function(x, ...) {
  depths <- roahd::multiMBD(x)
  depth_threshold <- max(depths) *
    optimize(depth_cost, c(0, 1), depth_data = depths)$minimum
  indices <- (depths >= depth_threshold)
  xmat <- x$fDList[[1]]$values[indices, ]
  ymat <- x$fDList[[2]]$values[indices, ]
  zmat <- x$fDList[[3]]$values[indices, ]
  volume <- 0
  for (i in 1:x$P) {
    points <- cbind(xmat[, i], ymat[, i], zmat[, i])
    area <- geometry::convhulln(points, "FA")$area
    volume <- volume + area
  }
  volume / (x$P - 1)
}

#' Subsetting operator for \code{\link{tract}} objects
#'
#' This method provides an easy and natural way to subset a tract
#' stored in a \code{\link{tract}} object, without having to deal with the inner
#' representation of the \code{\link{tract}} class.
#'
#' @param tract The input \code{\link{tract}}.
#' @param i A valid expression to extract subtract with fewer streamlines (could be an integer, a numeric vector, of a logical vector).
#'
#' @return A \code{\link{tract}} with the selected streamlines only.
#'
#' @examples
#' file <- system.file("extdata", "Case001_CST_Left.csv", package = "fdatractography")
#' cst_left <- read_tract(file)
#' t1 <- cst_left[1]
#' t2 <- cst_left[1:2]
#' n <- length(cst_left$data)
#' selected_streamlines <- sample(c(TRUE, FALSE), n, replace = TRUE)
#' t3 <- cst_left[selected_streamlines]
#' @export
"[.tract" <- function(tract, i) {
  if (missing(i))
    return(tract)
  tract$data <- tract$data[i]
  tract
}

#' Tract Simplification by Cluster Medians
#'
#' @param tract An object of class \code{\link{tract}}.
#' @param clusters A numeric or integer vector giving the cluster membership of
#'   each \code{\link{streamline}} in the \code{\link{tract}}.
#'
#' @return A \code{\link{tract}}.
#' @export
#'
#' @examples
#' file <- system.file("extdata", "Case001_CST_Left.csv", package = "fdatractography")
#' cst_left <- read_tract(file)
#' m <- simplify_tract(cst_left)
simplify_tract <- function(tract, clusters = NULL, validate = TRUE) {
  if (validate) {
    if (!is_tract(tract))
      stop("First argument should be an object of class tract.")
  }
  if (is.null(clusters) && tract$name == "CST") {
    cl <- cluster_cst(tract, validate = FALSE)
    writeLines(paste0("   --> ", cl$G, " clusters detected."))
    clusters <- cl$classification
  }
  all_clusters <- sort(unique(clusters))
  subtracts <- list()
  for (i in seq_along(all_clusters)) {
    cl <- all_clusters[i]
    writeLines(paste0(" - Simplyfing cluster #", cl, "..."))
    subtract <- tract[clusters == cl]
    if (length(subtract$data) == 1L) {
      subtracts[[i]] <- subtract
      next
    }
    subtract <- subtract %>%
      align_tract() %>%
      roahd::as.mfData() %>%
      roahd::median_mfData() %>%
      as_tract(metadata = tract)
    subtract$scan <- cl
    subtracts[[i]] <- subtract
  }
  bind_tracts(subtracts)
}


#' Clustering Method for the Cortico-Spinal Tract
#'
#' This function uses the \code{\link[mclust]{Mclust}} method from the
#' \code{mclust} package to fit a Gaussian mixture model to the end-points of
#' the CST on the precentral gyrus via Expectation-Maximization for selecting
#' the optimal clustering according to the Bayesian Information Criterion.
#'
#' @param cst An object of class \code{\link{tract}} representing the
#'   cortico-spinal tract.
#' @param validate A boolean specifying whether the input object should be
#'   checked or not (default: \code{TRUE}). For this function, it has to be of
#'   class \code{\link{tract}} and represent the cortico-spinal tract
#'   (\code{name} field of the \code{\link{tract}} should read \code{"CST"}).
#'
#' @return An object of class \code{\link[mclust]{Mclust}} providing the
#'   BIC-optimal Gaussian mixture model for clustering the CST streamlines based
#'   on cortical position of their end-point.
#' @export
#' @importFrom mclust mclustBIC
#'
#' @examples
#' file <- system.file("extdata", "Case001_CST_Left.csv", package = "fdatractography")
#' cst_left <- read_tract(file)
#' cl <- cluster_cst(cst_left)
cluster_cst <- function(cst, validate = TRUE) {
  if (validate) {
    if (!is_tract(cst))
      stop("The input should be an object of class tract.")
    if (cst$name != "CST")
      stop("The input tract should be a CST (name field should indicate CST).")
  }
  writeLines(" - Performing streamline clustering based on cortical position...")
  cst %>%
    tibble::as_tibble() %>%
    dplyr::transmute(
      x0 = purrr::map_dbl(data, ~ .$x[1]) %>% abs(),
      y0 = purrr::map_dbl(data, ~ .$y[1]),
      z0 = purrr::map_dbl(data, ~ .$z[1])
    ) %>%
    mclust::Mclust(G = seq_len(60L), prior = mclust::priorControl())
}
