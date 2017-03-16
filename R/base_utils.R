#' @importFrom dplyr %>%
#' @export
dplyr::`%>%`

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
