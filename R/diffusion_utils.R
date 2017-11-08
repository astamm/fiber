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
as_tensor <- function(x, validate = TRUE, ...) {
  UseMethod("as_tensor", x)
}

#' @export
as_tensor.matrix <- function(x, validate = TRUE, ...) {
  if (validate) x <- validate_tensor(x)
  class(x) <- c("tensor", class(x))
  x
}

#' @export
as_tensor.numeric <- function(x, validate = TRUE, ...) {
  .as_tensor_numeric(x, validate, ...)
}

.as_tensor_numeric <- function(x, validate = TRUE, twice = FALSE, ...) {
  if (length(x) != 6L)
    stop("Input vector should be of dimension 6.")

  xx <- x[1]
  yx <- x[2]
  yy <- x[3]
  zx <- x[4]
  zy <- x[5]
  zz <- x[6]

  if (twice) {
    yx <- yx / 2
    zx <- zx / 2
    zy <- zy / 2
  }

  output_matrix <- cbind(c(xx, yx, zx), c(yx, yy, zy), c(zx, zy, zz))
  as_tensor(output_matrix, validate)
}

is_tensor <- function(x) {
  "tensor" %in% class(x)
}

is.tensor <- is_tensor

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
get_microstructure <- function(x) {
  UseMethod("get_microstructure", x)
}

#' @export
get_microstructure.tensor <- function(x) {
  vals <- eigen(tensor, symmetric = TRUE, only.values = TRUE)$values
  m <- mean(vals)
  list(
    AD = vals[1],
    RD = (vals[2] + vals[3]) / 2,
    MD = m,
    FA = sqrt(1.5 * sum((vals - m)^2) / sum(vals^2))
  )
}

#' @export
get_microstructure.streamline <- function(x) {
  x %>%
    dplyr::mutate(
      microstructure = purrr::map(diffusion, get_microstructure)
    ) %>%
    as_streamline(validate = FALSE)
}

#' @rdname get-biomarker
#' @export
get_axial_diffusivity <- function(x) {
  UseMethod("get_axial_diffusivity", x)
}

#' @export
get_axial_diffusivity.tensor <- function(x) {
  l <- get_microstructure(x)
  l$AD
}

#' @export
get_axial_diffusivity.streamline <- function(x) {
  x %>%
    dplyr::mutate(
      AD = purrr::map_dbl(diffusion, get_axial_diffusivity)
    ) %>%
  as_streamline(validate = FALSE)
}

#' @rdname get-biomarker
#' @export
get_radial_diffusivity <- function(x) {
  UseMethod("get_radial_diffusivity", x)
}

#' @export
get_radial_diffusivity.tensor <- function(x) {
  l <- get_microstructure(x)
  l$RD
}

#' @export
get_radial_diffusivity.streamline <- function(x) {
  x %>%
    dplyr::mutate(
      RD = purrr::map_dbl(diffusion, get_radial_diffusivity)
    ) %>%
    as_streamline(validate = FALSE)
}

#' @rdname get-biomarker
#' @export
get_mean_diffusivity <- function(x) {
  UseMethod("get_mean_diffusivity", x)
}

#' @export
get_mean_diffusivity.tensor <- function(x) {
  l <- get_tensor_microstructure(x)
  l$MD
}

#' @export
get_mean_diffusivity.streamline <- function(x) {
  x %>%
    dplyr::mutate(
      MD = purrr::map_dbl(diffusion, get_mean_diffusivity)
    ) %>%
    as_streamline(validate = FALSE)
}

#' @rdname get-biomarker
#' @export
get_fractional_anisotropy <- function(x) {
  UseMethod("get_fractional_anisotropy", x)
}

#' @export
get_fractional_anisotropy.tensor <- function(x) {
  l <- get_tensor_microstructure(x)
  l$FA
}

#' @export
get_fractional_anisotropy.streamline <- function(x) {
  x %>%
    dplyr::mutate(
      FA = purrr::map_dbl(diffusion, get_fractional_anisotropy)
    ) %>%
    as_streamline(validate = FALSE)
}

add_biomarkers <- function(data, model = "None") {
  switch(
    model,
    None = data,
    DTI = data %>%
      dplyr::mutate(
        Tensors = purrr::pmap(list(`Tensors#0`, `Tensors#1`, `Tensors#2`,
                                   `Tensors#3`,`Tensors#4`, `Tensors#5`), c) %>%
          purrr::map(as_tensor) %>%
          purrr::map(get_tensor_microstructure),
        AD = purrr::map_dbl(Tensors, "AD"),
        RD = purrr::map_dbl(Tensors, "RD"),
        MD = purrr::map_dbl(Tensors, "MD"),
        FA = purrr::map_dbl(Tensors, "FA")
      ) %>%
      dplyr::select(-dplyr::starts_with("Tensors"))
  )
}

validate_tensor <- function(tensor) {
  d <- nrow(tensor)
  if (d != 3L)
    stop("Input tensor should be of dimension 3.")
  if (ncol(tensor) != d)
    stop("Input tensor should be a square matrix.")
  if (!all(dplyr::near(tensor, t(tensor))))
    stop("Input tensor should be symmetric.")
  if (!all(eigen(tensor, symmetric = TRUE, only.values = TRUE)$values > .Machine$double.eps)) {
    warning("Input tensor is not SDP. Applying Matrix::nearPD.")
    tensor <- Matrix::nearPD(tensor)
  }
  tensor
}
