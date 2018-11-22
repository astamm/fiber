#' Diffusion Biomarkers
#'
#' \code{get_fractional_anisotropy} computes the fractional anisotropy of a
#' diffusion tensor.
#' \code{get_mean_diffusivity} computes the mean diffusivity of a diffusion tensor.
#'
#' @param x An object of class \code{\link{tensor}}.
#'
#' @return A scalar giving the desired diffusion biomarker extracted from the input diffusion tensor.
#' @name get-biomarker
#' @examples
#' DT <- as_tensor(diag(c(1.71e-3, 3e-4, 1e-4)))
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
  x <- .ConvertTensorToMatrix(x)
  vals <- eigen(x, symmetric = TRUE, only.values = TRUE)$values
  m <- mean(vals)
  list(
    AD = vals[1],
    RD = (vals[2] + vals[3]) / 2,
    MD = m,
    FA = sqrt(1.5 * sum((vals - m)^2) / sum(vals^2))
  )
}

#' @rdname get-biomarker
#' @export
get_axial_diffusivity <- function(x) {
  if (!is_tensor(x))
    stop("Axial diffusivity is a tensor-specific feature but input is not of class tensor.")
  l <- get_microstructure(x)
  l$AD
}

#' @rdname get-biomarker
#' @export
get_radial_diffusivity <- function(x) {
  if (!is_tensor(x))
    stop("Radial diffusivity is a tensor-specific feature but input is not of class tensor.")
  l <- get_microstructure(x)
  l$RD
}

#' @rdname get-biomarker
#' @export
get_mean_diffusivity <- function(x) {
  if (!is_tensor(x))
    stop("Mean diffusivity is a tensor-specific feature but input is not of class tensor.")
  l <- get_microstructure(x)
  l$MD
}

#' @rdname get-biomarker
#' @export
get_fractional_anisotropy <- function(x) {
  if (!is_tensor(x))
    stop("Fractional anisotropy is a tensor-specific feature but input is not of class tensor.")
  l <- get_microstructure(x)
  l$FA
}

#' @export
add_biomarkers <- function(streamline, model = "None") {
  switch(
    model,
    None = streamline,
    SingleTensor = streamline %>%
      dplyr::mutate(
        t = list(`tensors#0`, `tensors#1`, `tensors#2`, `tensors#3`,`tensors#4`, `tensors#5`) %>%
          purrr::pmap(c) %>%
          purrr::map(as_tensor)
      ) %>%
      dplyr::select(s, x, y, z, t) %>%
      as_streamline(validate = FALSE)
  )
}
