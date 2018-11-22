#' Diffusion Biomarkers
#'
#' \code{get_fractional_anisotropy} computes the fractional anisotropy of a
#' diffusion tensor.
#' \code{get_mean_diffusivity} computes the mean diffusivity of a diffusion tensor.
#'
#' @param x A tensor or a streamline.
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
  l <- get_microstructure(x)
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
  l <- get_microstructure(x)
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
