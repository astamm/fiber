#' Tensor Invariants
#'
#' A set of functions calculating orthogonal tensor invariants as suggested in
#' the paper by Ennis, D. B., & Kindlmann, G. (2006). Orthogonal tensor
#' invariants and the analysis of diffusion tensor magnetic resonance images.
#' Magnetic resonance in medicine, 55(1), 136-146.
#'
#' \code{get_K1_invariant} computes the trace of the tensor,
#' \code{get_R2_invariant} computes the fractional anisotropy of the tensor and
#' \code{get_R3_invariant} computes the mode of the tensor.
#'
#' @param x A \code{3x3} symmetric definite positive matrix.
#' @param validate A boolean which is \code{TRUE} if the input type should be
#'   checked or \code{FALSE} otherwise (default: \code{TRUE}).
#'
#' @return A scalar giving the desired tensor invariant extracted from the input
#'   diffusion tensor.
#' @name get-invariant
#' @examples
#' DT <- diag(c(1.71e-3, 3e-4, 1e-4))
#' get_K1_invariant(DT)
#' get_R2_invariant(DT)
#' get_R3_invariant(DT)
NULL

#' @export
#' @keywords internal
get_tensor_invariants <- function(x, validate = TRUE) {
  if (validate) x <- validate_tensor(x)
  K1 <- sum(diag(x))
  MD <- K1 / 3
  Dtilde <- x - diag(MD, 3L)
  Dtilde_norm <- norm(Dtilde, "F")
  R2 <- sqrt(1.5) * Dtilde_norm / norm(x, "F")
  R3 <- 3 * sqrt(6) * det(Dtilde / Dtilde_norm)
  list(K1 = K1, R2 = R2, R3 = R3)
}

#' @rdname get-invariant
#' @export
get_K1_invariant <- function(x, validate = TRUE) {
  l <- get_tensor_invariants(x, validate)
  l$K1
}

#' @rdname get-invariant
#' @export
get_R2_invariant <- function(x, validate = TRUE) {
  l <- get_tensor_invariants(x, validate)
  l$R2
}

#' @rdname get-invariant
#' @export
get_R3_invariant <- function(x, validate = TRUE) {
  l <- get_tensor_invariants(x, validate)
  l$R3
}

#' @export
#' @keywords internal
get_LI_eigenvalues <- function (K1, R2, R3) {
  Coeff1 <- K1 / 3
  Coeff2 <- 2 * K1 * R2 / (3 * sqrt(3 - 2 * R2^2))
  lambda1 <- Coeff1 + Coeff2 * cos(acos(R3) / 3)
  lambda2 <- Coeff1 + Coeff2 * cos((acos(R3) - 2*pi) / 3)
  lambda3 <- Coeff1 + Coeff2 * cos((acos(R3) + 2*pi) / 3)
  c(lambda1, lambda2, lambda3)
}

#' Tensor Interpolation
#'
#' The function \code{approx_tensors} is the analog to the
#' \code{\link[stats]{approx}} function for \code{tensor}-valued data. It
#' performs tensor interpolation according to the paper by Gahm, J. K.,
#' Wisniewski, N., Kindlmann, G., Kung, G. L., Klug, W. S., Garfinkel, A., &
#' Ennis, D. B. (2012, October). Linear invariant tensor interpolation applied
#' to cardiac diffusion tensor MRI. In International Conference on Medical Image
#' Computing and Computer-Assisted Intervention (pp. 494-501). Springer, Berlin,
#' Heidelberg.
#'
#' @param x Numeric vector specifying the locations at which tensors are
#'   defined.
#' @param tensors List of \code{tensor} objects to be interpolated.
#' @param xout An optional set of numeric values specifying where interpolation
#'   is to take place.
#' @param method Specifies the interpolation method to be used. Choices are
#'   "linear" [default] or "constant".
#' @param n If \code{xout} is not specified, interpolation takes place at
#'   \code{n} equally spaced points spanning the interval \code{[min(x),
#'   max(x)]}. Default is \code{n = 50}.
#' @param yleft The value to be returned when input \code{x} values are less
#'   than \code{min(x)}. The default is defined by the value of \code{rule}
#'   given below.
#' @param yright The value to be returned when input \code{x} values are greater
#'   than \code{max(x)}. The default is defined by the value of \code{rule}
#'   given below.
#' @param rule An integer (of length 1 or 2) describing how interpolation is to
#'   take place outside the interval \code{[min(x), max(x)]}. If \code{rule} is
#'   1 then \code{NA}s are returned for such points and if it is 2, the value at
#'   the closest data extreme is used. Use, e.g., \code{rule = 2:1}, if the left
#'   and right side extrapolation should differ.
#' @param f For \code{method = "constant"}, a number between 0 and 1 inclusive,
#'   indicating a compromise between left- and right-continuous step functions.
#'   If \code{y0} and \code{y1} are the values to the left and right of the
#'   point then the value is \code{y0} if \code{f == 0}, \code{y1} if \code{f ==
#'   1}, and \code{y0*(1-f)+y1*f} for intermediate values. In this way, the
#'   result is right-continuous for \code{f == 0} and left-continuous for
#'   \code{f == 1}, even for non-finite \code{y} values.
#' @param ties Handling of tied \code{x} values. Either a function with a single
#'   vector argument returning a single number result or the string
#'   \code{"ordered"}.
#'
#' @return A \code{list} with components \code{x} and \code{y}, containing
#'   \code{n} coordinates which interpolate the given data points according to
#'   the \code{method} (and \code{rule}) desired.
#' @export
#'
#' @examples
#' theta <- c(0, pi/6, pi/4, pi/3, pi/2)
#' R <- NULL
#' for (i in seq_along(theta)) {
#'   th <- theta[i]
#'   R[[i]] <- cbind(
#'     c(cos(th), sin(th), 0),
#'     c(-sin(th), cos(th), 0),
#'     c(0, 0, 1)
#'   )
#' }
#' L <- diag(c(1.7, 0.3, 0.1))
#' D <- purrr::map(R, ~ . %*% L %*% t(.))
#' s <- seq(0, 1, length.out = length(theta))
#' approx_tensors(s, D)
approx_tensors <- function(x, tensors, xout, method = "linear", n = 50, yleft,
                           yright, rule = 1, f = 0, ties = mean) {
  out <- NULL
  for (i in 1:3) {
    for (j in 1:i) {
      y <- purrr::map_dbl(tensors, `[`, i, j)
      out[[length(out) + 1]] <- approx(x, y, xout, method, n, yleft, yright, rule, f, ties)$y
    }
  }

  R <- out %>%
    purrr::transpose() %>%
    purrr::map(purrr::flatten_dbl) %>%
    purrr::map(as_tensor, validate = FALSE) %>%
    purrr::map(~ eigen(., symmetric = TRUE)$vectors)

  tmp <- purrr::map(tensors, get_tensor_invariants, validate = FALSE)
  K1 <- purrr::map_dbl(tmp, "K1")
  K1 <- approx(x, K1, xout, method, n, yleft, yright, rule, f, ties)$y
  R2 <- purrr::map_dbl(tmp, "R2")
  R2 <- approx(x, R2, xout, method, n, yleft, yright, rule, f, ties)$y
  R3 <- purrr::map_dbl(tmp, "R3")
  tmp <- approx(x, R3, xout, method, n, yleft, yright, rule, f, ties)
  R3 <- tmp$y
  L <- purrr::pmap(list(K1, R2, R3), get_LI_eigenvalues)

  list(
    x = tmp$x,
    y = purrr::map2(R, L, ~ .x %*% diag(.y) %*% t(.x)) %>%
      purrr::map(as_tensor, validate = FALSE)
  )
}
