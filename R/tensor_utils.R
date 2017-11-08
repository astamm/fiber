#' Tensor Invariants
#'
#' \code{get_fractional_anisotropy} computes the fractional anisotropy of a
#' diffusion tensor.
#' \code{get_mean_diffusivity} computes the mean diffusivity of a diffusion tensor.
#'
#' @param x A \code{3x3} symmetric definite positive matrix.
#' @param validate A boolean which is \code{TRUE} if the input type should be
#'   checked or \code{FALSE} otherwise (default: \code{TRUE}).
#'
#' @return A scalar giving the desired tensor invariant extracted from the input diffusion tensor.
#' @name get-invariant
#' @examples
#' DT <- diag(c(1.71e-3, 3e-4, 1e-4))
#' get_K1_invariant(DT)
#' get_R2_invariant(DT)
#' get_R2_invariant(DT)
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
#' @param x
#' @param tensors
#' @param xout
#' @param method
#' @param n
#' @param yleft
#' @param yright
#' @param rule
#' @param f
#' @param ties
#'
#' @return
#' @export
#'
#' @examples
approx_tensors <- function(x, tensors, xout, method = "linear", n = 50, yleft,
                           yright, rule = 1, f = 0, ties = mean) {
  out <- NULL
  for (i in 1:3) {
    for (j in 1:i) {
      y <- purrr::map_dbl(tensors, `[`, i, j)
      out[[length(out) + 1]] <- approx(x, y, xout, method, n, yleft, yright, rule = 1, f = 0, ties)$y
    }
  }

  R <- out %>%
    purrr::transpose() %>%
    purrr::map(purrr::flatten_dbl) %>%
    purrr::map(as_tensor, validate = FALSE) %>%
    purrr::map(~ eigen(., symmetric = TRUE)$vectors)

  tmp <- purrr::map(tensors, get_tensor_invariants, validate = FALSE)
  K1 <- purrr::map_dbl(tmp, "K1")
  K1 <- approx(x, K1, xout, method, n, yleft, yright, rule = 1, f = 0, ties)$y
  R2 <- purrr::map_dbl(tmp, "R2")
  R2 <- approx(x, R2, xout, method, n, yleft, yright, rule = 1, f = 0, ties)$y
  R3 <- purrr::map_dbl(tmp, "R3")
  R3 <- approx(x, R3, xout, method, n, yleft, yright, rule = 1, f = 0, ties)$y
  L <- purrr::pmap(list(K1, R2, R3), get_LI_eigenvalues)

  purrr::map2(R, L, ~ .x %*% diag(.y) %*% t(.x)) %>%
    purrr::map(as_tensor, validate = FALSE)
}
