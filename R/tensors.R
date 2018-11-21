#' Tensor Constructor
#'
#' This function initialize a \code{\link{tensor}} object as a 3x3 matrix proportional to the identity.
#'
#' @param scale The proportionality constant.
#'
#' @return A diagonal \code{\link{tensor}} proportional to the identity.
#' @export
#'
#' @examples
#' tensor()
tensor <- function(scale = 1) {
  if (scale <= .Machine$double.eps)
    stop("The proportionality constant should be positive.")
  as_tensor(c(scale, 0, scale, 0, 0, scale), validate = FALSE)
}

#' Tensor Coercion
#'
#' This function transforms an input numeric vector or matrix into a proper
#' \code{\link{tensor}} object defined as a 3x3 definite positive symmetric
#' matrix.
#'
#' @param x A numeric vector of size 6 storing the 6 unique components of a
#'   tensor or a 3x3 numeric matrix.
#' @param validate A boolean which is \code{TRUE} if the input type should be
#'   checked or \code{FALSE} otherwise (default: \code{TRUE}).
#' @param twice A boolean that says whether off-diagonal tensor elements were
#'   doubled in vector representation (default: \code{FALSE}).
#'
#' @return A 3x3 symmetric postive definite matrix stored as a
#'   \code{\link{tensor}} object.
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
  if (validate) x <- .ValidateTensorMatrix(x)
  x <- .ConvertMatrixToTensor(x)
  class(x) <- c("tensor", class(x))
  x
}

#' @export
as_tensor.numeric <- function(x, validate = TRUE, ...) {
  if (validate) {
    x <- x %>%
      .ConvertTensorToMatrix() %>%
      .ValidateTensorMatrix() %>%
      .ConvertMatrixToTensor()
  }
  class(x) <- c("tensor", class(x))
  x
}

.ConvertTensorToMatrix <- function(x, twice = FALSE) {
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

  cbind(c(xx, yx, zx), c(yx, yy, zy), c(zx, zy, zz))
}

.ConvertMatrixToTensor <- function(x, twice = FALSE) {
  xx <- x[1,1]
  xy <- x[2,1]
  yy <- x[2,2]
  xz <- x[3,1]
  yz <- x[3,2]
  zz <- x[3,3]

  if (twice) {
    xy <- xy * 2
    xz <- xz * 2
    yz <- yz * 2
  }

  c(xx, xy, yy, xz, yz, zz)
}

.ValidateTensorMatrix <- function(x) {
  d <- nrow(x)
  if (d != 3L)
    stop("Input tensor should be of dimension 3.")
  if (ncol(x) != d)
    stop("Input tensor should be a square matrix.")
  if (!all(dplyr::near(x, t(x))))
    stop("Input tensor should be symmetric.")
  if (!all(eigen(x, symmetric = TRUE, only.values = TRUE)$values > .Machine$double.eps)) {
    warning("Input tensor is not SDP. Applying Matrix::nearPD.")
    x <- Matrix::nearPD(x)
  }
  x
}

#' Tensor Format Verification
#'
#' Check whether the input object is of class \code{\link{tensor}}.
#'
#' @param x An input R object.
#'
#' @return A boolean which is \code{TRUE} if the input object is of class
#'  \code{\link{tensor}} and \code{FALSE} otherwise.
#' @export
#'
#' @examples
#' is_tensor(tensor())
is_tensor <- function(x) {
  "tensor" %in% class(x)
}

#' @export
print.tensor <- function(x, ...) {
  .ConvertTensorToMatrix(x)
}

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
#' @param x An object of class \code{\link{tensor}}.
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
get_tensor_invariants <- function(x) {
  if (!is_tensor(x))
    stop("Input should be of class tensor.")
  x <- .ConvertTensorToMatrix(x)
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
get_K1_invariant <- function(x) {
  l <- get_tensor_invariants(x)
  l$K1
}

#' @rdname get-invariant
#' @export
get_R2_invariant <- function(x) {
  l <- get_tensor_invariants(x)
  l$R2
}

#' @rdname get-invariant
#' @export
get_R3_invariant <- function(x) {
  l <- get_tensor_invariants(x)
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
#' \code{\link[stats]{approx}} function for interpolating
#' \code{\link{tensor}}-valued data.
#'
#' If \code{method == "linear"}, it performs linear tensor interpolation
#' according to the paper by Gahm, J. K., Wisniewski, N., Kindlmann, G., Kung,
#' G. L., Klug, W. S., Garfinkel, A., & Ennis, D. B. (2012, October). Linear
#' invariant tensor interpolation applied to cardiac diffusion tensor MRI. In
#' International Conference on Medical Image Computing and Computer-Assisted
#' Intervention (pp. 494-501). Springer, Berlin, Heidelberg. If \code{method ==
#' "log"}, it performs log-Euclidean tensor interpolation according to the paper
#' by Arsigny, V., Fillard, P., Pennec, X. and Ayache, N. (2006) ‘Log-Euclidean
#' metrics for fast and simple calculus on diffusion tensors.’, Magnetic
#' resonance in medicine : official journal of the Society of Magnetic Resonance
#' in Medicine / Society of Magnetic Resonance in Medicine, 56(2), pp. 411–21.
#' doi: 10.1002/mrm.20965.
#'
#' @param x Numeric vector specifying the locations at which tensors are
#'   defined.
#' @param tensors List of \code{\link{tensor}} objects to be interpolated.
#' @param xout An optional set of numeric values specifying where interpolation
#'   is to take place.
#' @param method Specifies the interpolation method to be used. Choices are
#'   "linear" [default] or "log" for interpolating in the log-Euclidean space.
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
#' approx_tensors(s, D, method = "log")
approx_tensors <- function(x, tensors, xout, method = "linear", n = 50, yleft,
                           yright, rule = 1, ties = mean) {
  method <- match.arg(method, c("linear", "log"))
  if (is.na(method))
    stop("Invalid interpolation method. Choices are: linear or log.")

  stopifnot(is.numeric(rule), (lenR <- length(rule)) >= 1L, lenR <= 2L)

  if (lenR == 1)
    rule <- rule[c(1, 1)]

  if (missing(yleft))
    yleft <- if (rule[1L] == 1)
      NA
  else y[1L]

  if (missing(yright))
    yright <- if (rule[2L] == 1)
      NA
  else y[length(y)]

  if (method == "log")
    tensors <- purrr::map(tensors, tensor_log)

  out <- tensors %>%
    purrr::transpose() %>%
    purrr::simplify_all() %>%
    purrr::map(~ stats::approx(x, .x, xout, "linear", n, yleft, yright, rule, 0, ties))

  R <- out %>%
    purrr::map("y") %>%
    purrr::transpose() %>%
    purrr::simplify_all()

  if (method == "log") {
    return(list(
      x = out[[1]]$x,
      y = R %>%
        purrr::map(tensor_exp, validate = FALSE) %>%
        purrr::map(as_tensor, validate = FALSE)
    ))
  }

  R <- R %>%
    purrr::map(.ConvertTensorToMatrix) %>%
    purrr::map(eigen, symmetric = TRUE) %>%
    purrr::map("vectors")

  tmp <- purrr::map(tensors, get_tensor_invariants)
  K1 <- purrr::map_dbl(tmp, "K1")
  K1 <- stats::approx(x, K1, xout, "linear", n, yleft, yright, rule, 0, ties)$y
  R2 <- purrr::map_dbl(tmp, "R2")
  R2 <- stats::approx(x, R2, xout, "linear", n, yleft, yright, rule, 0, ties)$y
  R3 <- purrr::map_dbl(tmp, "R3")
  R3 <- stats::approx(x, R3, xout, "linear", n, yleft, yright, rule, 0, ties)$y
  L <- purrr::pmap(list(K1, R2, R3), get_LI_eigenvalues)

  list(
    x = out[[1]]$x,
    y = purrr::map2(R, L, ~ .x %*% diag(.y) %*% t(.x)) %>%
      purrr::map(.ConvertMatrixToTensor) %>%
      purrr::map(as_tensor, validate = FALSE)
  )
}

#' Tensor Exponentiation and Logarithm
#'
#' @param x An object of class \code{\link{tensor}}.
#' @param validate A boolean specifying wether the function should check that
#'   input is a tensor.
#'
#' @return A 3x3 symmetric matrix which is the exponential of the input in the
#'   case of \code{tensor_exp} or the logarithm of the input in the case of
#'   \code{tensor_log}. It is stored as a 6-dimensional vector.
#' @name tensor-logexp
#'
#' @examples
#' tensor_exp(tensor())
#' tensor_log(tensor())
NULL

#' @rdname tensor-logexp
#' @export
tensor_exp <- function(x, validate = TRUE) {
  if (validate) {
    if (!is_tensor(x))
      stop("Input should be of class tensor.")
  }
  x <- .ConvertTensorToMatrix(x)
  eig <- eigen(x, symmetric = TRUE)
  x <- eig$vectors %*% diag(exp(eig$values)) %*% t(eig$vectors)
  .ConvertMatrixToTensor(x)
}

#' @rdname tensor-logexp
#' @export
tensor_log <- function(x, validate = TRUE) {
  if (validate) {
    if (!is_tensor(x))
      stop("Input should be of class tensor.")
  }
  x <- .ConvertTensorToMatrix(x)
  eig <- eigen(x, symmetric = TRUE)
  x <- eig$vectors %*% diag(log(eig$values)) %*% t(eig$vectors)
  .ConvertMatrixToTensor(x)
}
