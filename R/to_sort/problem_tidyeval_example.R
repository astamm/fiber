library(tidyverse)

n <- 100L
df <- tibble(
  s = seq(0, 2*pi, len = n),
  x = cos(s),
  y = sin(s)
)

### This works but one variable is allowed only

resample <- function(df, var, n = 200L) {
  var <- enquo(var)
  var_name <- quo_name(var)
  df %>%
    do(tibble(
      s = modelr::seq_range(.$s, n),
      !!var_name := approx(.$s, .[[!!var_name]], n = n)$y
    ))
}

resample(df, x)

### With multiple variables to resample, I cannot make it work by splicing with !!!

resample_multi <- function(df, ..., n = 200L) {
  vars <- quos(...)
  var_names <- map_chr(vars, quo_name)
  l <- list()
  for (i in seq_along(var_names))
    l[[i]] <- quo(approx(df$s, df[[!!var_names[i]]], n = n)$y)
  names(l) <- var_names
  print(l)

  tibble(
      s = modelr::seq_range(df$s, n),
      !!!l
    )
}

resample_multi(df, x, y)

resample_multi2 <- function(df, ..., n = 200L) {
  vars <- quos(...)
  var_names <- map_chr(vars, quo_name)
  res <- df %>%
    do(tibble(
      s = modelr::seq_range(.$s, n),
      !!!(map(var_names, function(xx) quo(!!xx := approx(.$s, .[[!!xx]], n = n)$y)))
    ))
  res
}

resample_multi2(df, x, y)

### But it works by sequentially applying mutate once per variables requested in ...

resamp <- function(data, ..., n = 0L, range = NULL) {
  if (n == 0L) n <- nrow(data)
  if (is.null(range)) range <- data$s
  res <- tibble(s = modelr::seq_range(range, n = n))

  vars <- purrr::map_chr(dplyr::quos(...), dplyr::quo_name)
  for (var in vars) {
    res <- res %>%
      mutate(!!var := approx(data$s, data[[!!var]], xout = s)$y)
  }
  res
}

resamp(df, x, y, n = 200)

### Hadley's version

resample_hadley <- function(data, ..., n = 0L, range = NULL) {
  if (n == 0L) n <- nrow(data)
  if (is.null(range)) range <- data$s
  s_new <- modelr::seq_range(range, n = n)

  vars <- rlang::quos(...)
  new_vars <- purrr::map(vars, function(x) {
    rlang::expr(approx(!!data$s, dplyr::pull(!!data, !!x), xout = !!s_new)$y)
  })
  names(new_vars) <- purrr::map_chr(vars, rlang::quo_name)

  tibble::tibble(s_new = s_new, !!!new_vars)
}

resample_hadley(df, x, y, n = 200L, range = c(0, 0.06))
