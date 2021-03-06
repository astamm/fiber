% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/tensors.R
\name{tensor}
\alias{tensor}
\title{Tensor Constructor}
\usage{
tensor(scale = 1)
}
\arguments{
\item{scale}{The proportionality constant.}
}
\value{
A diagonal \code{\link{tensor}} proportional to the identity.
}
\description{
This function initialize a \code{\link{tensor}} object as a 3x3 matrix proportional to the identity.
}
\examples{
tensor()
}
