% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/streamlines.R
\name{streamline}
\alias{streamline}
\title{Streamline Constructor}
\usage{
streamline(..., validate = TRUE)
}
\arguments{
\item{...}{A set of name-value pairs. Arguments are evaluated sequentially,
so you can refer to previously created variables. To be a valid tract, the
set should contain at least the fields x, y, z.}

\item{validate}{A boolean specifying whether the class of the input object
should be checked (default: \code{TRUE}).}
}
\value{
A \code{\link{streamline}}.
}
\description{
\code{streamline} is the constructor for objects of class \code{streamline}.
}
\examples{
st <- streamline(X = 1:10, Y = 1:10, Z = 1:10, Pie = 1:10)
}
