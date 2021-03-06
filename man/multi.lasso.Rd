% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/multi.lasso.R
\name{multi.lasso}
\alias{multi.lasso}
\title{Multi-split Lasso}
\usage{
multi.lasso(X, Y, p.fac = NULL, fold.num)
}
\arguments{
\item{X}{An input matrix, of dimension nobs x nvars.}

\item{Y}{A numeric response vector, containing nobs variables.}

\item{p.fac}{A sequence of penalty factor applied on each variable.}

\item{fold.num}{The number of cross validation folds.}
}
\value{
A list of two numeric objects of index of (1) selected and (2) unselected variables.
}
\description{
Multi-splitted variable selection using Lasso
}
