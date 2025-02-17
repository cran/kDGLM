% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/methods.R
\name{coefficients.fitted_dlm}
\alias{coefficients.fitted_dlm}
\title{coefficients.fitted_dlm}
\usage{
\method{coefficients}{fitted_dlm}(object, ...)
}
\arguments{
\item{object}{A fitted_dlm object.}

\item{...}{Arguments passed to coef.}
}
\value{
A list containing:
\itemize{
   \item data data.frame: A table with the model evaluated at each observed time.
   \item theta.mean matrix: The mean of the latent states at each time. Dimensions are n x t, where t is the size of t.eval and n is the number of latent states.
   \item theta.cov array: A 3D-array containing the covariance matrix of the latent states at each time. Dimensions are n x n x t, where t is the size of t.eval and n is the number of latent states.
   \item lambda.mean matrix: The mean of the linear predictor at each time. Dimensions are k x t, where t is the size of t.eval and k is the number of linear predictors.
   \item lambda.cov array: A 3D-array containing the covariance matrix for the linear predictor at each time. Dimensions are k x k x t, where t is the size of t.eval and k is the number of linear predictors.
   \item log.like, mae, mase, rae, mse, interval.score: The metric value at each time.
   \item conj.param list: A list containing, for each outcome, a data.frame with the parameter of the conjugated distribution at each time.
}
}
\description{
This method is wrapper for the coef method.
}
\seealso{
\code{\link{coef.fitted_dlm}}
}
\keyword{internal}
