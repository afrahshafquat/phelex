% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/utils.R
\name{estimate_misclassification_probability}
\alias{estimate_misclassification_probability}
\title{Calculates misclassification probability for samples}
\usage{
estimate_misclassification_probability(misclassified.samples,
  standardize = TRUE)
}
\arguments{
\item{misclassified.samples}{Matrix of misclassification indicators (misclassified.samples or misclassified.controls)}

\item{standardize}{Flag to normalize probability values (only when they are not normalized). Recommended to keep default.}
}
\value{
Misclassification probabilities for samples
}
\description{
Function estimates the probability
that a case is not a true case or a control is not a true control with the
matrix of misclassified.cases (misclassified.samples or misclassified.controls)
}
\keyword{misclassification,probability}
