% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/phelex_mm.R
\name{phelex_mm}
\alias{phelex_mm}
\title{Estimate misclassification in phenotype}
\usage{
phelex_mm(x, y, iterations = 1e+05, alpha.prior = c(10, 1),
  lambda.prior = c(1, 1), link = "pnorm", beta.prior = "norm",
  beta.prior.params = c(0, 1), beta.initial.vec = NULL, mu.update = 0.5,
  verbose = TRUE, normalize = FALSE, stamp = 1000)
}
\arguments{
\item{x}{Genotype matrix with dimensions n x m.}

\item{y}{Phenotype vector with length n.}

\item{iterations}{Total number of iterations for Metropolis-Hastings Sampling}

\item{alpha.prior}{Beta prior parameters for true-positve rate}

\item{lambda.prior}{Beta prior parameters for false-positive rate.}

\item{link}{probit or logistic model (options pnorm and plogis respectively)}

\item{beta.prior}{'norm'(default): Normal prior or 'unif' Uniform prior on fixed effects}

\item{beta.prior.params}{if beta.prior is norm, then (mean, sd), if unif then (min, max)}

\item{beta.initial.vec}{Vector of initial values for beta parameters i.e. effect sizes for n snps.}

\item{mu.update}{Fractions of iterations to start updating mean fixed effects or mu. Default = 0.5}

\item{verbose}{Default TRUE. Prints progress information}

\item{normalize}{Default FALSE. Scales liability computed}

\item{stamp}{Iteration breakpoint to print time}
}
\value{
List containing
\itemize{
 \item betas: Matrix of estimated effect sizes for each SNP (SNPs[rows] x iterations[columns]).
 \item parameters: Matrix with estimated parameter values[rows] across iterations[columns].
 Order is c(mu, alpha, lambda) where alpha = true positive rate, lambda = false positive rate, mu = mean effect size value
 \item misclassified.cases: Matrix of alpha vectors where 1s represent false positives and 0s represent true positives as inferred at each iterations
 \item misclassified.controls: Matrix of lambda vectors where 1s represent false negatives and 0s represent true negatives as inferred at each iterations
 \item posterior: Vector of posterior probability across iterations
 \item accept: Vector of 1/0 values across iterations; 1 indicates proposal was accepted at iteration;0 o/w
 }
}
\description{
Estimates misclassification probabilities in observed GWAS phenotype y given genotypes dataset x.
The method follows the PheLEx-m algorithm to predict misclassification probabilities
using Adaptive Metropolis-Hastings defined by Shafquat et al.
}
\keyword{hastings}
\keyword{metropolis}
\keyword{phelex_mm,misclassification,GWAS,phenotype,adaptive}
