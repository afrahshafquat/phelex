#' Estimate misclassification in phenotype using Gibbs Sampler
#'
#' Estimates misclassification probabilities in observed GWAS phenotype y given genotypes dataset x.
#' The method follows the algorithm defined by Rekaya et. al (PMC5138056; 2016) to predict misclassification probabilities
#' using Gibbs Sampling algorithm.
#'
#' @param y Phenotype vector with length n.
#' @param x Genotype matrix with dimensions n x i.
#' @param pi1.prior hyperparameters for false positive rate pi1. Default Beta(1, 4).
#' @param pi2.prior hyperparameters for false negative rate pi2. Default Beta(1, 4).
#' @param beta.initial.vec Initial values for beta parameters in order c(Beta_a_1,...Beta_a_i). Default values are random values for all parameters.
#' @param iterations Number of iterations for sampling
#' @param stamp Iteration breakpoint to print time
#' @param verbose Default TRUE. Prints progress information
#'
#' @return  List containing
#' \itemize{
#'  \item betas: Matrix of estimated effect sizes for each SNP (SNPs[rows] x iterations[columns]).
#'  \item parameters: Matrix with estimated parameter values[rows] across iterations[columns].
#'  Order is c(mu, pi1, pi2) where mu is mean effect size value, pi1/pi2 false positive/false negative rates.
#'  \item misclassified.cases: Matrix of vectors where 1s represent false positives and 0s represent true positives as inferred at each iterations
#'  \item misclassified.controls: Matrix of vectors where 1s represent false negatives and 0s represent true negatives as inferred at each iterations
#'  }
#' @keywords keywords
#'
#' @import truncdist
#' @import stats
#' @export

rekaya = function(y,
                 x,
                 pi1.prior = c(1, 4),
                 pi2.prior = c(1, 4),
                 beta.initial.vec = NULL,
                 iterations = 1e4,
                 stamp = 1e3,
                 verbose = T) {

  if(verbose) print(paste('Entering function', date()))

  a1 = pi1.prior[1]; b1 = pi1.prior[2];
  a2 = pi2.prior[1]; b2 = pi2.prior[2];

  y = y;
  case.inds = which(y == 1)
  control.inds = which(y == 0)
  n1 = length(case.inds);
  n2 = length(control.inds);
  n = length(y)

  r = y
  lj = y

  x = x;
  markers = ncol(x)  # Number of betas to estimate

  iterations = iterations

  a1.tmp = a1;
  a2.tmp = a2;
  b1.tmp = b1;
  b2.tmp = b2;

  pi1 = rbeta(1, a1, b1);
  pi2 = rbeta(1, a2, b2);
  mu = rnorm(1, 0, .1)

  if (is.null(beta.initial.vec)) {  # Initialize parameter values
    beta.initial.vec = rnorm(markers, 0, 0.1)
  }

  beta = as.matrix(beta.initial.vec, nrow = markers, ncol = 1);

  alfa = matrix(0, nrow = n1)
  lambda = matrix(0, nrow = n2)

  num.params = 3 #c(mu, pi1, pi2)

  parameters = matrix(NA, nrow = num.params, ncol = iterations)
  betas = matrix(NA, nrow = markers, ncol = iterations)
  alfas = matrix(NA, nrow = n1, ncol = iterations)
  lambdas = matrix(NA, nrow = n2, ncol = iterations)

  parameters[, 1] = c(mu, pi1, pi2)
  betas[, 1] = beta[]
  alfas[, 1] = alfa
  lambdas[, 1] = lambda

  if(verbose) print(paste('Starting Gibbs', date()))

  for(iteri in 2:iterations) {
    if(!(iteri %% stamp) & verbose) {
      print(paste(iteri, date()))
    }

    beta.tmp = beta  # Update Betas
    for(j in 1:markers) {
      beta.tmp[j,] = 0
      x.tmp = x
      x.tmp[, j] = 0
      inv.x = ((t(x[, j]) %*% x[, j]) ^ -1);
      bj = inv.x * (t(x[, j]) %*% (lj - x.tmp %*% beta.tmp))
      beta.tmp[j, ] = rnorm(1, bj, sd = sqrt(inv.x));
    }
    beta = beta.tmp

    mu = sum(lj - x %*% beta) / n  # Update Mu
    mu = rnorm(1, mean = mu, sd = 1 / sqrt(n))

    l.mu = mu + x %*% beta  # Update liability
    lj = r
    lj[r==1] = rtrunc(sum(r), spec = 'norm', mean = l.mu[r==1], sd = 1, a = 0)
    lj[r==0] = rtrunc(sum(r==0), spec = 'norm', mean = l.mu[r==0], sd = 1, b = 0)

    pi.b = pnorm(lj)  # Update piB funciton of SNP effects

    a1.tmp = a1 + sum(alfa, na.rm = T)  # Update Priors for pi1, pi2
    b1.tmp = b1 + n1 - sum(alfa, na.rm = T)

    a2.tmp = a2 + sum(lambda, na.rm = T)
    b2.tmp = b2 + n2 - sum(lambda, na.rm=T)

    pi1 = rbeta(1, shape1 = a1.tmp, shape2 = b1.tmp)  # Update pi1, pi2
    pi2 = rbeta(1, shape1 = a2.tmp, shape2 = b2.tmp)

    p.a.1 = pi1 * (1 - pi.b[case.inds])  # Update Alpha
    p.a.0 = (1 - pi1) * (pi.b[case.inds])
    p.alfa = p.a.1 / (p.a.1 + p.a.0)
    alfa = rbinom(n1, 1, prob = p.alfa)

    p.l.1 = pi2 * (pi.b[control.inds])  # Update Lambda
    p.l.0 = (1 - pi2) * (1 - pi.b[control.inds])
    p.lambda = p.l.1 / (p.l.1 + p.l.0)
    lambda = rbinom(n2, 1, prob = p.lambda)

    flipped.cases = case.inds[(alfa == 1)]
    flipped.controls = control.inds[(lambda == 1)]

    r = y  # Update r (real phenotype)
    r[flipped.controls] = 1
    r[flipped.cases] = 0

    parameters[, iteri] = c(mu, pi1, pi2)  # Record all updates
    betas[, iteri] = beta
    alfas[, iteri] = alfa
    lambdas[, iteri] = lambda
  }
  results = list("parameters" = parameters,
                 "betas" = betas,
                 "misclassified.cases" = alfas,
                 "misclassified.controls" = lambdas)
  return(results)
}
