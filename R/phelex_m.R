#' Estimate misclassification in phenotype
#'
#' Estimates misclassification probabilities in observed GWAS phenotype y given genotypes dataset x.
#' The method follows the PheLEx-m algorithm to predict misclassification probabilities
#' using Adaptive Metropolis-Hastings defined by Shafquat et al.
#'
#' @param y Phenotype vector with length n.
#' @param x Genotype matrix with dimensions n x m.
#' @param iterations Total number of iterations for Metropolis-Hastings Sampling
#' @param beta.initial.vec Vector of initial values for beta parameters i.e. effect sizes for n snps.
#' @param alpha.prior Beta prior parameters for true-positve rate
#' @param lambda.prior Beta prior parameters for false-positive rate.
#' @param beta.prior 'norm'(default): Normal prior or 'unif' Uniform prior on fixed effects
#' @param beta.prior.params if beta.prior is norm, then (mean, sd), if unif then (min, max)
#' @param stamp Iteration breakpoint to print time
#' @param link probit or logistic model (options pnorm and plogis respectively)
#' @param verbose Default TRUE. Prints progress information
#' @param normalize Default FALSE. Scales liability computed
#' @param mu.update Fractions of iterations to start updating mean fixed effects or mu. Default = 0.5
#'
#' @return  List containing
#' \itemize{
#'  \item betas: Matrix of estimated effect sizes for each SNP (SNPs[rows] x iterations[columns]).
#'  \item parameters: Matrix with estimated parameter values[rows] across iterations[columns].
#'  Order is c(mu, pi1, pi2) where pi1/pi2 = false positive/false negative rates, mu = mean effect size value
#'  \item misclassified.cases: Matrix of alpha vectors where 1s represent false positives and 0s represent true positives as inferred at each iterations
#'  \item misclassified.controls: Matrix of lambda vectors where 1s represent false negatives and 0s represent true negatives as inferred at each iterations
#'  \item posterior: Vector of posterior probability across iterations
#'  \item accept: Vector of 1/0 values across iterations; 1 indicates proposal was accepted at iteration;0 o/w
#'  }
#'
#' @keywords phelex_m,misclassification,GWAS,phenotype,adaptive metropolis hastings
#'
#' @import MASS
#' @import utils
#' @import truncdist
#' @export
#'
phelex_m = function(x,
                y,
                iterations = 1e5,
                alpha.prior = c(1, 1),
                lambda.prior = c(1, 1),
                link = 'pnorm',
                beta.prior='norm',
                beta.prior.params = c(1, 3),
                beta.initial.vec = NULL,
                mu.update = 0.5,
                verbose = TRUE,
                normalize = FALSE,
                stamp = 1e3) {

  markers = ncol(x)  # Number of markers
  x = remove_na(x)  # Removes missing values
  x = x

  y = y
  case.inds = which(y == 1)
  control.inds = which(y == 0)
  n = length(y)

  link = ifelse(link == 'plogis', plogis, pnorm)
  beta.prior = ifelse(beta.prior == 'unif', dunif, dnorm)

  iterations = iterations
  num.params = 3  # mu, alpha, lambda

  accept.rate = rep(0, iterations)  # Acceptance of proposal sampled in MH
  energy = rep(0, iterations)  # Posterior
  betas = matrix(0, nrow = markers, ncol = iterations)  # Fixed effects
  parameters = matrix(0, nrow = num.params, ncol = iterations)  # Parameters mu, alpha, lambda
  alfas = matrix(0, nrow = length(case.inds), ncol = iterations) # Misclassification indicators in cases
  lambdas = matrix(0, nrow = length(control.inds), ncol = iterations) # Misclassification indicators in controls

  alpha.index = num.params - 1
  lambda.index = num.params
  mu.index = 1

  if (is.null(beta.initial.vec)) {  # Initialize parameter values
    beta.initial.vec = rnorm(markers, 0, .1)
  }

  lsi = 0
  beta.jump.sd = exp(2*-1.49)
  beta.prior.param1 = beta.prior.params[1]
  beta.prior.param2 = beta.prior.params[2]

  compute_posterior = function(plogis.x, betas, alpha, lambda){
    ll = sum(log((plogis.x * (alpha ^ y) * ((1 - alpha) ^ (1 - y)))+
                   ((1 - plogis.x) * (lambda ^ y) * ((1 - lambda) ^ (1 - y)))))
    prior = sum(c(beta.prior(betas, beta.prior.param1, beta.prior.param2, log = T),
                  dbeta(alpha, shape1 = alpha.prior[1], shape2 = alpha.prior[2], log = T),
                  dbeta(lambda, shape1 = lambda.prior[1], shape2 = lambda.prior[2], log = T)))
    energy = ll + prior
    return(energy)
  }

  compute_liability = function(betas, normalize=normalize){
    lj = x %*% betas + mu
    if(normalize) lj = scale(lj)
    return(lj)
  }

  mu.start = ceiling((1-mu.update)*iterations)
  mu = 0
  alpha = rbeta(1, 1, 1)
  lambda = rbeta(1, 1, 1)
  beta = beta.initial.vec[]
  liab.x = compute_liability(beta)
  plogis.x = link(liab.x)

  parameters[, 1] = c(mu, alpha, lambda)
  betas[, 1] = beta
  current.post = compute_posterior(plogis.x, beta, alpha, lambda)
  energy[1] = current.post

  for(i in 2:iterations) {  # Adaptive Metropolis Hastings sampling algorithm
    if ((! i %% stamp) & verbose) {
      print(paste(i, date()))
    }
    if (! i %% 1e2) { # Adaptive part of MH
      #update beta.jump.sd
      a.rate = (sum(accept.rate[(i-99):i])/100)
      delta.n = min(0.01, ((i/100)^-.5))
      if(a.rate < 0.2){
        lsi = lsi - delta.n
      }else{
        lsi = lsi + delta.n
      }
      beta.jump.sd = max(0.005, exp(2 * (-1.49 + lsi)))
    }

    # MCMC Proposal step
    proposal.beta = rnorm(markers, mean = beta, sd = beta.jump.sd)

    proposal.alpha = truncdist::rtrunc(1, 'norm', a = 0, b = 1, mean = alpha, sd = 0.1)
    proposal.lambda = truncdist::rtrunc(1, 'norm', a = 0, b = 1, mean = lambda, sd = 0.1)

    #Posterior Proposal
    proposal.liab = compute_liability(proposal.beta)
    proposal.plogis = link(proposal.liab)
    proposal.post = compute_posterior(proposal.plogis, proposal.beta, proposal.alpha, proposal.lambda) #ll1 + prior1

    #Accept/Reject
    p = exp(proposal.post - current.post)

    if (runif(1) < p) {
      accept.rate[i] = 1
      beta = proposal.beta
      alpha = proposal.alpha
      lambda = proposal.lambda
      plogis.x = proposal.plogis
      liab.x = proposal.liab
      current.post = proposal.post
    }

    if ((i > mu.start & (! i %% 10))){
    #update mu using gibbs
      mu = sum(liab.x - (x %*% beta)) / n
      mu = rnorm(1, mean = mu, sd = 1 / sqrt(n))
      liab.x = compute_liability(beta)
      plogis.x = link(liab.x)
      current.post = compute_posterior(plogis.x, beta, alpha, lambda)
    }

    flip.alpha = log(alpha) + log(plogis.x) - log(1 - plogis.x) - log(lambda)
    flip.alpha = 1 / (1 + exp(flip.alpha)) #probability for misclassification in observed cases

    flip.lambda = log(1 - plogis.x) + log(1 - lambda) - log(1 - alpha) - log(plogis.x)
    flip.lambda = 1 / (1 + exp(flip.lambda)) #probability for misclassification in observed controls

    alfa.i = rbinom(length(case.inds), 1, flip.alpha[case.inds])
    lambda.i = rbinom(length(control.inds), 1, flip.lambda[control.inds])

    alfas[, i] = alfa.i
    lambdas[, i] = lambda.i
    betas[, i] = beta
    energy[i] = current.post
    parameters[, i] = c(mu, alpha, lambda)
  }
  result = list("betas" = betas,
                "parameters" = parameters,
                "posterior" = energy,
                "accept" = accept.rate,
                "misclassified.cases" = alfas,
                "misclassified.controls" = lambdas)
  return(result)
}
