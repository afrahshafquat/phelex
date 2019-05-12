#' Function to remove missing values from matrix
#'
#' Replaces NA values from columns of matrix with column mean.
#'
#' @param x Genotype matrix with dimensions n x m
#'
#' @return Genotype matrix with dimension n x m with NA values replaced with
#' the corresponding column mean.
#'
#' @keywords replace,missing,impute
#'
#' @export
remove_na = function(x) {
  for (i in 1:ncol(x)) {
    x[is.na(x[, i]), i] = mean(x[, i], na.rm = T)
  }
  return(x)
}

#' Estimate modal values
#'
#' Estimates mode for each parameter in matrix given vector of values.
#'
#' @param betas Matrix with dimensions n x m (parameters x iterations).
#' @param method Default value is shorth. See more in @details .
#'
#' @details Available methods are "mfv", "lientz", "naive",
#' "venter", "grenander", "hsm", "hrm", "parzen", "tsybakov", and "asselin".
#' For more, see \code{\link[modeest]{mlv}} documentation.
#'
#' @return Vector of estimated mode values for each parameter (n)
#' @keywords mode,parameter,estimation
#'
#' @import modeest
#' @export
estimate_betas = function(betas,  method = 'shorth') {
  betas = betas
  final_betas = apply(betas, 1, modeest::mlv, method=method)
  return(final_betas)
}

#' Simulates misclassification in binary phenotype
#'
#' This function simulates misclassification of binary phenotype (vector y) by switching fraction lambda of 0 (controls)
#' to 1 (cases) and fraction c_alpha of 1 (cases) to 0 (controls).
#'
#' @param y Vector of values of 0s and 1s.
#' @param lambda Fraction of controls switched to cases. Default is NULL (no switch).
#' @param c_alpha Fraction of cases switched to controls. Default is NULL (no switch).
#' @param seed Seed number to ensure different switches take place everytime. Default is 1000.
#'
#' @return Vector with switched 0s and 1s according to defined lambda and alpha.
#' @keywords misclassification,phenotype,simulation
#'
#' @export
perturb_y = function(y,
                     lambda = NULL,
                     c_alpha = NULL,
                     seed = 1000) {
  ybar = y  # True y
  set.seed(seed = seed)
  n = length(ybar)

  case.inds = which(y == 1)
  control.inds = which(y == 0)

  if (! is.null(c_alpha)) {
    flipcases = floor(c_alpha * length(case.inds))
    flipcase.inds = sample(case.inds, flipcases)
    y[flipcase.inds] = 0
  }
  if (! is.null(lambda)) {
    flipcontrols = floor(lambda * length(control.inds))
    flipcontrol.inds = sample(control.inds, flipcontrols)
    y[flipcontrol.inds] = 1
  }
  return(y)
}

#' Log sum of values x and y
#'
#' Given a = log(x) and b = log(y), returns log(x + y).
#'
#' @param a Log value of x.
#' @param b Log value of y.
#'
#' @return Returns log(x + y).
#' @keywords keywords
#'
#' @export
logsum = function(a, b) {
  #log sum of log(a+b) where a and b are in log.
  if(a > b){
    p = a + log1p(exp(a) - exp(b))
  } else{
    p = b + log1p(exp(b) - exp(a))
  }
  return(p)
}

#' Calculates misclassification probability for samples
#'
#' Function estimates the probability
#' that a case is not a true case or a control is not a true control with the
#' matrix of misclassified.cases (misclassified.samples or misclassified.controls)
#'
#' @param misclassified.samples Matrix of misclassification indicators (misclassified.samples or misclassified.controls)
#' @param standardize Flag to normalize probability values (only when they are not normalized). Recommended to keep default.
#'
#' @return Misclassification probabilities for samples
#'
#' @keywords misclassification,probability
#'
#' @export
estimate_misclassification_probability = function(misclassified.samples,
                         standardize=TRUE) {
  n = ncol(misclassified.samples)
  misclassified.p.samples = rowSums(misclassified.samples)/n
  if(standardize & (mean(misclassified.p.samples) > 0.5)){
    misclassified.p.samples = 1-misclassified.p.samples
  }
  return(misclassified.p.samples)
}

#' Computes modified phenotype
#'
#' Function computes corrected phenotype using misclassification probability for each sample.
#'
#'
#' @param y Vector of observed phenotype
#' @param misclassified.p.cases vector of misclassification probabilities for cases. If absent, only controls are switched
#' @param misclassified.p.controls vector of misclassification probabilities for controls. If absent, only cases are switched
#' @param case.threshold Threshold = percentile cut-off on misclassification probabilities. All cases with Pr(misclassification) > percentile threshold are switched to controls.
#' @param control.threshold Threshold = percentile cut-off on misclassification probabilities. All controls with Pr(misclassification) > percentile threshold are switched to cases.
#'
#' @return Modified Phenotype according to misclassification probabilities provided
#'
#' @keywords corrected,phenotype,misclassification
#'
#' @export
get_phenotype = function(misclassified.p.cases=NULL,
                         misclassified.p.controls=NULL,
                         y,
                         case.threshold=NULL,
                         control.threshold=NULL) {
  case.inds = which(y == 1)
  control.inds = which(y == 0)
  corrected.phenotype = y

  # Switch cases to controls
  if(!is.null(misclassified.p.cases)){
    if(is.null(case.threshold)) case.threshold = .99
    flip.inds = which(misclassified.p.cases >= quantile(misclassified.p.cases, case.threshold))
    corrected.phenotype[case.inds[flip.inds]] = 0
  }

  # Switch controls to cases
  if(!is.null(misclassified.p.controls)){
    if(is.null(control.threshold)) control.threshold = .99
    flip.inds = which(misclassified.p.controls >= quantile(misclassified.p.cases, case.threshold))
    corrected.phenotype[control.inds[flip.inds]] = 1
  }

  return(corrected.phenotype)
}
