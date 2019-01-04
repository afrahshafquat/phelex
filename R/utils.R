#' Function to remove missing values from matrix
#'
#' Replaces NA values from columns of matrix with column mean.
#'
#' @param x Genotype matrix with dimensions s x n.
#'
#' @return Genotype matrix with dimension s x n with NA values replaced with
#' the corresponding column mean.
#'
#' @keywords keywords
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
#' @param betas Matrix with dimensions n x s (parameters x samples).
#' @param method Default value is shorth. See more in @details .
#'
#' @details Available methods are "mfv", "lientz", "naive",
#' "venter", "grenander", "hsm", "hrm", "parzen", "tsybakov", and "asselin".
#' For more, see \code{\link[modeest]{mlv}} documentation.
#'
#' @return Vector of estimated mode values for each parameter (n)
#' @keywords keywords
#'
#' @export
estimate_betas = function(betas,  method = 'shorth') {
  betas = betas
  final_betas = rep(0, nrow(betas))
  for(k in 1:nrow(betas)){
    final_betas[k] = get_beta(betas[k, ], method = method)
  }
  return(final_betas)
}

#' Switches binary status of values
#'
#' For vector y, fraction lambda of 0s are changed to 1s and fraction alpha of 1s to 0s.
#'
#' @param y Vector of values of 0s and 1s.
#' @param lambda Fraction of controls switched to cases. Default is NULL (no switch).
#' @param c_alpha Fraction of cases switched to controls. Default is NULL (no switch).
#' @param seed Seed number to ensure different switches take place everytime. Default is 1000.
#'
#' @return Vector with switched 0s and 1s according to defined lambda and alpha.
#' @keywords keywords
#'
#' @export
perturb_y = function(y,
                     lambda = NULL,
                     c_alpha = NULL,
                     seed = 1000) {
  ybar = y  # True y
  set.seed(seed = seed)
  n = length(ybar)

  case_inds = which(y == 1)
  control_inds = which(y == 0)

  if (! is.null(c_alpha)) {
    flipcases = floor((1 - c_alpha) * length(case_inds))
    flipcase_inds = sample(case_inds, flipcases)
    y[flipcase_inds] = 0
  }
  if (! is.null(lambda)) {
    flipcontrols = floor(lambda * length(control_inds))
    flipcontrol_inds = sample(control_inds, flipcontrols)
    y[flipcontrol_inds] = 1
  }
  return(y)
}

#' Return mode for a vector of values
#'
#' Given vector of values beta, the function returns the mode as
#' estimated by \code{\link[modeest]{mlv}}.
#'
#' @param beta Vector of values..
#' @param method Default value is shorth. See more in @details .
#'
#' @details Available methods are "mfv", "lientz", "naive",
#' "venter", "grenander", "hsm", "hrm", "parzen", "tsybakov", and "asselin".
#' For more, see \code{\link[modeest]{mlv}} documentation.
#'
#' @return Estimated mode value.
#' @keywords keywords
#'
#' @import modeest
#' @export
get_beta = function(beta, method = 'shorth') {
  est_beta = modeest::mlv(beta, method = 'shorth')
  # est_beta = est_beta$M
  return(est_beta)
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
#' Given model defined by Shafquat et al (...), the function estimates the probability
#' that a case is not a true case or a control is not a true control with the
#' matrix of flip indicators (flip.cases or flip.controls) from plvm_lmm/plvm/gibbs/gibbs_lmm
#'
#' @param flip.indicators Matrix of flip indicators (flip.cases or flip.controls)
#' @param standardize Flag to normalize probability values (only when they are not normalized). Recommended to keep default.
#'
#' @return Misclassification probabilities for samples
#'
#' @keywords keywords
#'
#' @export
estimate_flip_probability = function(flip.indicators,
                         standardize=TRUE) {
  n = ncol(flip.indicators)
  flip.p = rowSums(flip.indicators)/n
  if(standardize & (mean(flip.p) > 0.6)){
    flip.p = 1-flip.p
  }
  return(flip.p)
}

#' Computes modified phenotype
#'
#' Given model defined by Shafquat et al (...), function computes
#' modified phenotype using misclassification probability for each sample.
#'
#' @param y Vector of observed phenotype
#' @param flip.p.cases vector of misclassification probabilities for cases. If absent, only controls are switched
#' @param flip.p.controls vector of misclassification probabilities for controls. If absent, only cases are switched
#' @param case.threshold All cases with Pr(misclassification) >= threshold are switched to controls.
#' Default value is mean(probability) + 2*sd(probability).
#' @param control.threshold All controls with Pr(misclassification) >= threshold are switched to cases
#'
#' @return Modified Phenotype according to misclassification probabilities provided
#'
#' @keywords keywords
#'
#' @export
get_phenotype = function(flip.p.cases=NULL,
                         flip.p.controls=NULL,
                         y,
                         case.threshold=NULL,
                         control.threshold=NULL) {
  case.inds = (y == 1)
  control.inds = which(!case.inds)
  case.inds = which(case.inds)
  corrected_phenotype = y

  # Switch cases to controls
  if(!is.null(flip.p.cases)){

    if(is.null(case.threshold)){  # Provide threshold if not specified
      case.threshold = mean(flip.p.cases) + (2 * sd(flip.p.cases))
    }

    flip.inds = which(flip.p.cases >= case.threshold)
    corrected_phenotype[case.inds[flip.inds]]=0
  }

  # Switch controls to cases
  if(!is.null(flip.p.controls)){

    if(is.null(control.threshold)){  # Provide threshold if not specified
      control.threshold = mean(flip.p.controls) + (2 * sd(flip.p.controls))
    }

    flip.inds = which(flip.p.controls >= control.threshold)
    corrected_phenotype[control.inds[flip.inds]] = 1
  }

  return(corrected_phenotype)
}
