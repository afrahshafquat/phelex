load('~/Desktop/response/sim0419.RData')
load('~/Desktop/response/sim0419.1000.3000.50.RData')
load('~/Desktop/response/sim0419.1000.3000.6.3.50.3.input.RData')

betas = dataset$betas[50,which(dataset$disease.inds %in% c(22640,50015,71884,74615))]

x = X.subset
y = y
A = A
iterations = 1e5
alpha.prior = c(10,1)
lambda.prior = c(1,1)
u = matrix(dataset$u[phen.inds], ncol=1)
link = 'pnorm'

verbose = TRUE
stamp = 1e4

markers = ncol(x)  # Number of markers
x = remove_na(x)  # Removes missing values
x = x

y = y
case.inds = which(y == 1)
control.inds = which(y == 0)
n = length(y)

iterations = iterations
num.params = 3  # sigmaA.index alpha, lambda,

chol.A = chol(A)
invA = chol2inv(A)

lsi = 0

compute_posterior = function(plogis.x, alpha, lambda){
  ll = sum(log((plogis.x * (alpha ^ y) * ((1 - alpha) ^ (1 - y)))+
                 ((1 - plogis.x) * (lambda ^ y) * ((1 - lambda) ^ (1 - y)))))
  prior = sum(dbeta(alpha, shape1 = alpha.prior[1], shape2 = alpha.prior[2], log = T))#,
  #              dbeta(lambda, shape1 = lambda.prior[1], shape2 = lambda.prior[2], log = T))
  # )
  energy = ll + prior
  return(energy)
}

compute_liability = function(betas){
  lj = x %*% betas + u
  lj = scale(lj)
  return(lj)
}
  
parameters = matrix(NA,nrow=10,ncol=iterations)

liab.x = compute_liability(betas)
plogis.x = pnorm(liab.x)

rm(dataset)
  # u = rep(0, n)
for(ii in 1:3){
  print(ii)
  alpha = 1#rbeta(1, 1, 1)
  lambda = 0.03#rbeta(1, 1, 1)
  sigmaA = 0.01#sigmaA.initial
  
  current.post = compute_posterior(plogis.x, alpha, lambda)
  
  if(verbose) print(paste('Warmups for Beta, Alpha and Lambda', date()))
  
  # Adaptive Metropolis Hastings sampling algorithm for initial estimates for parameters
  for(i in 1:iterations) {
    if((!(i %% 1e4)) & verbose) {
      print(paste(i, date()))
    }
    # 
    # # MCMC Proposal step
    # proposal.alpha = alpha#truncdist::rtrunc(1, 'norm', a = 0, b = 1, mean = alpha, sd = 0.1)
    # proposal.lambda = lambda#truncdist::rtrunc(1, 'norm', a = 0, b = 1, mean = lambda, sd = 0.1)
    # 
    # #Posterior Proposal
    # proposal.post = compute_posterior(plogis.x, proposal.alpha, proposal.lambda) #ll1 + prior1
    # 
    # #Accept/Reject
    # p = exp(proposal.post - current.post)
    # 
    # if (runif(1) < p) {
    #   alpha = proposal.alpha
    #   # lambda = proposal.lambda
    #   current.post = proposal.post
    # }
    
    # u.tmp = u
    # lambda.A = (1 / sigmaA)
    # 
    # for(j in 1:n) {
    #   inv.u = 1 / (1 + (invA[j, j] * lambda.A))
    #   ci_i = invA[j,, drop = FALSE]
    #   ci_i[j] = 0
    #   u.i = u.tmp
    #   u.i[j] = 0
    #   uj = inv.u * ((liab.x[j] - (x[j, ] %*% beta)) - (lambda.A * ci_i %*% u.i))
    #   u.tmp[j] = rnorm(1, mean = uj, sd = sqrt(inv.u))
    # }
    # u = u.tmp
    sigmaA = 1/(truncdist::rtrunc(1,'gamma', shape=(n-2)/2, rate=(t(u) %*% invA %*% u)/2, a=0.01)) ##Truncated 0-1000
    parameters[ii,i] = sigmaA
  }
}

save(parameters, file='~/Desktop/response/sigmaAconstantu.estimates.RData')



parameters[, 1] = c(sigmaA, alpha, lambda)

lsi = 0
beta.jump.sd = exp(2*-1.49)

if(verbose) print(paste('Starting sampling for all parameters', date()))

# Adaptive Metropolis Hastings sampling algorithm within Gibbs
for(i in 2:iterations) {
  if ((!(i %% stamp)) & verbose) {
    print(paste(i, date()))
  }

  # MCMC Proposal step
  proposal.alpha = truncdist::rtrunc(1, 'norm', a = 0, b = 1, mean = alpha, sd = 0.1)
  proposal.lambda = truncdist::rtrunc(1, 'norm', a = 0, b = 1, mean = lambda, sd = 0.1)

  #Posterior Proposal
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

    #Gibbs Sampling for random effect
    # u.tmp = u
    # lambda.A = (1 / sigmaA)
    # 
    # for(j in 1:n) {
    #   inv.u = 1 / (1 + (invA[j, j] * lambda.A))
    #   ci_i = invA[j,, drop = FALSE]
    #   ci_i[j] = 0
    #   u.i = u.tmp
    #   u.i[j] = 0
    #   uj = inv.u * ((liab.x[j] - (x[j, ] %*% beta)) - (lambda.A * ci_i %*% u.i))
    #   u.tmp[j] = rnorm(1, mean = uj, sd = sqrt(inv.u))
    # }
    # u = u.tmp
    sigmaA = 1/(truncdist::rtrunc(1,'gamma', shape=(n-2)/2, rate=(t(u) %*% invA %*% u)/2, a=0.01)) ##Truncated 0-1000

    liab.x = compute_liability(beta)
    plogis.x = pnorm(liab.x)
    current.post = compute_posterior(plogis.x, beta, alpha, lambda)
  }

  betas[, i] = beta
  energy[i] = current.post
  parameters[, i] = c(sigmaA, alpha, lambda)

  flip.alpha = log(alpha) + log(plogis.x) - log(1 - plogis.x) - log(lambda)
  flip.alpha = 1 / (1 + exp(flip.alpha)) #probability for misclassification in observed cases

  flip.lambda = log(1 - plogis.x) + log(1 - lambda) - log(1 - alpha) - log(plogis.x)
  flip.lambda = 1 / (1 + exp(flip.lambda)) #probability for misclassification in observed controls

  alfa.i = rbinom(length(case.inds), 1, flip.alpha[case.inds])
  lambda.i = rbinom(length(control.inds), 1, flip.lambda[control.inds])

  alfas[, i] = alfa.i
  lambdas[, i] = lambda.i
}
result = list("betas" = betas,
              "parameters" = parameters,
              "posterior" = energy,
              "accept" = accept.rate,
              "misclassified.cases" = alfas,
              "misclassified.controls" = lambdas)