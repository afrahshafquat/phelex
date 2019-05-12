library(phelex)
library(PhenotypeSimulator)
source('~/Desktop/mplvm_pipeline/simulate_dataset.R')
glmApply=function(x,y){
  model = summary(glm(true.y~x, family = binomial()))
  return(model$coefficients[,4])
}

dataset = simulate_dataset(dsnps = 10,
                           n = 500,
                           nsnps = 20,
                           beta.a.mean = 5,
                           npheno = 1,
                           pop = F)
x = dataset$x
true.y = dataset$phenotypes[,1]
misclassified.y = perturb_y(true.y, lambda = .2)
flip.inds = which((misclassified.y - true.y)==1)
A = getKinship(N=nrow(x),X = x, standardise = F)

pvalues.orig = glmApply(x=x, y= true.y)
plot(-log10(pvalues.orig));points(dataset$disease.inds+1, -log10(pvalues.orig)[dataset$disease.inds+1], col='red')

pvalues.mod = glmApply(x=x, y= misclassified.y)
plot(-log10(pvalues.mod));points(dataset$disease.inds+1, -log10(pvalues.mod)[dataset$disease.inds+1], col='red')


phelex.input = x[, dataset$disease.inds]

phelex.results = phelex(x=phelex.input, y=misclassified.y, A=A, stamp=5e4, alpha.prior = c(10,1), iterations=2e5) #err
phelexm.results = phelex_m(x=phelex.input, y=misclassified.y, stamp=5e4, alpha.prior = c(100,1), iterations = 2e5, mu.update = .2) #done
phelexmh.results = phelex_mh(x=phelex.input, y=misclassified.y, A=A, stamp = 1e4, iterations = 5e4) #done
rekaya.results = rekaya(x=phelex.input, y=misclassified.y, stamp=1e4, iterations = 5e4) #done

misclass.p.phelex = estimate_misclassification_probability(phelex.results$misclassified.cases)
misclass.p.phelexm = estimate_misclassification_probability(phelexm.results$misclassified.cases)
misclass.p.phelexmh = estimate_misclassification_probability(phelexmh.results$misclassified.cases)
misclass.p.rekaya = estimate_misclassification_probability(rekaya.results$misclassified.cases)

corrected.phelex = get_phenotype(misclassified.p.cases = misclass.p.phelex, y=misclassified.y)
corrected.phelexm = get_phenotype(misclassified.p.cases = misclass.p.phelexm, y=misclassified.y)
corrected.phelexmh = get_phenotype(misclassified.p.cases = misclass.p.phelexmh, y=misclassified.y)
corrected.rekaya = get_phenotype(misclassified.p.cases = misclass.p.rekaya, y=misclassified.y)

z=rowSums(phelex.results$misclassified.cases[,1e5:2e5])/1e5; plot(z);points(flip.inds, z[flip.inds], col='red')
z=rowSums(phelexm.results$misclassified.cases[,1e5:2e5])/1e5; plot(z);points(flip.inds, z[flip.inds], col='red')
z=rowSums(phelexmh.results$misclassified.cases[,1e4:5e4])/4e4; plot(z);points(flip.inds, z[flip.inds], col='red')
z=rowSums(rekaya.results$misclassified.cases[,1e4:5e4])/4e4; plot(z);points(flip.inds, z[flip.inds], col='red')
