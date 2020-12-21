# Likelihood ratio test to test difference in asymmetries
# Inigo Martincorena - November 2020

# Toy scenario: 3 mutations, 2 from the major branch and 1 from the minor branch, from two bulk samples (i.e. colon and brain)

p1 = 0.85  # true asymmetry in the brain
p2 = 0.7 # true asymmetry in the colon
cov = rep(200, 6) # Coverage vector. If you increase the coverage a lot, you can see how the MLEs below converge to the true values

# Simulating the counts for an experiment

pvec = c(p1*0.5, p1*0.5, (1-p1)*0.5, p2*0.5, p2*0.5, (1-p2)*0.5) # Vector of true p for the 3 mutations in the 2 tissues
xvec = rbinom(n=6, size=cov, prob=pvec)

# Likelihood function to get the MLE from 1 tissue

ll0 = function(xv, cv) {
  pv = seq(0,1,by=0.001) # grid
  ll = sapply(pv, function(p) sum(dbinom(xv, size=cv, prob=c(p*0.5, p*0.5, (1-p)*0.5), log=T))) # likelihood surface
  ml = pv[which.max(ll)] # MLE
  return(ml)
}

ml1 = ll0(xvec[1:3], cov[1:3]) # MLE for the asymmetry in tissue 1
ml2 = ll0(xvec[4:6], cov[4:6]) # MLE for the asymmetry in tissue 2

# Likelihood function to get a single MLE for both tissues (null model)

ll1 = function(xv, cv) {
  pv = seq(0,1,by=0.001) # grid
  ll = sapply(pv, function(p) sum(dbinom(xv, size=cv, prob=c(p*0.5, p*0.5, (1-p)*0.5, p*0.5, p*0.5, (1-p)*0.5), log=T))) # likelihood surface
  ml = pv[which.max(ll)] # MLE
  return(ml)
}

mlboth = ll1(xvec, cov)

# Likelihood Ratio Test. Null model p1==p2. Alternative model p1!=p2.

loglik0 = sum(dbinom(xvec, size=cov, prob=c(mlboth*0.5, mlboth*0.5, (1-mlboth)*0.5, mlboth*0.5, mlboth*0.5, (1-mlboth)*0.5), log=T))
loglik1 = sum(dbinom(xvec, size=cov, prob=c(ml1*0.5, ml1*0.5, (1-ml1)*0.5, ml2*0.5, ml2*0.5, (1-ml2)*0.5), log=T))
pval = 1-pchisq(2*(loglik1-loglik0),1)