from rl_skat import *

K = loadtxt("kinship.txt")
Y = loadtxt("phenotypes.txt")
RL = RL_SKAT_Full_Kernel(K, fixed_covariates=None, add_intercept=True)
pvals = RL.test(Y, acc=1e-7)
print pvals

