# RL-SKAT

To use, simply import/execfile (or just %run in IPython) the file.

You construct an object which does all the required preprocessing, and can then use it to test multiple phenotypes.

All matrices should have the same number of rows (sample size).

## Full kernel

When you have a kernel matrix already computed, use this version:

```
RL_SKAT_Full_Kernel(kernel_matrix, fixed_covariates, add_intercept=True).test(phenotypes)
```

## Set testing

If you have a set of covariates to be treated as having random effects (that is, the square of the kernel matrix), use this:

```
RL_SKAT_Low_Rank(weighted_Z, fixed_covariates, add_intercept=True).test(phenotypes)
```
