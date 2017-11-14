# RL-SKAT

RL-SKAT is a method for the calculation of exact p-values for the score test in heritability, in the case of a single kernel and a continuous phenotype. RL-SKAT is described in the following [paper](http://biorxiv.org/content/early/2017/06/11/140889).

Code written by Regev Schweiger. For questions or comments, mail [schweiger@post.tau.ac.il](mailto:schweiger@post.tau.ac.il).

### Installation

You can simply download the code:

1. Download the ZIP file of the latest release here: [https://github.com/cozygene/RL-SKAT/releases](https://github.com/cozygene/RL-SKAT/releases).
2. Extract the ZIP file to a folder

RL-SKAT is implemented as a set of Python classes. To use, simply import/execfile (or just %run in IPython) the file, e.g.:
```python
from rl_skat import *
```

RL-SKAT requires [NumPy](http://www.numpy.org/), [SciPy](http://www.scipy.org) and [FaST-LMM](https://github.com/MicrosoftGenomics/FaST-LMM) (used for the implementation of the Davies method used in the code).

## Usage

An object is constructed which does all the required preprocessing, and can then use it to test multiple phenotypes. All matrices should have the same number of rows (*n*, the sample size).

The file `example.py` contains a short example (for the heritability use-case), along with example files.

## Full kernel

When you have a kernel matrix already computed, use this version:

```python
RL_SKAT_Full_Kernel(kernel_matrix, fixed_covariates=None, add_intercept=True).test(phenotypes, acc=1e-7)
```
Where:
* `kernel_matrix` - An *n by n* kernel matrix.
* `fixed_covariates` - A set of *n by p* optional additional covariates with fixed effects.
* `add_intercept` - Set to True if an intercept should be added as a covariate.
* `phenotypes` - A matrix of *n by N* of *N* phenotypes to test.
* `acc` - Optional parameter of p-value accuracy. More accurate p-values take longer to calculate.

The method returns *N* p-values.

## Set testing

If you have a set of covariates to be treated as having random effects (that is, the square of the kernel matrix), use this  (add covariates as needed):

```python
RL_SKAT_Low_Rank(weighted_Z, fixed_covariates=None, add_intercept=True).test(phenotypes, acc=1e-7)
```

Where:
* `weighted_Z` - An *n by m* matrix of (possibly weighted) covariates corresponding to random effects, e.g. rare variants.
* `fixed_covariates` - A set of *n by p* optional additional covariates with fixed effects.
* `add_intercept` - Set to True if an intercept should be added as a covariate.
* `phenotypes` - A matrix of *n by N* of *N* phenotypes to test.
* `acc` - Optional parameter of p-value accuracy. More accurate p-values take longer to calculate.

The method returns *N* p-values.

## Using in `SOLAR-Eclipse`

[SOLAR-Eclipse](http://solar-eclipse-genetics.org/) is a software package for genetic analysis. SOLAR deals mainly with explicit pedigree files, but it is possible to generate a kinship matrix from it for processing in FIESTA. Whenever you load a pedigree file, solar saves the kinship matrix it computes (see 8.3 [here](https://helix.nih.gov/Documentation/solar-6.6.2-doc/08.chapter.html)). After loading the pedigree and locating the matrix file, run the following Python code:
```python
import numpy as np
import gzip
import numpy.linalg

filename = "<path_to_file>/phi2.gz"
lines = np.vstack([np.array(line.strip().split()[:3]).astype(float) for line in gzip.open(filename).read().splitlines()[1:]])
n_samples = int(max(lines[:,0]))
kinship_matrix = np.zeros((n_samples, n_samples))
for x,y,c in lines:
    kinship_matrix[int(x)-1, int(y)-1] = c

vals, vecs = numpy.linalg.eigh(kinship_matrix)
savetxt("<path>/eigenvalues.txt", vals[::-1])
savetxt("<path>/eigenvectors.txt", vecs[::-1])
```

The output `txt` can now be used in FIESTA.
