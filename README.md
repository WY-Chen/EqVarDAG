# Causal Discovery with Equal Variance Assumption 

## Abstract

Prior work has shown that causal structure can be uniquely identified from observational data when these follow a structural equation model whose error terms have equal variances. We show that this fact is implied by an ordering among (conditional) variances. We demonstrate that ordering estimates of these variances yields a simple yet state-of-the-art method for causal structure learning that is readily extendable to high-dimensional problems.

This repository maintains the code for this project. 

## Installation

- The `EqVarDAG` package:

  The package contains the four methods (top-down and bottom-up approach, for low-dimensional settings, and high-dimensional settings.) Among the methods, high-dimensional bottom-up is implemented with the CLIME estimator, as proposed in [Ghoshal and Honorio (2018)][http://proceedings.mlr.press/v84/ghoshal18a/]. The package can be installed by `devtools` 

  ```R
  # install.packages("devtools")
  devtools::install_github("WY-Chen/EqVarDAG")
  ```

## Methods

Let X be an n-by-p matrix generated according to an SEM with equal variance, we can estimate the topological ordering and the corresponding DAG structure using any of the following:

```R
# low-dimensional, top-down
EqVarDAG_TD(X)
# low-dimensional, bottom-up
EqVarDAG_BU(X)
# high-dimensional, top-down, assuming maximum in-degree <= J
EqVarDAG_HD_TD(X,J)
# high-dimensional, bottom-up, using CLIME, as in [Ghoshal and Honorio (2018)]
EqVarDAG_HD_CLIME(X)
```

Specifially, the arguments to the function is
```R
EqVarDAG_TD(
  X,
  mtd = "ztest",
  alpha = 0.05,
  threshold = 0.1,
  FCD = NULL,
  precmtd = NULL
)
```
- X: n by p data matrix
- mtd: DAG inference method, we implemented the following options:
  - "ztest": (only for p<n) Simultaneous tests of partial correlations, see [Drton and Perlman.2007];
  - "dlasso": debiased lasso (default with FCD=True and precmtd="sqrtlasso", see below for details); 
  - "lasso": lasso with fixed lambda, see [Shojaie and Michailidis. 2010]; 
  - "adalasso": adaptive lasso with fixed lambda, see [Shojaie and Michailidis. 2010]; 
  - "cvlasso": cross-validated lasso; 
  - "scallasso": scaled lasso;
  - "rls": (only for p<n) fit recursive least squares and threshold the regression coefs 
  - "chol": (only for p<n) perform cholesky decomposition and threshold the regression coefs 
- alpha: level for DAG inference step
- threshold: threshold for "rls" or "chol" methods
- FCD: only used in debiased lasso
  - TRUE: use the FCD procedure [Javanmard and Montanari. 2018]
  - FALSE: use individual tests to select
- precmtd: only used in debiased lasso, choose how to compute debiasing matrix 
  - "cv": node-wise lasso with joint 10 fold cv 
  - "sqrtlasso": square-root lasso (no tune)

We proposed a two step procedure for DAG learning under equal error-variance. First, a topological ordering is learned via the iterative algorithm; Then, the DAG is inferred using the learned ordering. 
Since the main contribution of this paper is the first step (the second step is a well-studied problem, see e.g.,  [Shojaie and Michailidis. 2010]), in the original version of this paper/package we only implemented/reported results of the "cvlasso" method. 
However, as we realized people not only use this package to learn orderings but also DAGs, we implemented a variety of DAG learning algorithms to choose from. In particular, we recommend using the simultaneous test of partial correlation (__"ztest"__) in low dimensional problems; for high dimensional problems, we recommend using the FCD procedure with debiased lasso whose debiasing matrix is computed by square-root lasso (__"dlasso"__ with FCD=TRUE and precmtd="sqrtlasso"). These two methods both consistently recovers the correct gragh given true ordering. 


## Simulation studies

The details of the simulation settings are specified in our paper. 

In low dimensional settings, we compare our algorithms with the Greedy DAG Search algorithm (GDS) proposed by [Peters and B&uuml;hlmann (2014)][https://doi.org/10.1093/biomet/ast043]. The GDS algorithm is implemented in the supplement of said paper. A slightly modified version is maintained in this repository. 

  ```R
  # download this repository
  setwd(dir = "/some/path/")
  download.file(url = "https://github.com/WY-Chen/EqVarDAG/archive/master.zip",
                destfile = "EqVarDAG.zip")
  # unzip the .zip file
  unzip(zipfile = "EqVarDAG.zip")
  setwd(dir = "/some/path/EqVarDAG-master")
  ```

 Low dimensional simulation with N replicates, each problem with p variables and n data can be performed by calling

```R
# run simulations
source("/some/path/EqVarDAG/Experiments/sims_low.R")
simulation_lowD(p,n,N) 			# low-dimensional settings
```

High dimensional simulation is implemented with varies graph types, specified in the paper. We can call `type`  with choices `'er'`, `'chain'`, or `'hub'`. 

```R
# run simulations
source("/some/path/EqVarDAG/Experiments/sims_high.R")
simulation_highD(p,n,N,type) 	# high-dimensional settings
```
