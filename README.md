# Causal Discovery with Equal Variance Assumption 

## Abstract

Prior work has shown that causal structure can be uniquely identified from observational data when these follow a structural equation model whose error terms have equal variances. We show that this fact is implied by an ordering among (conditional) variances. We demonstrate that ordering estimates of these variances yields a simple yet state-of-the-art method for causal structure learning that is readily extendable to high-dimensional problems.

This repository maintains the code for this project. 

## Installation

- The `EqVarDAG` package:

  The package contains the four methods (top-down and bottom-up approach, for low-dimensional settings, and high-dimensional settings.) Among the methods, high-dimensional bottom-up is implemented with the CLIME estimator, as proposed in [Ghoshal and Honorio (2018)][http://proceedings.mlr.press/v84/ghoshal18a/]. The package can be installed by `devtools` 

  ```
  # install.packages("devtools")
  devtools::install_github("WY-Chen/EqVarDAG")
  ```

- Simulation studies

  In low dimensional settings, we compare our algorithms with the Greedy DAG Search algorithm (GDS) proposed by [Peters and B&uuml;hlmann (2014)][https://doi.org/10.1093/biomet/ast043]. The GDS algorithm is implemented in the supplement of said paper. A slightly modified version is maintained in this repository. 

  ```r
  # download this repository
  setwd(dir = "/some/path/")
  download.file(url = "https://github.com/WY-Chen/EqVarDAG/archive/master.zip",
                destfile = "EqVarDAG.zip")
  # unzip the .zip file
  unzip(zipfile = "EqVarDAG.zip")
  setwd(dir = "/some/path/EqVarDAG")
  ```

## Methods

Let X be an n-by-p matrix generated according to an SEM with equal variance, we can estimate the topological ordering and the corresponding DAG structure using any of the following:

```R
# low-dimensional, top-down
EqVarDAG_TD(X)
# low-dimensional, bottom-up
EqVarDAG_BU(X)
# high-dimensional, top-down, 
#		assuming maximum in-degree = J
EqVarDAG_HD_TD(X,J)
# high-dimensional, bottom-up, 
#		using CLIME with cv-ed lambda
EqVarDAG_HD_CLIME(X,cv=TRUE,NULL)
# high-dimensional, bottom-up, 
#		using CLIME with given lambda
EqVarDAG_HD_CLIME(X,cv=FALSE,lambdafix)
```

## Simulation studies

The details of the simulation settings are specified in our paper. 

 Low dimensional simulation with N replicates, each problem with p variables and n data can be performed by calling

```R
# run simulations
source("/some/path/EqVarDAG/Experiments/sim_script_lowD.R")
simulation_lowD(p,n,N) 			# low-dimensional settings
```

High dimensional simulation is implemented with varies graph types, specified in the paper. We can call `type`  with choices `'er'`, `'chain'`, or `'hub'`. 

```R
# run simulations
source("/some/path/EqVarDAG/Experiments/sim_script_highD.R")
simulation_highD(p,n,N,type) 	# low-dimensional settings
```

