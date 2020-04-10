# News for Package maotai

### changes in version 0.1.5
  * Following functions are added,
    - `bmds`     : Bayesian Multidimensional Scaling.
    - `cmds`     : Classical Multidimensional Scaling.
    - `kmeanspp` : k-means++ Algorithm for Clustering.
    - `nem`      : Negative Eigenvalue Magnitude.
    - `nef`      : Negative Eigenfraction.
    - `tsne`     : t-Stochastic Neighbor Embedding.
  * `distgmm` is removed for better composition of the package.

### changes in version 0.1.4
  * Update `README` for listing the available functions in the package.
  * Corrected an Armadillo type-casting error.
  * Corrected example visualization parameter settings.
  
### changes in version 0.1.3
  * Following functions are added,
    - `boot.mblock`     : generate index for Moving Block Bootstrapping.
    - `boot.stationary` : generate index for Stationary Bootstrapping.
    - `epmeans`         : EP-means algorithm for clustering empirical distributions.
    - `weiszfeld`       : Weiszfeld algorithm for computing L1-median.
    
### changes in version 0.1.2
  * Following functions are added,
    - `mmd2test`        : kernel two-sample test with Maximum Mean Discrepancy.
    
### changes in version 0.1.1
  * Following functions are added,
    - `dpmeans`         : DP-means algorithm for clustering Euclidean data.
    - `distgmm`         : distance measures between multiple samples using Gaussian mixture models.
    - `ecdfdist`        : distance measures between multiple empirical cumulative distribution functions.
    
### changes in version 0.1.0
  * Package is first deployed.
  * Initialize the following documentation,
    - NEWS for keeping record of updates.
    - README to briefly introduce the method.
  * Initial set of functions are listed as follows,
    - `lgpa`            : large-scale generalized Procrustes analysis.
    - `lyapunov`        : solve Lyapunov equation.
    - `matderiv`        : numerical approximation to gradient of a function with matrix argument.
    - `pdeterminant`    : calculate the pseudo-determinant of a matrix.
    - `shortestpath`    : find shortest path length using Floyd-Warshall algorithm.
    - `sylvester`       : solve Sylvester equation.
    - `trio`            : trace ratio optimization.