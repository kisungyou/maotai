
<!-- README.md is generated from README.Rmd. Please edit that file -->
Tools for Matrix Algebra, Optimization and Inference Problems
=============================================================

<!-- badges: start -->
[![Travis build status](https://travis-ci.org/kyoustat/maotai.svg?branch=master)](https://travis-ci.org/kyoustat/maotai) [![CRAN status](https://www.r-pkg.org/badges/version/maotai)](https://CRAN.R-project.org/package=maotai) <!-- badges: end -->

`maotai` is an acronym for **M**atrix **A**lgebra, **O**p**T**imization, **A**nd **I**nference problems - though I can't deny motivation from one of [my father's favorite](https://en.wikipedia.org/wiki/Maotai) for the namesake. More detailed introduction will be added later.

Installation
------------

`maotai` released version can be obtained from [CRAN](https://CRAN.R-project.org/package=maotai) with:

``` r
install.packages("maotai")
```

or the up-to-date development version from github:

``` r
## install.packages("devtools")
## library(devtools)
devtools::install_github("kyoustat/maotai")
```

Available Functions
-------------------

Current version of `maotai` supports following functions in the table.

<table style="width:83%;">
<colgroup>
<col width="19%" />
<col width="63%" />
</colgroup>
<thead>
<tr class="header">
<th>Function</th>
<th>Description</th>
</tr>
</thead>
<tbody>
<tr class="odd">
<td><code>bmds</code></td>
<td>Bayesian Multidimensional Scaling</td>
</tr>
<tr class="even">
<td><code>boot.mblock</code></td>
<td>Generate Index for Moving Block Bootstrapping</td>
</tr>
<tr class="odd">
<td><code>boot.stationary</code></td>
<td>Generate Index for Stationary Bootstrapping</td>
</tr>
<tr class="even">
<td><code>cmds</code></td>
<td>Classical Multidimensional Scaling</td>
</tr>
<tr class="odd">
<td><code>dpmeans</code></td>
<td>DP-means Algorithm for Clustering Euclidean Data</td>
</tr>
<tr class="even">
<td><code>ecdfdist</code></td>
<td>Distance Measures between Multiple Empirical CDFs</td>
</tr>
<tr class="odd">
<td><code>ecdfdist2</code></td>
<td>Pairwise Measures for Two Sets of Empirical CDFs</td>
</tr>
<tr class="even">
<td><code>epmeans</code></td>
<td>EP-means Algorithm for Clustering Empirical Distributions</td>
</tr>
<tr class="odd">
<td><code>kmeanspp</code></td>
<td>K-Means++ Clustering Algorithm</td>
</tr>
<tr class="even">
<td><code>lgpa</code></td>
<td>Large-scale Generalized Procrustes Analysis</td>
</tr>
<tr class="odd">
<td><code>lyapunov</code></td>
<td>Solve Lyapunov Equation</td>
</tr>
<tr class="even">
<td><code>matderiv</code></td>
<td>Numerical Approximation to Gradient of a Function with Matrix Argument</td>
</tr>
<tr class="odd">
<td><code>mmd2test</code></td>
<td>Kernel Two-sample Test with Maximum Mean Discrepancy</td>
</tr>
<tr class="even">
<td><code>nem</code></td>
<td>Negative Eigenvalue Magnitude</td>
</tr>
<tr class="odd">
<td><code>nef</code></td>
<td>Negative Eigenfraction</td>
</tr>
<tr class="even">
<td><code>pdeterminant</code></td>
<td>Calculate the Pseudo-Determinant of a Matrix</td>
</tr>
<tr class="odd">
<td><code>shortestpath</code></td>
<td>Find Shortest Path using Floyd-Warshall Algorithm</td>
</tr>
<tr class="even">
<td><code>sylvester</code></td>
<td>Solve Sylvester Equation</td>
</tr>
<tr class="odd">
<td><code>tsne</code></td>
<td>t-Stochastic Neighbor Embedding</td>
</tr>
<tr class="even">
<td><code>trio</code></td>
<td>Trace Ratio Optimization</td>
</tr>
<tr class="odd">
<td><code>weiszfeld</code></td>
<td>Weiszfeld Algorithm for L1-median</td>
</tr>
</tbody>
</table>
