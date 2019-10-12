#ifndef _MAOTAI_EVALUATIONS_H
#define _MAOTAI_EVALUATIONS_H

#define ARMA_NO_DEBUG

#include <RcppArmadillo.h>

using namespace Rcpp;
using namespace arma;

double eval_gaussian(arma::vec x, arma::vec mu, arma::mat cov);
double eval_gmm(arma::vec x, arma::mat mus, arma::cube covs, arma::vec weight);
arma::vec eval_gaussian_data(arma::mat X, arma::vec mu, arma::mat cov);
arma::vec eval_gmm_data(arma::mat X, arma::mat mus, arma::cube covs, arma::vec weight);

#endif
