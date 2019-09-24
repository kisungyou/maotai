#include <RcppArmadillo.h>

// [[Rcpp::depends(RcppArmadillo)]]

using namespace arma;

// [[Rcpp::export]]
arma::mat cppsub_2007Wang(arma::mat V0, int mm, int d, arma::mat Spu, arma::mat Stu, int maxiter, double eps){
  // 1. preliminary setup
  double abstol = std::sqrt((static_cast<double>(mm*d))*eps);
  
  arma::mat Vold = V0;
  arma::mat Vnew(mm,d,fill::zeros);
  
  double lbdn = 0.0;
  double incV = 0.0;
  
  arma::vec Vval(mm,fill::zeros);
  arma::mat Vvec(mm,mm,fill::zeros);
  arma::mat Vtmp(mm,d,fill::zeros);
  
  arma::mat Stv(mm,mm,fill::zeros);
  
  // 2. do the iteration
  for (int i=0;i<maxiter;i++){
    // 2-1. compute the lambda
    lbdn = arma::trace(Vold.t()*Spu*Vold)/arma::trace(Vold.t()*Stu*Vold);
    // 2-2. solve the eigenvalue problem
    eig_sym(Vval,Vvec,(Spu-(lbdn*Stu)));
    // 2-3. extract Vtmp
    Vtmp = Vvec.tail_cols(d);
    // 2-4. update by readjusting
    Stv  = Vtmp*Vtmp.t()*Stu*Vtmp*Vtmp.t();
    eig_sym(Vval,Vvec,Stv);
    Vnew = Vvec.tail_cols(d);
    // 2-5. updating info and update
    incV = arma::norm(Vold-Vnew,"fro");
    Vold = Vnew;
    // 2-6. stop if the criterion is met
    if (incV < abstol){
      break;
    }
    // no need to print; Rcpp::Rcout << "* iteration " << i+1 << " complete.." << std::endl;
  }
  
  // 3. return the value
  return(Vold);
}
