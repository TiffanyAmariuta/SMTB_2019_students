#include <Rcpp.h>
#include <RcppEigen.h>

using namespace Rcpp;
// [[Rcpp::depends(RcppEigen)]]


//

// [[Rcpp::export]]
NumericVector SparseRowVar2(Eigen::SparseMatrix<double> mat,
                            NumericVector mu, 
                            int power){
  mat = mat.transpose();
  NumericVector allVars = no_init(mat.cols());
  for (int k=0; k<mat.outerSize(); ++k){
    double colSum = 0;
    int nZero = mat.rows();
    for (Eigen::SparseMatrix<double>::InnerIterator it(mat,k); it; ++it) {
      nZero -= 1;
      colSum += pow(it.value() - mu[k], power);
    }
    colSum += pow(mu[k], power) * nZero;
    allVars[k] = colSum / (mat.rows() - 1);
  }
  return(allVars);
}


// [[Rcpp::export]]
NumericVector SparseRowVarStd(Eigen::SparseMatrix<double> mat,
                              NumericVector mu,
                              NumericVector sd,
                              double vmax){
  mat = mat.transpose();
  NumericVector allVars(mat.cols());
  for (int k=0; k<mat.outerSize(); ++k){
    if (sd[k] == 0) continue;
    double colSum = 0;
    int nZero = mat.rows();
    for (Eigen::SparseMatrix<double>::InnerIterator it(mat,k); it; ++it)
    {
      nZero -= 1;
      colSum += pow(std::min(vmax, (it.value() - mu[k]) / sd[k]), 2);
    }
    colSum += pow((0 - mu[k]) / sd[k], 2) * nZero;
    allVars[k] = colSum / (mat.rows() - 1);
  }
  return(allVars);
}

// [[Rcpp::export]]
NumericVector SparseRowVarStd2(Eigen::SparseMatrix<double> mat,
                               NumericVector mu,
                               NumericVector var,
                               double varmax, 
                               int power){
  mat = mat.transpose();
  NumericVector allVars(mat.cols());
  for (int k=0; k<mat.outerSize(); ++k){
    if (var[k] == 0) continue;
    double colSum = 0;
    int nZero = mat.rows();
    for (Eigen::SparseMatrix<double>::InnerIterator it(mat,k); it; ++it)
    {
      nZero -= 1;
      colSum += std::min(varmax*var[k], pow(it.value() - mu[k], power));
    }
    colSum += pow(mu[k], power) * nZero;
    allVars[k] = colSum/((mat.rows() - 1)*var[k]);
  }
  return(allVars);
}

// You can include R code blocks in C++ files processed with sourceCpp
// (useful for testing and development). The R code will be automatically 
// run after the compilation.
//
