// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
using namespace Rcpp;
using namespace arma;

// [[Rcpp::export]]
NumericVector diffmeanC(NumericVector x, IntegerMatrix g){
  // g must be encoded as -1, 1!!
  mat gt = as<mat>(g);
  vec xs = as<vec>(x);
  mat ts = xs.t() * gt;
  return( wrap(ts) );
}


// [[Rcpp::export]]
NumericVector sumdiffC(NumericVector x, IntegerMatrix g){
  mat gt = as<mat>(g);
  mat gc = ones(g.nrow(),g.ncol())-gt;
  vec xs = as<vec>(x);
  mat ts = xs.t() * gt, cs = xs.t() * gc;
  return( wrap(ts-cs) );
}

// [[Rcpp::export]]
NumericVector meandiffC(NumericVector x, IntegerMatrix g){
  mat gt = as<mat>(g);
  mat gc = ones(g.nrow(),g.ncol())-gt;
  vec xs = as<vec>(x);
  mat nt = ones(g.nrow()).t() * gt;
  mat nc = g.nrow() - nt;
  mat ts = xs.t() * gt, cs = xs.t()  * gc;
  return( wrap(ts/nt - cs/nc) );
}


// NumericVector colVarsC(NumericMatrix x){
//   mat xs = as<mat>(x);
  

NumericVector pooled_varianceC(NumericVector x, IntegerMatrix g){
  mat gt = as<mat>(g);
  mat gc = ones(g.nrow(),g.ncol())-gt;
  vec xs = as<vec>(x);
  mat nt = ones(g.nrow()).t() * gt;
  mat nc = g.nrow() - nt;
  mat ts = xs.t() * gt, cs = xs.t()  * gc;
  return( wrap(ts/nt - cs/nc) );
}

// // [[Rcpp::export]]
// NumericVector tstatC(NumericVector x, IntegerMatrix g){
  
