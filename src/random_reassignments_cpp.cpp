// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
using namespace Rcpp;
using namespace arma;


// In-place permutation of x
static inline int isample(const int n) {return floor(n*unif_rand());}
static inline IntegerVector sample(IntegerVector g) {
    std::random_shuffle(g.begin(),g.end(),isample);
    return g;
}
static inline NumericVector sampleNum(NumericVector xs, int k) {
  NumericVector rs(k);
  std::random_shuffle(xs.begin(),xs.end(),isample);
  rs = head(xs,k);
  return xs;
}


// [[Rcpp::export]]
IntegerMatrix random_reassignments_cpp(const IntegerVector g, const int nperm=1000) {
    IntegerVector gg = clone(g);
    int i, p=1, n=gg.size();
    IntegerMatrix ans(n,nperm);
    for (i=0; i<nperm; i++) ans(_,i) = sample(gg);
    return ans;
}

// [[Rcpp::export]]
NumericMatrix random_samples_cpp(const NumericVector xs, const int k, const int nsam=1000){
  NumericVector xxs = clone(xs), rs(k);
  int i, p=1;
  NumericMatrix ans(k,nsam);
  for (i=0; i<nsam; i++) ans(_,i) = sampleNum(xxs,k);
  return ans;
}


mat combinations_rec(const int n, const int k, vec vs){
  if(k == 1){
    mat rmx(1,n);
    rmx.row(0) = vs.t();
    return rmx;
  } else if(k == n) {
    mat rmx(k,1);
    rmx.col(0) = vs;
    return rmx;
  } else {
    int nr = Rf_choose(n-1,k-1);
    mat v0mx(1,nr);
    vec vts = vs.tail(n-1);
    v0mx.fill(vs(0));
    return( join_rows(join_cols(v0mx,combinations_rec(n - 1,k - 1,vts)),
                      combinations_rec(n-1,k,vts))) ;
  }
}


// [[Rcpp::export]]
IntegerMatrix combinations_cpp(const int n, const int k){
  vec vs = linspace<vec>(1, n, n);
  return( wrap(combinations_rec(n,k,vs)) );
}

// [[Rcpp::export]]
NumericMatrix subsamples_cpp(NumericVector xs, const int k){
  vec vs = as<vec>(xs);
  int n = xs.size();
  return( wrap(combinations_rec(n,k,xs)) );
}

mat bincombinations_rec(const int k){
  if(k == 1){
    mat vmx =linspace<mat>(0,1,2);
    return vmx.t();
  } else {
    int k2 = pow(2,k-1);
    mat v0(1,k2);
    mat v1(1,k2);
    v0.fill(0);
    v1.fill(1);
    return( join_rows(join_cols(v0,bincombinations_rec(k-1)),
                      join_cols(v1,bincombinations_rec(k-1))) );
  }
}

// [[Rcpp::export]]
IntegerMatrix bincombinations_cpp(const int p){
  IntegerMatrix retval(pow(2,p),p);
  int n;
  for (n=1;n<p+1;n++) {
    int p2 = pow(2,p);
    int k2 = p2/pow(2,n);
    IntegerVector v(2*k2);
    std::fill(v.begin(),v.end(),1);
    std::fill_n(v.begin(),k2,0);
    retval(_, n-1) = rep(v,p2);
  }
  return retval;
}
