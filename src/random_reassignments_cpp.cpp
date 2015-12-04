#include <Rcpp.h>
using namespace Rcpp;


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
    
