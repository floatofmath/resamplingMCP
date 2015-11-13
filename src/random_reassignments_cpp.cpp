#include <Rcpp.h>
using namespace Rcpp;


// In-place permutation of x
static inline int isample(const int n) {return floor(n*unif_rand());}
static inline IntegerVector sample(IntegerVector g) {
    std::random_shuffle(g.begin(),g.end(),isample);
    return g;
}


// [[Rcpp::export]]
IntegerMatrix random_reassignments_cpp(const IntegerVector g, const int nperm=1000) {
    IntegerVector gg = clone(g);
    int i, p=1, n=gg.size();
    IntegerMatrix ans(n,nperm);
    for (i=0; i<nperm; i++) ans(_,i) = sample(gg);
    return ans;
}
