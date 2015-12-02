#include <Rcpp.h>
#include <stdio.h>
using namespace Rcpp;
using namespace arma;

unsigned nChoosek( unsigned n, unsigned k )
{
    if (k > n) return 0;
    if (k * 2 > n) k = n-k;
    if (k == 0) return 1;

    int result = n;
    for( int i = 2; i <= k; ++i ) {
        result *= (n-i+1);
        result /= i;
    }
    return result;
}


// [[Rcpp::export]]
IntegerVector nextCombination(IntegerVector combination,int m)
{
  int i = combination.length()-1;
  while (i >=0 && combination(i) == m)
    {
      printf("%d\n",i);
      i -= 1;
      m -= 1;
    }
  if(i >= 0){
    combination(i) = combination(i)+1;
  }
  return(combination);
}


mat armaComb(int n, int k, mat vs){
  int r = nChoosek(n,k);
  mat ans;v0(0,0), vr(1,n), vc(n,1), ans(r,n);
  if(r == 0){
    return v0;
  } else if (r == 1) {
    mat.insert_cols(0,vs);
    return mat;
  } else if (k == n) {
    vr(_,1) = vr;
    return vr;
  } else {
    v.fill(vs(0));
    ans(_,1) = v1rep;
    ans( Range(0,n-1), Range(1,k) ) = recComb(n-1,k-1,vs(Range(1,k)));
    ans( Range(n,r-1), Range(0,k) ) = recComb(n-1,k,tail(vs,n-1));
    return ans;
  }
}
      

// [[Rcpp::export]]
IntegerMatrix recComb(int n, int k, IntegerVector v){
  vec vs = as<vec>v;
  return armaComb(n,k,vs)
}

// [[Rcpp::export]]
IntegerMatrix combinationsC(const int n, const int k)
{
  int i, r = nChoosek(n,k);
  IntegerVector cc = seq_len( k );
  IntegerMatrix ans(r,k);
  for(i = 0;i<r;++i)
    {
      ans(i,_) = cc;
      cc = nextCombination(cc,n);
    }
  return(ans);
  }


