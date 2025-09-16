#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
IntegerVector cxx_merge_ranges(const IntegerVector chrom, 
                               const IntegerVector start, 
                               const IntegerVector end, 
                               const IntegerVector strand, 
                               const IntegerVector grp, 
                               const IntegerVector tol) {
  
  const int n=chrom.size();
  const int y=grp.size();
  IntegerVector idx(n);
  int g=0;
  if (n) { 
    g=1;
    idx[0]=1;
  }
  
  for (int i=1; i<n; ++i){
    if (chrom[i]!=chrom[i-1]
        ||(start[i] - end[i-1] - 1)>tol[0]
        ||strand[i]!=strand[i-1]
        ||grp[i]!=grp[i-1]
    ) {
       ++g;
    }
    idx[i]=g;
  }
  return(idx);
}