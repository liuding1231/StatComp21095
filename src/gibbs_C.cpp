#include <Rcpp.h>
using namespace Rcpp;
//' @title A Gibbs sampler using Rcpp
//' @description A Gibbs sampler using Rcpp
//' @param N the number of samples
//' @return a random matrix of nrow \code{N} and ncol 2
//' @examples
//' \dontrun{
//' x<-gibbs_C(5000)
//' }
//' @export
//[[Rcpp::export]]
NumericMatrix gibbs_C(int N){
  NumericMatrix mat(N, 2);
  int x0=5;
  double y0=0.5;
  int n=9;
  int a=2;
  int b=3;
  for (int i=0; i<N;i++) {
    int x = rbinom(1, n, y0)[0];
    mat(i,0)=x;
    x0=x;
    double y = rbeta(1, x0+a, n-x0+b)[0];
    mat(i,1)=y;
    y0=y;
  }
  return(mat);
}
