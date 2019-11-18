// [[Rcpp::depends(RcppArmadillo)]]
#include<RcppArmadillo.h>
using namespace Rcpp;
using namespace arma;

// This is a simple example of exporting a C++ function to R. You can
// source this function into an R session using the Rcpp::sourceCpp 
// function (or via the Source button on the editor toolbar). Learn
// more about Rcpp at:
//
//   http://www.rcpp.org/
//   http://adv-r.had.co.nz/Rcpp.html
//   http://gallery.rcpp.org/
//

// [[Rcpp::export]]

arma::vec dvec_C(const arma::mat &L){
  // This functin generates p * (p +1) / 2 which stores the lower triangular of
  // matrix S by stucking diagonals and subdiagonals sequentiall
  int p  = L.n_rows;
  int start = 0;   
  int end = p - 1;
  arma::vec x; x.zeros(p * (p + 1) / 2);
  for(int j=0; j<(p);j++){
    int ind = p - j ;
    x(span(start, end)) = L.diag(-j);
    start = start + ind;
    end = end + ind - 1;
  }
  return(x);
}


// [[Rcpp::export]]

arma::vec offsubsum(const List &A,const arma::vec &x,const int &i, const arma::mat& S){
  int p = S.n_rows;
  arma::vec vi = A[i-1];
  int li = vi.n_elem;
  arma:: vec summ; summ.zeros(li);
  for (int j=1; j<(p + 1);j++){
    arma::uvec vj = A[j - 1];
    int lj = vj.n_elem;
    arma::vec xj = x((vj - 1));
    if(j > i){
      if(j != (p)){
        summ(span((j - i), (li - 1))) = summ(span((j-i), (li - 1))) + S(span((j - i), (li - 1)), span(0, (lj - 1))).diag() % xj;
      }else{
        summ(span((j-i), (li - 1))) = summ(span((j-i), (li - 1))) + 
              S((j-i), (lj - 1)) * xj;
      }
    } else if (j < i){
      if (i != p){
        summ = summ + S(span((i - j), (lj - 1)), span(0, (li - 1))).diag() % xj(span(i - j, lj - 1));
      }else{
       summ = summ + S(span((i - j), (lj - 1)), span(0, (li - 1))) * xj(span(i - j, lj - 1));
      }
    }

  }
  return summ;
}

