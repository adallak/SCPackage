// [[Rcpp::depends(RcppArmadillo)]]
#include<RcppArmadillo.h>
using namespace Rcpp;
using namespace arma;
using namespace std;


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

////////////////////////////////////////////////////
////////////////////////////////////////////////////
///////////////////////////////////////////////////

//Diagonal Update 

// [[Rcpp::export]]
arma::vec diag_update(const List A, arma::vec x,arma::mat S,int i=1){
  // return sqrt(square(offsubsum(A,x,1,S))+4*S.diag());
  arma::vec sum = offsubsum(A,x,i,S);
  return ((-1 * sum) + sqrt(square(sum) + 4 * S.diag())) / (2 * S.diag());
}


/////Lasso soft thresholding function

// [[Rcpp::export]]

arma::vec soft_threshold(arma::colvec y, double lambda1){
  // the last element is the same as beta(unpenalized)
  int r = y.n_elem;       //y.n_elem;
  arma::vec result(r);//(r ) = y(r );
  // This performs the soft-threshold for beta
  // for the first r - 1 elements
  for (int l = 0; l < (r); l++){
    // element-wise soft-threshold
    if (y(l) > lambda1)
      result(l) = y(l) - lambda1;
    else if (y(l) < -lambda1)
      result(l) = y(l) + lambda1;
    else
      result(l) = 0.0;
  }
  return result;
}

// [[Rcpp::export]]

arma::vec weighted_soft_threshold(arma::vec &y,const arma::vec &weight, const double &lambda1){
  // the last element is the same as beta(unpenalized)
  int r = y.n_elem;
  arma::vec weight_lambda = weight * lambda1;
  arma::vec result(r);
  // arma::vec weightedsoft(r);//(r ) = y(r );
  // This performs the soft-threshold for beta
  // for the first r - 1 elements
  for (int l = 0; l < (r); l++){
    // element-wise soft-threshold
    if (y(l) > weight_lambda(l))
      result(l) = y(l) - weight_lambda(l);
    else if (y(l) < -weight_lambda(l))
      result(l) = y(l) + weight_lambda(l);
    else
      result(l) = 0.0;
  }
  return result;
}


// [[Rcpp::export]]

List matGenerateC(int p){
  List A(p) ;
  //  std::vector< std::vector<int> > ND; 
  int start ; start = 0;
  int end ; end = p;
  for (int i = 0; i < (p); i++){
    int ind = p - i ;
    arma::vec temp = zeros(p - i);
    for(int j = 0; j < (p - i + 1); j++){
      temp[j] = start + j + 1;
    }
    start = start + ind ;
    end = end + ind ;
    A[i] = temp;
    //    std::generate(temp.begin(), temp.end(), [n = start]() mutable { return n++; });
  }
  return A;
}


// [[Rcpp::export]]

arma::colvec fused_coef(std::vector<double> y, const double lambda2){
  Environment pkg = Environment::namespace_env("flsa");
  Function f = pkg["flsa"];
  Rcpp::NumericVector fs = wrap(f(y, Named("lambda2", lambda2)));
  return  fs;
}




// [[Rcpp::export]]
arma::colvec trend_coef(arma::vec y, double lambda2){
  
  Environment pkg = Environment::namespace_env("l1tf");
  Function f = pkg["l1tf"];
  Rcpp::NumericVector fs = f(y, Named("lambda", lambda2));
  return  fs;
}

// [[Rcpp::export]]
arma::mat getd1d(int j){
  
  Environment pkg = Environment::namespace_env("genlasso");
  Function g = pkg["getD1d"];
  arma::mat D = as<mat>(g(j));
  //  NumericMatrix res = g(fs,Named("lambda",lambda2));
  // Executing Matrix( m, sparse = TRIE )
  return  D;
}

// [[Rcpp::export]]
arma::mat getdtf(const int j, const int ord){
  
  Environment pkg = Environment::namespace_env("genlasso");
  Function g = pkg["getDtf"];
  arma::mat D = as<mat>(g(j, Named("ord",ord)));
  //  NumericMatrix res = g(fs,Named("lambda",lambda2));
  // Executing Matrix( m, sparse = TRIE )
  return  D;
}


// [[Rcpp::export]]
arma::vec fused_update(arma::vec x, arma::mat S, const double lambda1, const double lambda2,int band, List A){
  int p= S.n_rows;
  arma::vec Bii;
  arma::vec sqrt_Bii;
  arma::vec sqrt_Bii_inv;
  std::vector<double> y;
  arma::vec temp;
  arma::vec weight;
  arma::vec temp_y;
  //Rcpp::NumericVector y;
  //  arma::mat z=x;
  for (int i = 2; i < (band + 1); i++){
    arma::vec tm = A[(i-1)];
    int li = tm.n_elem;
    Bii = S.diag();
    Bii = Bii.rows(0, (li - 1));
    //    Bi = diagmat(sqrt(Bii.diag()));
    sqrt_Bii = sqrt(Bii);
    sqrt_Bii_inv = 1 / sqrt_Bii ;
    if(i > (p-1)){
      arma::uvec ind = A[i-1];
      x.elem(ind - 1).zeros();
    }else{
      arma::uvec ind = A[i - 1];
      temp_y = (-1) * sqrt_Bii_inv % offsubsum(A, x, i, S);
      y = arma::conv_to< std::vector<double>  >::from(temp_y);// wrap(temp_y);
      arma::vec x_i = fused_coef(y, lambda2);
      x.elem(ind - 1) = sqrt_Bii_inv % x_i;
    }
    if (lambda1 > 0){
      arma::uvec ind = A[i - 1];
      temp = x.elem(ind - 1);
      weight = 1 / (2 * Bii);
      x.elem(ind-1) = weighted_soft_threshold(temp, weight, lambda1);
      //       temp.print();
    }
  }
  arma::uvec ind1 = A[0];
  x.elem(ind1-1)= diag_update(A,x,S);
  return x;
}
// 

// [[Rcpp::export]]
Rcpp::List iter_fused(arma::vec x,arma::mat S,const double lambda1, const double lambda2,int band, List A,const double max_iter,double  ABSTOL   ){
  //  double  ABSTOL   = 1e-3;
  // double  RELTOL   = 1e-4;
  arma::vec history(max_iter);
  arma::vec oldL;
  arma::vec x_temp;
  arma::vec vecL;
  //int p = S.n_rows;
  for (int iter =0; iter<(max_iter);iter++)
  {
    oldL = x;
    x_temp = fused_update(x, S, lambda1, lambda2, band, A);
    x = x_temp;
    vecL = x;
    history(iter)  = arma::norm((vecL - oldL),"inf");
    //   Rcpp::Rcout <<history(iter) << std::endl;
    //  double current_iter = arma::as_scalar(history(iter));
    if (history(iter) <= ABSTOL){
      break;
    }
    if(iter==(max_iter - 1))
    {
      Rcpp::Rcout << "SSC fails to converge" << std::endl;
    }
  }
  
  return(Rcpp::List::create(Rcpp::Named("history")=history,
                            Rcpp::Named("x")=x));
  //     return x;
}


// 
// [[Rcpp::export]]

arma::vec trend_update(arma::vec x, arma::mat S, const double lambda1, const double lambda2, int band, List A){
  int p= S.n_rows;
  arma::vec Bii;
  arma::vec sqrt_Bii;
  arma::vec sqrt_Bii_inv;
  int li;
  arma::vec weight;
  arma::vec temp;
  arma::vec y;
  //  arma::mat z=x;
  for (int i = 2; i < (band + 1); i++){
    arma::vec tm = A[(i-1)];
    li = tm.n_elem;
    Bii = S.diag();
    Bii = Bii.rows(0, (li - 1));
    //    Bi = diagmat(sqrt(Bii.diag()));
    sqrt_Bii = sqrt(Bii);
    sqrt_Bii_inv = 1 / sqrt_Bii ;
    if(i >= (p-1)){
      arma::uvec ind = A[i-1];
      //     ind.print();
      x.elem(ind - 1).zeros();
    }else{
      arma::uvec ind=A[i-1];
      y= (-1) * sqrt_Bii_inv % offsubsum(A, x, i, S);
      x.elem(ind -1) = trend_coef(y, lambda2);
      x.elem(ind - 1) = sqrt_Bii_inv % x.elem(ind - 1);
      //      x.elem(ind-1).print();
    }
    if (lambda1 > 0){
      arma::uvec ind = A[i - 1];
      temp = x.elem(ind - 1);
      weight = 1 / (2 * Bii);
      x.elem(ind-1) = weighted_soft_threshold(temp, weight, lambda1);
      //       temp.print();
    }
  }
  arma::uvec ind1 = A[0];
  x.elem(ind1 - 1) = diag_update(A, x, S);
  
  return x;
}

// [[Rcpp::export]]

Rcpp::List iter_trend(arma::vec x,arma::mat S, const double lambda1, const double lambda2,int band, List A,const double max_iter,double  ABSTOL   ){
  //  double  ABSTOL   = 1e-3;
  // double  RELTOL   = 1e-4;
  arma::vec    history(max_iter);
  arma::vec oldL;
  arma::vec x_temp;
  arma::vec vecL;
  //int p = S.n_rows;
  for (int iter = 0; iter < (max_iter); iter++)
  {
    oldL = x;
    x_temp = trend_update(x, S, lambda1, lambda2, band, A);
    x = x_temp;
    vecL = x;
    history(iter)  = arma::norm((vecL - oldL),"inf");
    //   Rcpp::Rcout <<history(iter) << std::endl;
    //  double current_iter = arma::as_scalar(history(iter));
    if (history(iter - 1) <= ABSTOL){
      break;
    }
    if(iter == (max_iter))
    {
      Rcpp::Rcout << "SSC fails to converge" << std::endl;
    }
  }
  
  return(Rcpp::List::create(Rcpp::Named("history")=history,
                            Rcpp::Named("x")=x));
  //     return x;
}
// 
//// HP update
// [[Rcpp::export]]
arma::vec hp_coef(arma::mat Bii, arma::mat D,arma::vec y, double lambda1,double lambda2){
  arma::vec z = solve((2 * Bii + 2 * lambda2 * D.t() * D), y);
  //  arma::vec soft=soft_threshold(y,lambda1);
  return soft_threshold(z, lambda1);
}


// [[Rcpp::export]]
arma::vec hp_update(arma::vec &x, arma::mat S, const double lambda1, const double lambda2, int band, List A){
  int p = S.n_rows;
  arma::mat Bii;
  arma::vec diag;
  int li;
  arma::mat y;
  //  arma::mat z=x;
  for (int i=2; i<(band+1);i++){
    arma::vec tm = A[(i-1)];
    li = tm.n_elem;
    diag = S.diag();
    diag = diag.rows(0, (li - 1));
    Bii = diagmat(diag);
    y=(-2 * (offsubsum(A, x, i, S)));
    //   Bi = diagmat(sqrt(Bii.diag()));
    //    sqrt_Bi = diagmat(1/(Bi.diag()));
    if(i==(p)){
      arma::uvec ind = A[i - 1];
      int D = 1;
      //     ind.print();
      x.elem(ind - 1)= soft_threshold(y, lambda1)/(2 * Bii + 2 * lambda2 * std::pow(D, 2));
    }else if(i == (p - 1)){
      arma::uvec ind = A[i - 1];
      arma::mat D = getd1d(ind.n_elem);
      x.elem(ind - 1) = hp_coef(Bii, D, y, lambda1, lambda2);
    }else{
      arma::uvec ind = A[i-1];
      arma::mat D = getdtf(ind.n_elem, 1);
      x.elem(ind - 1) = hp_coef(Bii, D, y, lambda1, lambda2);
    }
    // if (lambda1>0){
    //    arma::uvec ind=A[i-1];
    //  arma::vec weight = 1/(2*Bii.diag());
    //arma::vec temp = x.elem(ind-1);
    //x.elem(ind-1)= weighted_soft_threshold(temp,weight,lambda1);
    //       temp.print();
  }
  arma::uvec ind1 = A[0];
  x.elem(ind1 - 1)= diag_update(A, x, S);
  
  return x;
}

// [[Rcpp::export]]
Rcpp::List iter_hp(arma::vec &x,arma::mat S,const double lambda1, const double lambda2,int band, List A,const double max_iter,double  ABSTOL   ){
  //  double  ABSTOL   = 1e-3;
  // double  RELTOL   = 1e-4;
  arma::vec    history(max_iter);
  //int p = S.n_rows;
  for (int iter =0; iter < (max_iter); iter++)
  {
    arma::vec oldL = x;
    arma::vec x_temp = hp_update(x, S, lambda1, lambda2, band, A);
    x = x_temp;
    arma::vec vecL = x;
    history(iter)  = arma::norm((vecL - oldL),"inf");
    //   Rcpp::Rcout <<history(iter) << std::endl;
    //  double current_iter = arma::as_scalar(history(iter));
    if (history(iter) <= ABSTOL){
      break;
    }
    if(iter == (max_iter - 1))
    {
      Rcpp::Rcout << "SSC fails to converge" << std::endl;
    }
  }
  
  return(Rcpp::List::create(Rcpp::Named("history") = history,
                            Rcpp::Named("x") = x));
  //     return x;
}


