
#include <RcppArmadillo.h>

// [[Rcpp::depends(RcppArmadillo)]]

//' Variance-covariance matrix of random components in LMM
//'
//' @description
//' \eqn{V_lambda_tau} returns the variance-covariance matrix
//' of the sum of polygenic effects and observational noise
//' in a genetic LMM
//'
//' @param K The kinship matrix in genetic LMM
//' @param log_lambda_tau The vector \eqn{(log \lambda, log \tau) }, where \eqn{\lambda}
//' and \eqn{\tau} are the two parameters used in the variance components 
//' of the genetic LMM.
//' 
//' @details
//' The variance covariance matrix is used in finding the 
//' RMLE of \eqn{\lambda} and \eqn{\tau}. This eventually 
//' helps us to find the MLE of SNP effect size. For more details,
//' see our article. 
//' 
//' @return \eqn{(1/\tau)(\lambda K + I)    }
//'
//' @export
// [[Rcpp::export]]


arma::mat V_lambda_tau(arma::mat K, arma::vec log_lambda_tau) {
  int R=K.n_rows ;
  return   pow(exp(log_lambda_tau[1]),-1)*(exp(log_lambda_tau[0])*K  + arma:: eye(R,R) )   ;
}

// [[Rcpp::depends(RcppArmadillo)]]
//' MLE of fixed effect sizes under known variance components
//'
//' @description
//' \eqn{gamma_tilda} returns the MLE of covariate and SNP effect sizes
//' in a genetic LMM when the variance components of polygenic
//' effect and observational noise are known.
//'
//' @param K The kinship matrix in genetic LMM
//' @param X The augmented matrix of fixed effects which is obtained by
//' augmenting the genotype vector of the SNP to the covariate matrix
//' @param y The phenotype vector (should be a column matrix)
//' @param log_lambda_tau The vector \eqn{(log \lambda, log \tau) }, where \eqn{\lambda}
//' and \eqn{\tau} are the two parameters used in the variance components 
//' of the genetic LMM.
//' 
//' @details
//' The MLE of covariate and SNP effect sizes for known values of \eqn{\lambda}
//' and \eqn{\tau} is used in finding the RMLE of \eqn{(\lambda, \tau)}
//' by maximizing the marginal likelihood. This in turn helps to find the unrestricted MLE
//' of SNP effect size. For more details, see our article. 
//' 
//' @return MLE of \eqn{\gamma}, the vector of fixed effect sizes 
//' for known value of \eqn{(\lambda, \tau)}
//'
//' @export
// [[Rcpp::export]]
arma::vec gamma_tilda(arma::mat K,arma::mat X,arma::vec y, arma::vec log_lambda_tau ){
  arma::mat A= V_lambda_tau(K, log_lambda_tau);
  arma::mat B= arma::inv_sympd(A);
  return( arma::inv_sympd(X.t() * B * X)  * (X.t() * B *y)           );
}

// [[Rcpp::depends(RcppArmadillo)]]
//' Marginal likelihood of variance parameters in genetic LMM
//'
//' @description
//' \eqn{neglr_lambda_tau} returns the (negated) marginal likelihood
//' of the parameters used in the variance component of the 
//' genetic LMM
//'
//' @param log_lambda_tau The vector \eqn{(log \lambda, log \tau) }, where \eqn{\lambda}
//' and \eqn{\tau} are the two parameters used in the variance components 
//' of the genetic LMM.
//' @param K The kinship matrix in genetic LMM
//' @param X The augmented matrix of fixed effects which is obtained by
//' augmenting the genotype vector of the SNP to the covariate matrix
//' @param y The phenotype vector (should be a column matrix)
//' 
//' @details
//' The marginal likelihood obtained is minimized (since it is negated)
//' to find the RMLE of \eqn{(\lambda, \tau)}.
//' For more details, see our article. 
//' 
//' @return Marginal likelihood of \eqn{(\lambda, \tau)}
//'
//' @export
// [[Rcpp::export]]
double neglr_lambda_tau(arma::vec log_lambda_tau, arma::mat K,arma::mat X,arma::vec y ){
  arma::mat A= V_lambda_tau(K, log_lambda_tau);
  arma::mat B= arma::inv_sympd(A);
  arma::vec b= gamma_tilda(K,X,y,log_lambda_tau);
  arma::mat C= 0.5*(arma::log_det_sympd(A)+ log_det_sympd(X.t() * B *X) + (y- X*b ).t() * B  * (y- X*b ) )  ;
  return(C(0,0));
}

// [[Rcpp::depends(RcppArmadillo)]]
//' MLE and standard deviation of SNP effect size in genetic LMM
//'
//' @description
//' \eqn{beta_hat_sigma_hat} returns the maximum likelihood estimate of 
//' SNP effect size and its standard deviation. 
//'
//' @param K The kinship matrix in genetic LMM
//' @param X The augmented matrix of fixed effects which is obtained by
//' augmenting the genotype vector of the SNP to the covariate matrix
//' @param y The phenotype vector (should be a column matrix)
//' @param log_lambda_tau_hat The \eqn{2}-vector consisting of the log of 
//' RMLE estimates of \eqn{\lambda} and \eqn{\tau} 
//' 
//' @details
//' The MLE and standard deviation of SNP effect size is used in performing 
//' the hypothesis testing of no association. For more details, see our article. 
//' 
//' @return A vector of two elements where the first element is the MLE,
//' and the second element is the standard deviation of the SNP effect size
//'
//' @export
// [[Rcpp::export]]
arma::vec beta_hat_sigma_hat(arma:: mat K, arma::mat X, arma::vec y, arma::vec log_lambda_tau_hat ){
  int R=X.n_cols;
  arma::mat A= V_lambda_tau(K, log_lambda_tau_hat);
  arma::mat B= arma::inv_sympd(A);
  arma::mat C=X.t()*B* X ;
  arma::mat D=arma::inv_sympd(C);
  arma::vec E= D * X.t() * B * y;
  int S=E.n_elem;
  arma::vec out={E(S-1), pow( D(R-1, R-1),0.5)};
  return out;
}