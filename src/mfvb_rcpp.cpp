//#include <Rcpp.h>
#include<RcppArmadillo.h>
//[[Rcpp::depends(RcppArmadillo)]]
//using namespace Rcpp;
//using namespace arma;

// This is a simple example of exporting a C++ function to R. You can
// source this function into an R session using the Rcpp::sourceCpp
// function (or via the Source button on the editor toolbar). Learn
// more about Rcpp at:
//
//   http://www.rcpp.org/
//   http://adv-r.had.co.nz/Rcpp.html
//   http://gallery.rcpp.org/
//

// potential issue with template
template <typename T>
Rcpp::NumericVector arma2vec(const T& x) {
  return Rcpp::NumericVector(x.begin(), x.end());
}

// Rcpp::NumericVector arma2vec(arma::vec x) {
//   return Rcpp::NumericVector(x.begin(), x.end());
// }


// // [[Rcpp::export]]
arma::mat create_Z(arma::mat fir_l, arma::mat sec_l) {
  int p = fir_l.n_cols, q = sec_l.n_cols, n = fir_l.n_rows;
  arma::mat out(n,p*q);
  int tmp =0;
  for (int i = 0; i < p; i++) {
    for (int j = 0; j < q; j++){
      out.col(tmp) = fir_l.col(i) % sec_l.col(j);
      tmp= tmp+1;
    }

  }
  return out;
}


// // [[Rcpp::export]]
arma::vec update_tau(int p, int q,
               double a_tau, double b_tau,
               arma::vec phi_old, arma::vec e_b_sq_s1, arma::vec e_tau_old){
  arma::mat e_b_sq_new = (phi_old % e_b_sq_s1 + (1-phi_old)/e_tau_old);
  e_b_sq_new.reshape(q,p);
  arma::vec sum_e_b_sq = sum(e_b_sq_new,1);
  arma::vec e_tau_j = a_tau/(0.5*sum_e_b_sq + b_tau);
  arma::mat e_tau_new_tmp;
  for (int i=0; i<p;i++){
    e_tau_new_tmp.insert_cols(i, e_tau_j);
  }
  return arma::vectorise(e_tau_new_tmp);

}


// // [[Rcpp::export]]
arma::vec update_pi(double a_pi, double b_pi,
              arma::vec phi_old){
  Rcpp::NumericVector b = arma2vec(phi_old);
  Rcpp::NumericVector e_logit_pi_new_tmp = Rcpp::digamma(arma2vec(b + a_pi))-
    Rcpp::digamma(arma2vec(b_pi - b));
  return Rcpp::as<arma::vec>(e_logit_pi_new_tmp);
}


// // [[Rcpp::export]]
double update_lambda(arma::mat ztz,
                     arma::vec e_beta, arma::vec var_beta,
                     int n, double len_res_sq){
  // NOTE: This function call R-function!!!
  double e_zbeta_sq = sum(arma::diagvec(ztz) % var_beta) + arma::as_scalar(e_beta.t() * ztz * e_beta);
  Rcpp::Function update_lambda_r("update_lambda_r");
  return Rcpp::as<double>(update_lambda_r(Rcpp::Named("n") = n, Rcpp::Named("e_zbeta_sq")=e_zbeta_sq,
                                    Rcpp::Named("len_res_sq") = len_res_sq));

}


// // [[Rcpp::export]]
arma::vec update_mu(int p, int q, double e_inv_lambda_new,
              arma::vec phi_old, arma::vec mu_old, arma::vec e_beta, arma::vec sigma_sq_new, arma::vec response,
              arma::mat Z){
  arma::vec e_m_all = Z * e_beta;
  arma::vec e_m_nok = e_m_all;
  arma::vec mu_new = mu_old;
  for (int k=0; k<p*q;k++){
    e_m_nok = e_m_all - (phi_old[k]*mu_old[k])*Z.col(k);
    mu_new[k] = -sigma_sq_new[k] *
      arma::as_scalar((response + e_inv_lambda_new * e_m_nok).t() * Z.col(k));
    e_m_all = e_m_nok + (phi_old[k]*mu_new[k])*Z.col(k);
  }
  return mu_new;
}


// // [[Rcpp::export]]
arma::vec update_phi(arma::vec e_logit_pi_new, arma::vec sigma_sq_new, arma::vec e_tau_new, arma::vec mu_new){
  arma::vec phi_logit =  e_logit_pi_new + 0.5 * log(sigma_sq_new) +
    0.5 * log(e_tau_new) + 0.5 * square(mu_new)/sigma_sq_new;
  return 1/(1+1/exp(phi_logit));
}

// [[Rcpp::export]]
Rcpp::List mfvb(arma::vec response, arma::mat fir_l, arma::mat sec_l, double a_tau, double b_tau,
          double a_pi, double b_pi, int max_iter, double max_tol){
  // Fixed values
  int p = fir_l.n_cols, q = sec_l.n_cols, n = fir_l.n_rows, iter =1;
  double tol = 1;
  a_tau = a_tau + p/2;
  b_pi = b_pi + 1;
  //b_tau = a_pi = max_iter=max_tol =0;
  arma::mat Z = create_Z(fir_l, sec_l);
  arma::mat ztz = Z.t() * Z;
  arma::vec len_Z_sq = ztz.diag();
  double len_res_sq = arma::as_scalar(response.t() * response);

  // Initialize
  arma::vec phi_new = arma::randu<arma::vec>(p*q);
  arma::vec e_tau_new = 2.0*arma::randu<arma::vec>(p*q);
  arma::vec sigma_sq_new = arma::randu<arma::vec>(p*q);
  arma::vec mu_new = arma::randu<arma::vec>(p*q);
  arma::vec e_logit_pi_new = arma::randu<arma::vec>(p*q);
  double e_inv_lambda_new = 1.0;

  while((iter<max_iter) && (tol > max_tol)){
    // change value: old <== new
    arma::vec phi_old = phi_new;
    arma::vec e_tau_old = e_tau_new;
    arma::vec sigma_sq_old = sigma_sq_new;
    arma::vec mu_old = mu_new;
    arma::vec e_logit_pi_old = e_logit_pi_new;
    double e_inv_lambda_old = e_inv_lambda_new;

    arma::vec e_b_sq_s1 = sigma_sq_new + square(mu_new);
    arma::vec e_beta = phi_new % mu_new;
    arma::vec var_beta = phi_new % e_b_sq_s1 - square(e_beta);

    //  update tau
    e_tau_new = update_tau(p, q, a_tau, b_tau,
                           phi_old, e_b_sq_s1, e_tau_old);

    // update pi
    e_logit_pi_new = update_pi(a_pi, b_pi, phi_old);


    // update lambda
    e_inv_lambda_new = update_lambda(ztz, e_beta, var_beta,n,len_res_sq);

    // update sigma_sq
    sigma_sq_new = 1/(e_tau_new + e_inv_lambda_new * len_Z_sq);

    // update mu
    mu_new = update_mu(p, q, e_inv_lambda_new,
                       phi_old, mu_old, e_beta, sigma_sq_new, response,
                       Z);

    // update phi
    phi_new = update_phi(e_logit_pi_new, sigma_sq_new, e_tau_new, mu_new);

    // tol and inter
    iter = iter + 1;
    //arma::max
    arma::vec tmp_max = max(abs(phi_new - phi_old), abs(e_logit_pi_old - e_logit_pi_new));
    arma::vec tmp_max2 = max(tmp_max,abs(mu_old - mu_new));
    tmp_max= max(tmp_max2, abs(e_tau_old - e_tau_new));
    tmp_max2 = max(tmp_max, abs(sigma_sq_old - sigma_sq_new));
    tol = max(tmp_max2);
  }

  arma::mat beta_output = phi_new % mu_new;
  beta_output.reshape(q,p);
  beta_output = beta_output.t();
  arma::mat phi_output = phi_new;
  phi_output.reshape(q,p);
  phi_output = phi_output.t();


  return Rcpp::List::create(Rcpp::Named("beta") = beta_output ,
                      Rcpp::Named("phi") = phi_output,
                      Rcpp::Named("lambda") = 1/e_inv_lambda_new,
                      Rcpp::Named("mu") = mu_new,
                      Rcpp::Named("sigma_sq") = sigma_sq_new,
                      Rcpp::Named("tol") = tol,
                      Rcpp::Named("tau") = e_tau_new);

}
// You can include R code blocks in C++ files processed with sourceCpp
// (useful for testing and development). The R code will be autoarma::matically
// run after the compilation.
//
