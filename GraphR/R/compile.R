# library(Rcpp)
# library(RcppArmadillo)
#
# library(ghyp)

update_lambda_r <- function(n, e_zbeta_sq,len_res_sq){
  e_inv_lambda_new = Egig((n+2)/2,e_zbeta_sq,len_res_sq,
                          func = "1/x")
  if (!is.finite(e_inv_lambda_new)){
    e_inv_lambda_new = sqrt(len_res_sq/e_zbeta_sq)*
      ((n+2)/2 + sqrt(((n+2)/2)^2 + e_zbeta_sq*len_res_sq))/(sqrt(e_zbeta_sq*len_res_sq)) -
      (n+2)/e_zbeta_sq
  }
  return(e_inv_lambda_new)
}




#sourceCpp('mfvb_rcpp.cpp')


