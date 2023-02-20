# source("compile.R")
#
# library(dplyr)
# library(reshape2)
######## pre_function for prediction
abs_max_fcn <- function(x) {
  return(x[which(abs(x) == max(abs(x)))])
}

abs_min_fcn <- function(x) {
  return(x[which(abs(x) == min(abs(x)))])
}


#### function for each individuals
pred_ind <- function(external_ind,
                     beta_ind,phi_ind, omega_diag_ind){

  id <- external_ind[1]
  external_ind <- external_ind[-1]

  p_ind = dim(beta_ind)[1]
  q_ind = dim(beta_ind)[3]

  cor_list <- matrix(0,nrow = p_ind,ncol = p_ind)
  node_phi_list <- matrix(0,nrow = p_ind,ncol = p_ind)
  x_beta_list <- array(0,c(p_ind,p_ind,q_ind))

  for (i in 1:q_ind){
    cor_list = cor_list + external_ind[i] * beta_ind[,,i]
    x_beta_list[,,i] <- external_ind[i] * beta_ind[,,i]
  }


  ####### phi
  x_beta_list_sq <- x_beta_list^2
  x_beta_list_sum <- apply(x_beta_list_sq, c(1,2), sum)
  x_beta_list_pro <- x_beta_list_sq

  #### need help here
  for (q_this in 1:q_ind){
    x_beta_list_pro[,,q_this] <-
      x_beta_list_pro[,,q_this] / x_beta_list_sum
  }
  for (q_this in 1:q_ind){
    node_phi_list <- node_phi_list +
      x_beta_list_pro[,,q_this] * phi_ind[,,q_this]
  }

  node_phi_list <- melt(node_phi_list) %>%
    filter(Var1 != Var2) %>%
    mutate(ind_max = pmax(Var1,Var2),
           ind_min = pmin(Var1,Var2)) %>%
    select(ind_max, ind_min, value) %>%
    group_by(ind_max, ind_min) %>%
    summarise(phi_min = min(value),
              phi_max = max(value),
              .groups = "drop")


  ####### correlation
  cor_list <- apply(cor_list,2, function(x) -x/omega_diag_ind)
  cor_list <- melt(cor_list)
  cor_list <- cor_list %>% mutate(ind_max = pmax(Var1, Var2),
                                  ind_min = pmin(Var1, Var2)) %>%
    select(ind_max,ind_min,value) %>%
    filter(ind_max != ind_min) %>%
    group_by(ind_max,ind_min) %>%
    summarise(cor_max = abs_max_fcn(value),
              cor_min = abs_min_fcn(value),
              .groups = "drop")

  cor_phi <- merge(node_phi_list,cor_list)
  cor_phi <- mutate(cor_phi, id = id)
  return(cor_phi)

}









##### estimation function
#' @title GraphR Model Estimation
#' @description Estimate the graphical regression coefficients and inclusion probabilities of external
#' covariates for the GraphR models.
#' @param features Nodes of the graphs among which edges are built (e.g. a gene expression matrix of dimensions n \eqn{\times} p).
#' @param cont_external,dis_external Continuous and discrete external covariates (n \eqn{\times} \eqn{\text{q}_1} and n \eqn{\times} \eqn{\text{q}_2} matrix, with \eqn{\text{q}_1} + \eqn{\text{q}_2} = q).
#' @param a_pi,b_pi \eqn{\pi} ~ Beta(\eqn{a_\pi, b_\pi}). By default \eqn{a_\pi} = 1, \eqn{b_\pi} = 4.
#' @param a_tau,b_tau \eqn{\tau} ~ Gamma(\eqn{a_\tau, b_\tau}). By default \eqn{a_\tau} = 0.005, \eqn{b_\tau} = 0.005.
#' @param standardize_feature Standardize features. Default as FALSE
#' @param standardize_external Standardize continuous external covariates. Default as FALSE
#' @param max_iter Maximum iterations. Default as 2,000.
#' @param max_tol Maximum tolerance. Default as 0.01.
#' @return
#' \item{beta}{A p \eqn{\times} p \eqn{\times} q array of coefficients for external
#' covariates. The \eqn{[i,j,k]} element represents the effect of k-th
#' external covariates on regression of j-th node on i-th node.}
#' \item{phi}{A p \eqn{\times} p \eqn{\times} q array storing posterior inclusion probability (PIP)
#' of external covariates. The \eqn{[i,j,k]} elements represents the PIP of k-th
#' external covariates on regression of j-th node on i-th node.}
#' \item{omega_diag}{A p vector with i-th element representing the inverse
#' variance of error.}
#' @examples
#' set.seed(100)
#' data("Pam50")
#'
#' features <- apply(Pam50$features,2,scale) %>% as.matrix()
#' dim(features)
#' features[c(1:5),c(1:5)]
#'
#' external <- Pam50$external %>% as.matrix()
#' dim(external)
#' external[c(1:5),]
#'
#'
#' system.time(res <- GraphR_est(
#'   features,
#'   dis_external = external,
#'   a_pi = 1,
#'   b_pi = 4,
#'   a_tau = 0.005,
#'   b_tau = 0.005,
#'   standardize_feature = FALSE,
#'   standardize_external =FALSE,
#'   max_iter = 2000,
#'   max_tol = 0.001
#' ))
#'



GraphR_est <- function(features, cont_external = NULL, dis_external = NULL, # input
                   a_pi = 1, b_pi = 4,  #hyperparameter
                   a_tau = 0.005, b_tau = 0.005, #hyperparameter
                   standardize_feature = FALSE, standardize_external =FALSE, #standardization , default as FALSE
                   max_iter =2000, max_tol= 0.001 #implementation
                   ){
  p <- ncol(features)
  n <- nrow(features)

  ext_tmp <- cbind(cont_external, dis_external)
  q <- ncol(ext_tmp)

  if (n <= p*q){
    warning("Low n/pq ratio may lead to inaccurate estimation")
  }
  beta_list = phi_list = array(0,c(p,p,q))
  # variance: \omega_ii
  lambda_vec <- NULL
  # tolerance
  tol_vec <- NULL
  # variance of b
  tau_vec <- matrix(0,nrow = p,ncol = q)

  if (standardize_feature){
    features <- apply(features,2,scale) %>% as.matrix()
  }
  if (standardize_external){
    cont_external <- apply(cont_external,2,scale) %>% as.matrix()
    external <- cbind(cont_external, dis_external)
  } else{
    external <- cbind(cont_external, dis_external) %>% as.matrix()
  }


  for (i in 1:p) {
    response <- features[, i]
    fir_l <- features[, -i]
    tmp <- mfvb(
      response,
      fir_l,
      external,
      max_tol = max_tol,
      max_iter = max_iter,
      a_tau = a_tau,
      b_tau = b_tau,
      a_pi=a_pi, b_pi=b_pi
    )

    for (k in 1:q){
      phi_this <- append(tmp$phi[, k], 0, i-1)
      phi_list[i,,k] <- phi_this

      beta_this <- append(tmp$beta[, k], 0, i-1)
      beta_list[i,,k] <- beta_this

      tau_vec[i,k] <- tmp$tau[k]
    }

    lambda_vec <- c(lambda_vec, tmp$lambda)
    tol_vec <- c(tol_vec, tmp$tol)
  }

  return(list(
    beta = beta_list,
    pip = phi_list,
    omega_diag = lambda_vec))
}


##### prediction function
#' @title GraphR Model Predictions
#' @description Prediction of partial correlation between two nodes and the
#' corresponding inclusion probabilities from the results of GraphR model along with
#' Bayesian FDR-adjusted p-values.
#' @param new_df A matrix of new external covarites based on which predictions
#' are made.
#'
#' Note: Please ensure that the order and scale of new external covariates are same as those used in the estimation.
#'
#' @param graphR_est_res Results from `GraphR_est` function.
#' If graphR_est_res = NULL, then beta, phi and omega_diag are needed simultaneously.
#' @param beta A p \eqn{x} p \eqn{x} q array storing coefficients of external
#' covariates. The \eqn{[i,j,k]} elements represents the effect of k-th
#' external covariates on regression of j-th node on i-th node.
#' @param pip A p \eqn{x} p \eqn{x} q array storing posterior inclusion probability (PIP)
#' of external covariates. The \eqn{[i,j,k]} elements represents the PIP of k-th
#' external covariates on regression of j-th node on i-th node.
#' @param omega_diag A p vector with i-th element representing the inverse
#' variance of error.
#' @return
#' \item{feature_id1, feature_id2}{Indices of nodes.}
#'
#' \item{Pr_inclusion}{Posterior inclusion probability of connections between two nodes
#' based on "And" rules.}
#'
#' \item{Correlation}{Partial correlation between two nodes. Values with maximum magnitudes are provided.}
#'
#' \item{FDR_p}{Bayesian FDR-adjusted p values.}
#' @examples
#' set.seed(100)
#' data("Pam50")
#'
#' features <- apply(Pam50$features,2,scale) %>% as.matrix()
#' features[c(1:5),c(1:5)]
#'
#' external <- Pam50$external %>% as.matrix()
#' external[c(1:5),]
#'
#'
#' system.time(res <- GraphR_est(
#'   features,
#'   external,
#'   a_pi = 1,
#'   b_pi = 4,
#'   a_tau = 0.005,
#'   b_tau = 0.005,
#'   max_iter = 2000,
#'   max_tol = 0.001
#' ))
#'
#' ####### prediction
#' new_df <- diag(3)
#' colnames(new_df) <- colnames(external)
#'
#' pred <- GraphR_pred(new_df, res)
#' head(pred)
#'
#'
#'
#'
#'

GraphR_pred <- function(new_df,  ### new external covariates
                        graphR_est_res = NULL,  ### results obtained from graphR_est
                        beta = NULL, phi = NULL, omega_diag = NULL){
  #graphR_est_res <- res
  if (is.null(graphR_est_res) ==FALSE){
    beta = graphR_est_res[[1]]
    phi = graphR_est_res[[2]]
    omega_diag = graphR_est_res[[3]]
  }
  ##### Missing
  try(if(is.null(beta) | is.null(phi) | is.null(omega_diag))
      stop("Missing in parameters"))

  #### Dimension of parameters
  try(
    if(length(dim(beta)) == 3 &
       length(dim(phi)) == 3 &
       dim(beta)[1] == dim(phi)[1] &
       dim(beta)[2] == dim(phi)[2] &
       dim(beta)[3] == dim(phi)[3] &
       dim(beta)[1] == dim(beta)[2] &
       dim(beta)[1] == length(omega_diag)){
      p = dim(beta)[1]
      q = dim(beta)[3]
    } else{
      stop("Errors in parameters dimension")
    })

  # p = dim(beta)[1]
  # q = dim(beta)[3]

  ##### Check dimension of new dataframe
  try(if(ncol(new_df) != q) stop("Errors in new_df dimension"))

  #### add id
  n <- nrow(new_df)
  new_df <- cbind(c(1:n), new_df)
  colnames(new_df)[1] <- "id"
  ####function(external_ind, beta_ind,phi_ind, omega_diag_ind)

  cor_phi_list <- apply(new_df,1,function(x) pred_ind(x,beta,phi,omega_diag))
  cor_phi_list2 <- Reduce(rbind, cor_phi_list)
  cor_phi_list2 <- full_join(as.data.frame(new_df),
                             cor_phi_list2, by ="id") %>%
    select(-id)

  cor_phi_list2 <- mutate(cor_phi_list2, q = 1-phi_min) %>%
    arrange(q) %>%
    mutate(fdr_p = cummean(q)) %>%
    select(-q,-phi_max,-cor_min)

  # write.csv(cor_phi_list, "location_scale_res_all.csv")
  cor_phi_list2$cor_max <- ifelse(cor_phi_list2$cor_max > 1,1,cor_phi_list2$cor_max)
  cor_phi_list2$cor_max <- ifelse(cor_phi_list2$cor_max < -1,-1,cor_phi_list2$cor_max)

  colnames(cor_phi_list2)[c((q+1):(q+5))] <- c("feature_id1","feature_id2",
                                               "Pr_inclusion","Correlation",
                                               "FDR_p")

  return(cor_phi_list2)

}








