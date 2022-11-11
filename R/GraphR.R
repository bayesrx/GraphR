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
#' @description Estimate coefficients and inclusion probabilities of external
#' covariates in GraphR models
#' @param features Known as nodes in graphs among which connections are built.
#' A \eqn{n \times p} matrix.
#' @param external External covariates. A \eqn{n \times q} matrix.
#' @param a_pi,b_pi \eqn{\pi \sim} Beta (\eqn{a_\pi,b_\pi})
#' @param a_tau,b_tau \eqn{\tau \sim} Gamma (\eqn{a_\tau,b_\tau})
#' @param max_iter maximum iterations
#' @param max_tol maximum tolerance
#' @return
#' 1. beta: A \eqn{p \times p \times q} array storing coefficients of external
#' covariates. The \eqn{[i,j,k]} elements represents the effect of \eqn{k^{th}}
#' external covariates on regression of \eqn{j^{th}} node on \eqn{i^{th}} node
#'
#' 2. phi: A \eqn{p \times p \times q} array storing posterior inclusion probability (PIP)
#' of external covariates. The \eqn{[i,j,k]} elements represents the PIP of \eqn{k^{th}}
#' external covariates on regression of \eqn{j^{th}} node on \eqn{i^{th}} node
#'
#' 3. omega_diag: A p vector with \eqn{i^{th}} element representing the inverse
#' variance of error



graphR_est <- function(features, external, # input
                   a_pi = 1, b_pi = 4,  #hyperparameter
                   a_tau = 0.005, b_tau = 0.005, #hyperparameter
                   max_iter =2000, max_tol= 0.001 #implementation
                   ){
  p <- ncol(features)
  n <- nrow(features)
  q <- ncol(external)

  beta_list = phi_list = array(0,c(p,p,q))
  # variance: \omega_ii
  lambda_vec <- NULL
  # tolerance
  tol_vec <- NULL
  # variance of b
  tau_vec <- matrix(0,nrow = p,ncol = q)


  for (i in 1:p) {
    response <- features[, i]
    fir_l <- features[, -i]
    tmp <- mfvb(
      response,
      fir_l,
      ex,
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
    phi = phi_list,
    omega_diag = lambda_vec))
}


##### prediction function
#' @title GraphR Model Predictions
#' @description Predict partial correlation between two nodes and the
#' corresponding inclusion probabilities from the results of GraphR model.
#' Bayesian FDR-adjusted p-values are given.
#' @param new_df A matrix of new external covarites based on which predicitons
#' are made
#' @param graphR_est_res Results of `graphR_est`
#' @param beta: A \eqn{p \times p \times q} array storing coefficients of external
#' covariates. The \eqn{[i,j,k]} elements represents the effect of \eqn{k^{th}}
#' external covariates on regression of \eqn{j^{th}} node on \eqn{i^{th}} node
#' @param pip \eqn{p \times p \times q} array storing posterior inclusion probability (PIP)
#' of external covariates. The \eqn{[i,j,k]} elements represents the PIP of \eqn{k^{th}}
#' external covariates on regression of \eqn{j^{th}} node on \eqn{i^{th}} node
#' @param omega_diag: A p vector with \eqn{i^{th}} element representing the inverse
#' variance of error
#' @return
#' 1. feature_id1, feature_id2: Iindices of nodes
#'
#' 2. Pr_inclusion: Posterior inclusion probability of connections between two nodes
#' based on "And" rules.
#'
#' 3. Correlation: Partial correlation between two nodes. Values with maximum magnitudes are provided
#'
#' 4. FDR_p: Bayesian FDR-adjusted p values

graphR_pred <- function(new_df,  ### new external covariates
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



# ##############################
# data <- read.csv("C:/Users/chenliy/Desktop/PhD/porject I code/application1_brca/data/pan_gy.csv")
# ###############
# features <- data[,c(6:10)] %>% as.matrix()
# features <- apply(features,2,scale)
# ex <- data[,c(2:5)] %>% as.matrix()
# res <- graphR_est(features, ex)
# new_df <- ex[c(1:6),]
# new_pred <- graphR_pred(new_df = new_df, graphR_est_res = res)
