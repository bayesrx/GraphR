utils::globalVariables(c("id","phi_min","phi_max","cor_min","cor_max","FDR_p",
                "Correlation","feature1","feature2","cor_sum1","cor_sum2",
                "feature","name","x","y","Var1","Var2","ind_max","ind_min","value"))

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
    filter(Var1 != Var2)
  feature_f <- unique(node_phi_list$Var1)
  node_phi_list$Var1 <- factor(node_phi_list$Var1,
                               levels = feature_f)
  node_phi_list$Var2 <- factor(node_phi_list$Var2,
                               levels = feature_f)

  node_phi_list <- node_phi_list %>%
    mutate(ind_max = pmax(as.numeric(Var1),as.numeric(Var2)),
           ind_min = pmin(as.numeric(Var1),as.numeric(Var2))) %>%
    group_by(ind_max, ind_min) %>%
    summarise(phi_min = min(value),
              phi_max = max(value),
              feature1 = Var1[which.min(as.numeric(Var1))],
              feature2 = Var2[which.max(as.numeric(Var2))],
              .groups = "drop")


  ####### correlation
  cor_list <- apply(cor_list,2, function(x) -x/omega_diag_ind)
  cor_list <- melt(cor_list)  %>%
    filter(Var1 != Var2)

  cor_list$Var1 <- factor(cor_list$Var1,
                               levels = feature_f)
  cor_list$Var2 <- factor(cor_list$Var2,
                               levels = feature_f)

  cor_list <- cor_list %>%
    mutate(ind_max = pmax(as.numeric(Var1),as.numeric(Var2)),
           ind_min = pmin(as.numeric(Var1),as.numeric(Var2))) %>%
    group_by(ind_max, ind_min) %>%
    summarise(cor_max = abs_max_fcn(value),
              cor_min = abs_min_fcn(value),
              feature1 = Var1[which.min(as.numeric(Var1))],
              feature2 = Var2[which.max(as.numeric(Var2))],
              .groups = "drop")

  cor_phi <- merge(node_phi_list,cor_list)
  cor_phi <- mutate(cor_phi, id = id)
  cor_phi <- select(cor_phi,-ind_max,-ind_min)
  return(cor_phi)

}


#' @title wrap function
#' @param n sample size
#' @param e_zbeta_sq expectation
#' @param len_res_sq length
#' @noRd
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
#' @param max_tol Maximum tolerance. Default as 0.001.
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
#' library(magrittr)
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
#'   max_iter = 5,
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
  features_colid <- colnames(features) ## gene_id
  if (is.null(features_colid)){
    features_colid <- paste0("feature",c(1:p),recycle0 = TRUE)
  }

  if (standardize_feature){
    features <- apply(features,2,scale)
  }
  features <- as.matrix(features)


  #external_colid <- c(colnames(cont_external),colnames(dis_external))
  if ((is.null(cont_external)) & (is.null(dis_external))){
    stop("External covariates are required")
  }

  if (standardize_external){
    if (is.null(cont_external)) {
      stop("Continuous external covariates are required")
    } else{
      cont_external <- apply(cont_external,2,scale) %>% as.matrix()
    }
  }

  external <- cbind(cont_external,dis_external)
  q <- ncol(external)
  external_colid <- c(colnames(cont_external),colnames(dis_external))
  # if (!is.null(cont_external)){
  #   if (standardize_external){
  #     cont_external <- apply(cont_external,2,scale) %>% as.matrix()
  #   }
  #   if (!is.null(dis_external)){
  #     external <- cbind(cont_external,dis_external)
  #   } else{
  #     external <- cont_external
  #   }
  # } else{
  #   if (!is.null(dis_external)){
  #     external <- dis_external
  #   } else{
  #     stop("External covariates are required")
  #   }
  # }

  if (is.null(external_colid)){
    external_colid <- paste0("external",c(1:q),recycle0 = TRUE)
  }
  external <- as.matrix(external)


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
  dimnames(beta_list)[[1]] = dimnames(beta_list)[[2]] = dimnames(phi_list)[[1]] =
    dimnames(phi_list)[[2]] = features_colid
  dimnames(beta_list)[[3]] = dimnames(phi_list)[[3]] = external_colid
  names(lambda_vec) <- features_colid

  return(list(
    beta = beta_list,
    phi = phi_list,
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
#' @param phi A p \eqn{x} p \eqn{x} q array storing posterior inclusion probability (PIP)
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
#' library(magrittr)
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
#'   max_iter = 5,
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
  new_df <- cbind(c(1:n), new_df) %>% as.data.frame()
  colnames(new_df)[1] <- "id"
  new_df <- as.matrix(new_df)
  ####function(external_ind, beta_ind,phi_ind, omega_diag_ind)
  cor_phi_list <- apply(new_df,1,function(x) pred_ind(x,beta,phi,omega_diag))
  cor_phi_list2 <- Reduce(rbind, cor_phi_list)
  cor_phi_list2 <- full_join(as.data.frame(new_df),
                             cor_phi_list2, by ="id",multiple = "all") %>%
    select(-id)
  colnames(cor_phi_list2)[c(1:q)] <- dimnames(beta)[[3]]


  cor_phi_list2 <- mutate(cor_phi_list2, q = 1-phi_min) %>%
    arrange(q) %>%
    mutate(fdr_p = cummean(q)) %>%
    select(-q,-phi_max,-cor_min)

  # write.csv(cor_phi_list, "location_scale_res_all.csv")
  cor_phi_list2$cor_max <- ifelse(cor_phi_list2$cor_max > 1,1,cor_phi_list2$cor_max)
  cor_phi_list2$cor_max <- ifelse(cor_phi_list2$cor_max < -1,-1,cor_phi_list2$cor_max)

  colnames(cor_phi_list2)[c((q+3):(q+5))] <- c("Pr_inclusion","Correlation",
                                               "FDR_p")

  return(cor_phi_list2)

}



##### Network plotting function
#' @title Network Visualization.
#' @description Plot of networks based on the given new external covariates vector
#' and thresholds for FDR-p values and magnitudes of partial correlation.
#' @param new_vec A vector of new external covarites based on which plots
#' are made.
#'
#' Note: Please ensure that the order and scale of new external covariates are same as those used in the estimation.
#'
#' @param graphR_est_res Results from `GraphR_est` function.
#' If graphR_est_res = NULL, then beta, phi and omega_diag are needed simultaneously.
#' @param beta A p \eqn{x} p \eqn{x} q array storing coefficients of external
#' covariates. The \eqn{[i,j,k]} elements represents the effect of k-th
#' external covariates on regression of j-th node on i-th node.
#' @param phi A p \eqn{x} p \eqn{x} q array storing posterior inclusion probability (PIP)
#' of external covariates. The \eqn{[i,j,k]} elements represents the PIP of k-th
#' external covariates on regression of j-th node on i-th node.
#' @param omega_diag A p vector with i-th element representing the inverse
#' variance of error.
#' @param fdr_thre FDR-controlled significance level. Default as 0.01
#' @param magnitude_thre Threshold for magnitude of selected partial correlations. Default as 0.4
#' @return
#' \item{Plot}{Plot of circular networks. Node sizes represent connectivity degrees
#' of the corresponding features while edge widths are proportional to the partial
#' correlation between two features. Sign of the partial correlations are represented
#' by the color}
#'
#' @examples
#' data("Pam50")
#' library(magrittr)
#' features <- apply(Pam50$features,2,scale) %>% as.matrix()
#' features[c(1:5),c(1:5)]
#'
#' external <- Pam50$external %>% as.matrix()
#' external[c(1:5),]
#'
#'
#' res <- GraphR_est(
#'   features,
#'   external,
#'   a_pi = 1,
#'   b_pi = 4,
#'   a_tau = 0.005,
#'   b_tau = 0.005,
#'   max_iter = 5,
#'   max_tol = 0.001
#' )
#'
#' new_vec <- c(1,0,0)
#'
#'
#' GraphR_visualization(new_vec, graphR_est_res = res)
#'
#'
#'
#'
#'
GraphR_visualization <- function(new_vec,  ### new external covariates
                        graphR_est_res = NULL,  ### results obtained from graphR_est
                        beta = NULL, phi = NULL, omega_diag = NULL,
                        fdr_thre = 0.01, magnitude_thre = 0.4){

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


  ##### Check dimension of new dataframe
  new_vec <- matrix(new_vec, nrow = 1)

  cor_phi_list <- GraphR_pred(new_vec, beta = beta,
                              phi = phi, omega_diag = omega_diag)

  cor_phi_list <- cor_phi_list %>%
    filter(FDR_p <= fdr_thre) %>%
    filter(abs(Correlation) >= magnitude_thre)

  sum_f1 <- cor_phi_list %>% group_by(feature1) %>%
    summarise(cor_sum1 = sum(abs(Correlation)))
  sum_f2 <- cor_phi_list %>% group_by(feature2) %>%
    summarise(cor_sum2 = sum(abs(Correlation)))
  colnames(sum_f1)[1] = colnames(sum_f2)[1] = "feature"
  sum_all <- full_join(sum_f1,sum_f2)
  sum_all[is.na(sum_all)] <- 0
  sum_all <- mutate(sum_all, node_size = cor_sum1 + cor_sum2) %>%
    select(feature,node_size) %>%
    mutate(id = 1:nrow(sum_all))

  cor_phi_list <- select(cor_phi_list, "feature1","feature2","Correlation")
  colnames(cor_phi_list) <- c("from","to","Correlation")

  g <- graph_from_data_frame(cor_phi_list, directed=FALSE, vertices=sum_all)
  node_range = c(min(sum_all$node_size),
                      max(sum_all$node_size))
  edge_abs_range = c(min(abs(cor_phi_list$Correlation)),
                          max(abs(cor_phi_list$Correlation)))
  node_size <- sum_all$node_size
  p3 <- ggraph(g, layout = "circle") +
    #geom_edge_link() +
    ggraph::geom_edge_arc(aes(color = factor(sign(Correlation),
                                     levels = c(1,-1)),
                      width = abs(Correlation)),
                  strength = 0.2,
                  alpha = 0.8)+
    scale_edge_width(name = "Partial Correlation Magnitude",
                     breaks = c(0.2,0.4,0.6,0.8),
                     range = 1*edge_abs_range)+
    scale_edge_color_discrete(name = "Partial Correlation Sign", labels = c("Positive", "Negative")) +
    #geom_edge_link(aes(edge_width = width,color = as.factor(sign)))+
    geom_node_point(aes(size=node_size),color = "darkgrey") +
    scale_size(name = "Connectivity Degree",
               breaks = c(0.5,1,1.5,2),
               range = 2*node_range) +
    theme_void() +
    #guides(size = "none")+
    #theme(legend.position = "none") +
    #ggtitle("brca-like") +
    #theme(plot.title = element_text(face = "bold",size = 25,hjust = 0.5)) +
    coord_fixed() +
    geom_node_text(aes(label=name,
                       x = x * 1.05, y = y* 1.05,
                       angle = ifelse(atan(-(x/y))*(180/pi) < 0,
                                      90 + atan(-(x/y))*(180/pi),
                                      270 + atan(-x/y)*(180/pi)),
                       hjust = ifelse(x > 0, 0 ,1)),
                   size=2.5, alpha=10, fontface = "bold")+
    expand_limits(x = c(-2, 2), y = c(-2, 2))+
    guides(size = guide_legend(order = 1))+
    #edge_width = guide_legend(order = 2),
    #edge_color = guide_legend(order = 3)) +
    theme(
      #legend.title = element_text(size=4,face = "bold"),
      legend.title=element_text(size = 6,face = "bold"),
      legend.text = element_text(size = 5),
      legend.position = "bottom",
      legend.box="vertical",
      legend.margin = margin(0, 0, -0.5, -0.5),
      legend.spacing.x = unit(0, "mm"),
      legend.spacing.y = unit(0, "mm"))

  return(p3)
}


