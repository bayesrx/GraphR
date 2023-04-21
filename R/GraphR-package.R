#' The 'GraphR' package.
#'
#' @description Graphical Regression (GraphR) is a flexible approach to infer heterogeneous structure
#' of Gaussian graphical models (GGM) which allows for covariate-dependent graphs by
#' incorporating metrics of sample heterogeneity. Our regression-based method
#' provides a functional mapping from the covariate space to precision matrix
#' for different types of heterogeneous graphical model settings. GraphR imposes sparsity
#' in both edge and covariate selection and is computationally efficient via usage
#' of variational Bayes algorithms.
#'
#' The 'GraphR' package contains three main functions: GraphR_est, GraphR_pred and
#' GraphR_visualization. The first two functions allow users to obtain the estimation
#' results of GraphR model and the predicted values given the new external covariates
#' matrix. GraphR_visualization provides the circular network based on
#' a given new external covariates vector and thresholds for FDR-p values and
#' magnitudes of partial correlations.
#'
#' @docType package
#' @name GraphR-package
#' @aliases GraphR
#' @useDynLib GraphR, .registration = TRUE
#' @importFrom Rcpp evalCpp
#' @import dplyr
#' @import reshape2
#' @import igraph
#' @import ggraph
#' @importFrom ghyp Egig
#' @exportPattern '^[[:alpha:]]+'
#' @references
#' Stan Development Team (2022). RStan: the R interface to Stan. R package version 2.21.7. https://mc-stan.org
#'
NULL




