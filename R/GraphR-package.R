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
#' The 'GraphR' package contains two main functions: GraphR_est and GraphR_pred
#' which provides users with the estimation results of GraphR model and the predicted
#' values given the new matrix.
#'
#' @docType package
#' @name GraphR-package
#' @aliases GraphR
#' @useDynLib GraphR, .registration = TRUE
#' @importFrom Rcpp evalCpp
#' @import dplyr
#' @import reshape2
#' @importFrom ghyp Egig
#' @exportPattern '^[[:alpha:]]+'
#' @references
#' Stan Development Team (2022). RStan: the R interface to Stan. R package version 2.21.7. https://mc-stan.org
#'
NULL




