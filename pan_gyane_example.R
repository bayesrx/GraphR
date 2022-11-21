set.seed(1)
library(GraphR)
library(ghyp)
library(dplyr)
library(reshape2)

load("C:/Users/chenliy/Desktop/PhD/porject I code/code_package/example/pan_gyane.RData")
gene_expr <- apply(gene_expr,2,scale)

res <- GraphR_est(
  gene_expr,
  external,
  a_pi = 1,
  b_pi = 4,
  a_tau = 0.005,
  b_tau = 0.005,
  max_iter = 2000,
  max_tol = 0.001
)

saveRDS(res, "pan_gyane_res.rds")

####### prediction
new_df <- diag(4)
colnames(new_df) <- colnames(external)

pred <- GraphR_pred(new_df, res)
