
<!-- README.md is generated from README.Rmd. Please edit that file -->

# GraphR overview

<center>
<img src="GraphR_Method_plot_final.png" width="600" height="600" >
</center>

The GraphR (Graphical Regression) is a flexible approach which
incorporates sample heterogenity and enables covariate-dependent graphs.
Our regression-based method provides a functional mapping from the
covariate space to precision matrix for different types of heterogeneous
graphical model settings. GraphR imposes sparsity in both edge and
covariate selection and computationally efficient via use of variational
Bayes algorithms. The method is versatile to incorporate different type
of covariates such as <br> (I) **binary** (control and disease specific
graphs), <br> (II) **categorical** (category specific graphs such as
cancer subtypes), <br> (III) **univariate continuous** (time varying
graphs for single cell data), <br> (IV) **categorical + univariate
continuous** (graphs changing over category such as cancer sub-types and
continuous scale as biomarkers), <br> (V) **multivariate continuous**
(spatial transcriptomics co-expression networks). <br> GraphR is
implemented as an open-source R package.

## Installation

You can install the released version of GraphR from
(<https://github.com/bayesrx/GraphR>) with:

``` r
devtools::install_github("bayesrx/GraphR")
library(GraphR)
```

## Example

An example code with one of the existing datasets to demonstrate how to
run the functions and obtain inference.

The **GraphR\_est** can be used for estimation. The inputs of estimation
of function include features (nodes) and external covariates which are
$n \times p$ and $n \times q$ matrix respectively. Please note that our
function doesn’t provide standardization of features and external
covariates, and thus **please standardize features and external
covariates before plugging into the function**. It is suggested to
maintain $n/pq >1$ and efficacy of the method increase with high values
of $n/pq$ ratio. For priors, we assume $\pi \sim Beta(a_\pi, b_\pi)$ and
$\tau \sim \Gamma(a_\tau, b_\tau)$. **GraphR\_est** function returns a
list containing <br> (1) the coefficient (Beta), <br> (2) posterior
inclusion probability (PIP) of external covaraites and <br> (3) diagonal
elements of precision matrix.

We also provide prediction function **GraphR\_pred** which requires a
new matrix based on which prediction are made and the output of
**GraphR\_est** function. Optionally, one can also enter the the
coefficient (Beta) and posterior inclusion probability (PIP) of external
covaraites, and diagonal elements of precision matrix obtained from
other functions as the input. The output contains following information:
<br> (1) Corresponding external covariates values, <br> (2) Indices of
features, <br> (3) Partial correlations between two features and PIPs of
edges and, <br> (4) Bayesian FDR-adjusted p-values.

``` r
set.seed(100)
library(GraphR)
library(dplyr)
#> 
#> Attaching package: 'dplyr'
#> The following objects are masked from 'package:stats':
#> 
#>     filter, lag
#> The following objects are masked from 'package:base':
#> 
#>     intersect, setdiff, setequal, union
data("Pam50")

features <- apply(Pam50$features,2,scale) %>% as.matrix()
features[c(1:5),c(1:5)]
#>      X1433EPSILON     X4EBP1 X4EBP1_pS65 X4EBP1_pT37T46    X53BP1
#> [1,]   -0.9298711 -1.0325344  -0.1814837      0.3870419 -1.125110
#> [2,]   -1.2265417 -0.8121828  -0.9249897     -0.4834352  1.084052
#> [3,]   -0.9250730 -0.1882466   0.8258123     -0.4022269  0.289943
#> [4,]   -0.6566337 -0.2473042   0.2114522      0.9897723  2.134105
#> [5,]   -0.9476849  1.7654120   2.7128204      2.1739453  1.378139

external <- Pam50$external %>% as.matrix()
external[c(1:5),]
#>      basal_like her2_enriched luminal_ab
#> [1,]          0             1          0
#> [2,]          0             0          1
#> [3,]          0             0          1
#> [4,]          0             0          1
#> [5,]          0             0          1


system.time(res <- GraphR_est(
  features,
  external,
  a_pi = 1,
  b_pi = 4,
  a_tau = 0.005,
  b_tau = 0.005,
  max_iter = 2000,
  max_tol = 0.001
))
#>    user  system elapsed 
#> 233.934   9.091 245.203

####### prediction
new_df <- diag(3)
colnames(new_df) <- colnames(external)

system.time(pred <- GraphR_pred(new_df, res))
#>    user  system elapsed 
#>   1.568   0.058   1.639
head(pred)
#>   basal_like her2_enriched luminal_ab feature_id1 feature_id2 Pr_inclusion
#> 1          1             0          0          10           9            1
#> 2          1             0          0         103         102            1
#> 3          1             0          0         105          54            1
#> 4          1             0          0         113         112            1
#> 5          1             0          0         123          22            1
#> 6          1             0          0         129         128            1
#>   Correlation FDR_p
#> 1   0.6869106     0
#> 2   0.7938597     0
#> 3   0.4530799     0
#> 4   0.8136693     0
#> 5   0.4157183     0
#> 6   0.7698744     0
```

## Shiny App

(link to the Shiny App tab only)

## Tutorial website

(link to the tutorial website)

## Paper

(Link to the paper)

## Supplementary file

(link to the bookdown file)
