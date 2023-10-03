
<!-- README.md is generated from README.Rmd. Please edit that file -->

# STew <img width="25%" align = "right" src="https://github.com/fanzhanglab/STew/blob/main/STewlogo.png">

[![R-CMD-check](https://github.com/fanzhanglab/STew/actions/workflows/check-standard.yaml/badge.svg)](https://github.com/fanzhanglab/STew/actions/workflows/check-standard.yaml)
![](https://komarev.com/ghpvc/?username=fanzhanglab&style=flat-square&color=green)

<!-- badges: start -->
<!-- badges: end -->
<p align="justify">
We introduce STew, a multi-view representation learning method for
spatial transcriptomic data, to jointly characterize the gene expression
variation and spatial information in the shared low-dimenion space in a
scalable manner. STew will output distinct spatially informed cell
gradients, robust clusters, and statistical goodness of model fit to
reveal significant genes that reflect subtle spatial niches in complex
tissues.
</p>

<img width="100%" align = "center" src="https://github.com/fanzhanglab/STew/blob/main/man/figures/Figure1.png">

</br>

## Installation

You can install the STew Package from
[GitHub](https://github.com/fanzhanglab/STew/) using the devtools as
follows:

``` r
# install.packages("devtools")
devtools::install_github("fanzhanglab/STew")
```

(OR)

``` r
remotes::install_github("fanzhanglab/STew")
```

<br/>

### Dependencies / Other required packages

- R (\>= 4.2)
- MERINGUE (\>= 1.0)
- loe (\>= 1.1)
- Matrix (\>= 1.5.4)
- ggplot2 (\>= 3.4.2)
- ggpubr (\>= 0.6.0)
- Seurat (\>= 4.3.0)
- future.apply (\>= 1.10.0)
- RANN (\>= 2.6.1)
- sctransform
- tibble

<br/>

## Tutorials

**Step-by-step notebook** of applying STew on identifying spatially
informed low-dimensional embeddings and spatially aware clusters on the
10X Visium Human Brain Data (DLPFC):

- <a href="https://htmlpreview.github.io/?https://github.com/fanzhanglab/STew/blob/main/vignettes/Tutorial_STew_DLPFC.html">
  Tutorial of applying STew on DLPFC data </a>
- <a href = 'https://htmlpreview.github.io/?https://github.com/fanzhanglab/STew/blob/main/vignettes/count_modeling_tutorial_dlfcp.nb.html'>
  Tutorial of count data modelling </a>

<br/>

#### Below are several major steps of running STew:

``` r
# Create a new STew object for the loaded spatial transcriptomic data
STew = STew_Obj(count = dlpfc$count_exp,
                  spatial = dlpfc$spatial)
```

… (skip several preprocessing steps) …

``` r
# permute optimal penalty parameters
STew <- parallel_cca_permute(x = STew$exp_adj_matrix, z = STew$adj_matrix, obj = STew, nperms=50, niter=3)
```

``` r
# Perform sparse CCA based on the optimal penalty parameters
STew <- cca_main(x = STew$exp_adj_matrix, z = STew$adj_matrix, obj = STew, K=20, penaltyx=STew$bestpenaltyx, penaltyz=STew$bestpenaltyz, v=STew$v.init)
```

``` r
gradient_plot <- spatial_gradient(STew)
gradient_plot[1:5]
```

<img width="85%" align = "center" src="https://github.com/fanzhanglab/STew/blob/main/man/figures/cca_vis.png">

<br/>

``` r
cluster_plot <- plot_cluster(coordis = spatial, label = cluster$res_0.30, colors = colors, t="Cell clusters based on STew")
cluster_plot
```

<img width="25%" align = "center" src="https://github.com/fanzhanglab/STew/blob/main/man/figures/README-unnamed-chunk-14-1.png">

<br/>

``` r
# Save the main results into the STew object
saveRDS(STew, file="STew_10x_human_dlpfc_no.rds")
```

<br/>

#### Benchmarcking STew with other algorithms:

<img width="85%" align = "center" src="https://github.com/fanzhanglab/STew/blob/main/man/figures/ARI_DLPFC.png">

<br/>

## Citations

Guo, N., Vargas, J., Fritz, D., Krishna, R., Zhang, F. Uncover spatially
informed shared variations underlying single-cell spatial
transcriptomics with STew, [*bioRxiv*](link), 2023

<br/>

## Help, Suggestion and Contribution

Using github [**issues**](https://github.com/fanzhanglab/STew/issues)
section, if you have any question, comments, suggestions, or to report
coding related issues of STew is highly encouranged than sending emails.

- Please **check the GitHub
  [issues](https://github.com/fanzhanglab/STew/issues)** for similar
  issues that has been reported and resolved. This helps the team to
  focus on adding new features and working on cool projects instead of
  resolving the same issues!
- **Examples** are required when filing a GitHub issue. In certain
  cases, please share your STew object and related codes to understand
  the issues.

<br/>

## Contact

Please contact [fanzhanglab@gmail.com](fanzhanglab@gmail.com) for
further questions or protential collaborative opportunities!
