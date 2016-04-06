MixRF
=====

### A random-forest-based approach for imputing clustered incomplete data

This `R` package offers random-forest-based functions to impute clustered incomplete data. The package is tailored for but not limited to imputing multitissue expression data, in which a gene's expression is measured on the collected tissues of an individual but missing on the uncollected tissues.

### Installation
- For the stable version from [CRAN](https://cran.r-project.org/web/packages/MixRF/index.html):
```r
install.packages('MixRF')
```
- For the development version (requiring the `devtools` package):
```r
devtools::install_github('randel/MixRF')
```

### Reference

Wang, J., Gamazon, E. R., Pierce, B. L., Stranger B. E., Im, H. K., Gibbons, R. D., Cox, N. J., Nicolae, D. L., & Chen, L. S. (2016). Imputing Gene Expression in Uncollected Tissues Within and Beyond GTEx. *American Journal of Human Genetics*. http://dx.doi.org/10.1016/j.ajhg.2016.02.020
