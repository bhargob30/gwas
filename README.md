# gwas
## Title
This `R` package contains some functions to calculate the maximum likelihood estimates of the parameters in a genetic LMM
## Description 
The package was written to perform classical GWAS for one-dimensional trait. It contains
some functions written using `RcppArmadillo` for faster implementation of multiple hypothesis
testing in GWAS. The functions are used in calculating restricted maximum likelihood estimates
of the variance parameters as well as MLE of SNP effect size under a typical genetic Linear Mixed Model. 
For more details on the use of the functions, use our pdf `gwas_manual`. The names and notations used in the manual are same as those used in our article 
'**Pool-GWAS methods for complex traits**'. Description of the functions and the GWAS algorithm can be found in the subsection '*Bench-marking of methods for single trait GWAS*' of the article. 
## Installation 
`devtools::install_github("bhargob30/gwas")`

`library(gwas)`
