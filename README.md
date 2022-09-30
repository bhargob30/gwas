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

## Usage in our simulation
The package is applied in performing GWAS for a simulated instance of phenotype using the '*Algorithm $1$*' of the aforementioned article. 
The `R` script `classic_gwas_sim.R` contains the codes for the simulation as well as the whole GWAS procedure. The script corresponds to a particular case of simulated phenotype where we took $100$ individuals, $10000$ bi-allelic SNPs in total out of which $10$ were causal with total genetic heritability=$60 \%$ and SNP heritability=$10 \%$. There were two extra covariates too. For more details on the simulation conditions, see our article

## Special note to users
For faster implementation of large dimensional matrix multiplication we have written the functions of the package using `RcppArmadillo`. But `RcppArmadillo` might be little bit slow on some computers. All the simulation studies we present in our article for classical GWAS were performed in **Google Colab** and **Kaggle** notebooks. And, the computation time we report in our article is based on the implementation of our package in the above two environments. So, if someone want to replicate our results, it is advisable to use the package in the above two platform. We tried to use the package in our personal Windows PC, and the runtime appears to be more. We have, however, developed another package which uses `RcppEigen` instead of `RcppArmadillo`, and it appears to work faster in our personal PC. However, with `RcppEigen` the implementation becomes slower in **Kaggle** and **Colab**. The other package is available on request.
