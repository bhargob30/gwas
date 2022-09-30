# Classical GWAS using our algorithm

#(1) Install and Load required packages:-

if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
BiocManager::install("snpStats", force = T)
install.packages('PhenotypeSimulator')
install.packages(c('foreach', 'doParallel'))
install.packages("devtools")
library(devtools)
install_github("bhargob30/gwas")

#(2) Generate phenotype based on the simulation conditions given in our article:-

library(PhenotypeSimulator)
phenotype=runSimulation(N=100, P=1, tNrSNP = 10000,
                        cNrSNP = 10, SNPfrequencies=c(0.05,0.1,0.2),
                        pIndependentGenetic = 0,pTraitIndependentGenetic = 0,
                        NrFixedEffects = 2,
                        NrConfounders = c(1,1),pIndependentConfounders=c(0,0),
                        pTraitIndependentConfounders = c(0,0),
                        distConfounders = c("bin","norm"),
                        probConfounders = 0.5, mConfounders=2, genVar = 0.6, h2s=0.1/0.6,
                        theta=1,eta=1, phi=0.6, delta=0.4, alpha=1, gamma=1, seed=4)


# phenotype vector
y=as.matrix(phenotype$phenoComponentsFinal$Y)

#matrix of all genotypes
X=as.matrix(phenotype$rawComponents$genotypes$genotypes)

#kinship matrix
K=as.matrix(phenotype$rawComponents$kinship)

# the matrix of covariates
W=as.matrix(phenotype[["phenoComponentsIntermediate"]][["noiseFixed"]][["cov"]], nrow=nrow(y), ncol=2)

# Run parallelized code to test 10000 hypotheses of association. Result of each hypothesis
# test is a p value. The calculation of test statistic for the hypotheses requires finding 
# RMLE and MLE of the parameters in the genetic LMM. There is a optimization step in between which is 
# done by using `optim` function of R. All the intermediate functions to calculate the MLEs 
# and the test statistic are available in our package 'gwas'. For more details please 
# see algorithm 1 in our article as well as the `gwas_manual.pdf` in our github repository.
# The manual contains detailed description of each of the functions.
library(parallel)
library(foreach)
library(doParallel)
library(gwas)
numCores=detectCores()
registerDoParallel(numCores)

# change upper limit of i, if splitting simulations
pvalues_cov=foreach(i=1:ncol(X), .combine=rbind, .packages ="gwas" ) %dopar% {
  X_augmented=cbind(W, as.vector(X[,i]))
  # change iterlim below based on initial analysis
  rmle=nlm(f=gwas::neglr_lambda_tau, p=c(0,0), K=K, X=X_augmented, y=y, steptol=1e-4, gradtol=1e-4, iterlim=30)
  beta_hat_sigma_hat=gwas::beta_hat_sigma_hat(K=K, X=X_augmented, y=y, log_lambda_tau_hat= rmle$estimate)
  z=beta_hat_sigma_hat[1,1]/beta_hat_sigma_hat[2,1]
  2*min( pnorm(z), pnorm(-z) )
}
stopImplicitCluster()

df=as.data.frame(pvalues_cov)
colnames(df)="p_values"
rownames(df)= phenotype[["setup"]][["id_snps"]] #take subset if splitting simulations
write.csv(df, 'pvalues.csv')