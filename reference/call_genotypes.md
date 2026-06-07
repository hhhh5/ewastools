# Genotype calling

Detect SNP probes which do not fit into on of the three categories
(AA,AB,BB). A mixture model (3 Beta distributions, 1 uniform
distribution for outliers) is fitted to all SNP probes. After learning
the model parameters via EM algorithm, the probability of being an
outlier is computed for each SNP.

## Usage

``` r
call_genotypes(snpmatrix, learn = FALSE, maxiter = 50)

mxm_(genotypes)

snp_outliers(genotypes)

eBeta(x, w)
```

## Arguments

- snpmatrix:

  Matrix of beta-values for SNP probes. Provide SNPs probes as rows and
  samples as columns.

- maxiter:

  Maximal number of iterations of the Expectation-Maximization algorithm
  learning the mixture model

- genotypes:

  Output of `call_genotypes`

## Value

For `call_genotypes`, a list containing

- par:

  Parameters of the mixture model

- loglik:

  Log-likelihood in each iteration of the EM algorithm

- outliers:

  A-posteriori probability of SNP being an outlier

- gamma:

  A-posteriori probabilities for each of the three genotypes

For `snp_outliers`, a metric assessing the outlierness of the SNP
beta-values. High values may indicate either contaminated or failed
samples.

For `mxm_`, a histogram showing the distribution of beta-values for SNP
probes with the density function of the mixture model overlaid.

## Author

Jonathan A. Heiss
