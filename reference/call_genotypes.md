# Genotype calling

`call_genotypes()` calls genotypes from the SNP probes on the array and
detects outlying observations. A mixture model (three Beta distributions
representing the genotypes AA, AB and BB, and one uniform distribution
representing outliers) is fitted to the beta-values of all SNP probes.
Either predefined model parameters are used (default), or the parameters
are learned from the data via the EM algorithm. Subsequently, the
a-posteriori probabilities of each genotype and of being an outlier are
computed for every SNP observation.

`snp_outliers()` summarizes the outlier probabilities per sample as the
average log2 odds of being an outlier across all SNP probes. Samples
with a score greater than -4 should be flagged as they are likely
contaminated or failed assays.

## Usage

``` r
call_genotypes(snpmatrix, learn = FALSE, maxiter = 50)

snp_outliers(genotypes)
```

## Arguments

- snpmatrix:

  Matrix of beta-values for SNP probes. Provide SNP probes as rows and
  samples as columns.

- learn:

  Logical. If `TRUE`, model parameters are learned from the dataset
  using the EM algorithm. If `FALSE` (default), predefined model
  parameters are used. The predefined parameters are recommended for
  small datasets or datasets containing many outliers.

- maxiter:

  Maximum number of iterations of the EM algorithm. Only used if
  `learn = TRUE`.

- genotypes:

  Output of `call_genotypes()`.

## Value

For `call_genotypes()`, a list containing:

- snps:

  The input `snpmatrix`.

- outliers:

  A matrix of a-posteriori probabilities of each SNP observation being
  an outlier.

- gamma:

  A list of three matrices containing the a-posteriori probabilities for
  each of the three genotypes (AA, AB, BB) respectively.

- par:

  Parameters of the mixture model (class priors `pi`, Beta distribution
  parameters `shapes1` and `shapes2`, and outlier prior `alpha`).

- loglik:

  Log-likelihood in each iteration of the EM algorithm (`NULL` if
  `learn = FALSE`).

For `snp_outliers()`, a numeric vector with one outlier score per sample
(average log2 odds of SNP observations being outliers). High values
indicate contaminated or failed samples.

## See also

[`check_snp_agreement`](https://hhhh5.github.io/ewastools/reference/check_snp_agreement.md)
and
[`enumerate_sample_donors`](https://hhhh5.github.io/ewastools/reference/check_snp_agreement.md)
for using the called genotypes to track samples.

## Author

Jonathan Heiss
