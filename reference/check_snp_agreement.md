# Check agreement of genetic fingerprint between samples

These functions check the agreement of genetic fingerprints between
samples (after genotype calling with
[`call_genotypes`](https://hhhh5.github.io/ewastools/reference/call_genotypes.md)).
In case of `check_snp_agreement`, the user supplies for each sample the
assumed donor ID and the function returns all conflicts, meaning either
samples that are supposed to come from the same donor but have distinct
genetic fingerprints, or samples from supposedly different donors with
the same fingerprint. `enumerate_sample_donors` on the other hand infers
donor IDs in case they are unknown (useful, for example, to detect
technical replicates in a public dataset).

## Usage

``` r
check_snp_agreement(genotypes, donor_ids, sample_ids)

enumerate_sample_donors(genotypes)
```

## Arguments

- genotypes:

  Output of
  [`call_genotypes`](https://hhhh5.github.io/ewastools/reference/call_genotypes.md).

- donor_ids:

  Vector of donor IDs corresponding to the samples (columns of the SNP
  matrix). Contains duplicates if several samples come from the same
  donor.

- sample_ids:

  Vector of unique sample IDs corresponding to the samples.

## Value

`check_snp_agreement` returns `NULL` if no conflicts were found.
Otherwise a list of `data.table`s, each representing one connected
component of the conflict graph (i.e., a group of related conflicts),
with columns `donor1`, `sample1`, `donor2`, `sample2`, and `agreement`.

`enumerate_sample_donors` returns an integer vector of inferred donor
IDs, one for each sample. Samples with the same ID share the same
genetic fingerprint.

## Details

Two samples are considered to have the same genetic fingerprint if their
agreement score (the probability of sharing the same genotype, averaged
across SNP probes and weighted by the probability of the observations
not being outliers) exceeds 0.9.

## Author

Jonathan Heiss
