# Agreement of SNPs

These function check the agreement of genetic fingerprints between
samples (after genotype calling with `call_genotypes`). In case of
`check_snp_agreement`, the user supplies for each sample the assumed
donor ID and the function returns all conflicts, meaning either samples
that are supposed to come from the same donor but have distinct genetic
fingerprints, or samples from supposedly different donors with the same
fingerprint. `enumerate_sample_donors` on the other hand is inferring
donor IDs in case they are unknown (useful for example to detect
technical replicates in a public dataset).

## Usage

``` r
check_snp_agreement(genotypes, donor_ids, sample_ids)

agreement_(genotypes, donor_ids, sample_ids, ...)

enumerate_sample_donors(genotypes)
```

## Arguments

- genotypes:

  Output of `call_genotypes`

- donor_ids:

  Vector of donor IDs, should include duplicates

- sample_ids:

  Vector of unique sample IDs

## Value

`check_snp_agreement` returns all conflicts found between user-supplied
donor IDs and those derived from the genetic fingerprint. Either `NULL`
if not conflicts were found, or a list of `data.table`s, each
representing one connected component of the conflict graph.

`enumerate_sample_donors` returns for each sample the donor IDs.

## Author

Jonathan A. Heiss
