# Functions for preprocessing

Functions for preprocessing

Convert fluorescence intensities to beta-values

## Usage

``` r
correct_dye_bias(raw)

correct_dye_bias2(raw)

beta_values(raw)
```

## Arguments

- raw:

  Output of calling
  [`read_idats`](https://hhhh5.github.io/ewastools/reference/read_idats.md)

- tissue:

  Optional. If set to "blood", predefined reference values are used for
  normalization. Recommened when cell proportions are to be estimated.

## Value

A modified `raw` object with dye-bias corrected intensities using RELIC
for `correct_dye_bias`. A matrix of beta-values, either normalized (for
`normalize`) or not (for `dont_normalize`).

A matrix of beta-values. If IDATs were imported with option
`on_disk = TRUE`, a {codeff object.

## References

Xu Z, Langie SA, De Boever P, Taylor JA, Niu L. RELIC: a novel dye-bias
correction method for Illumina Methylation BeadChip. BMC Genomics. 2017
Jan 3;18(1):4.

Heiss JA, Brenner H. Between-array normalization for 450K data.
Frontiers in Genetics. 2015;6.

## Author

Jonathan Heiss
