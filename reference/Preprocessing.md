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

## Value

A modified `raw` object with dye-bias corrected intensities using RELIC
for `correct_dye_bias`. A matrix of beta-values, either normalized (for
`normalize`) or not (for `dont_normalize`).

A matrix of beta-values. If IDATs were imported with option
`on_disk = TRUE`, a {codeff object. If the data was imported using
`on_disk = TRUE`, it is recommended to transpose the returned ff object
for fast row-based access.

## References

Xu Z, Langie SA, De Boever P, Taylor JA, Niu L. RELIC: a novel dye-bias
correction method for Illumina Methylation BeadChip. BMC Genomics. 2017
Jan 3;18(1):4.

Heiss JA, Brenner H. Between-array normalization for 450K data.
Frontiers in Genetics. 2015;6.

## Author

Jonathan Heiss
