# Dye-bias correction and beta-values

`correct_dye_bias()` corrects the dye bias between the red and green
color channel. It implements an improved version of RELIC using a robust
Theil-Sen estimator to regress the intensities of the red on those of
the green normalization control probes.

`beta_values()` converts the methylated (M) and unmethylated (U)
intensities into beta-values (M / (M + U)). Intensities below 1 are set
to 1 before the conversion.

## Usage

``` r
correct_dye_bias(raw)

beta_values(raw)
```

## Arguments

- raw:

  Output of
  [`read_idats`](https://hhhh5.github.io/ewastools/reference/read_idats.md).

## Value

`correct_dye_bias` returns a modified `raw` list object with dye-bias
corrected intensities in `M`, `U`, and the out-of-band signals.

`beta_values` returns a matrix of beta-values with probe IDs as row
names and sample IDs as column names. If the data was imported with
`on_disk = TRUE`, an `ff` object is returned instead. In this case, it
is recommended to transpose the returned ff object for fast row-based
access.

## References

Xu Z, Langie SA, De Boever P, Taylor JA, Niu L. RELIC: a novel dye-bias
correction method for Illumina Methylation BeadChip. BMC Genomics. 2017
Jan 3;18(1):4.

## Author

Jonathan Heiss
