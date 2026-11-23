# Detection p-values

Detection p-values assess whether the signal from a given probe is
significantly different from background noise.

`detectionP()` is the recommended approach that generates realistic
p-values as described in Heiss and Just, 2019. It estimates the
background signal based on the opposite signal at completely methylated
and completely unmethylated probes.

`detectionP.neg()` follows the approach used in GenomeStudio, which
computes p-values based on the intensities of negative control probes
(yielding unrealistically small p-values).

`detectionP.minfi()` provides an implementation for `RGChannelSet`
objects as used in the `minfi` package.

`eval_detP_cutoffs()` generates a plot showing the number of undetected
chrY probes among male and female subjects for various p-value
thresholds, in order to empirically choose a threshold.

Finally, `mask()` masks all observations with detection p-values above a
specified threshold, i.e., sets the corresponding intensities in `M` and
`U` to `NA`.

## Usage

``` r
detectionP(raw)

detectionP.neg(raw)

mask(raw, threshold)

eval_detP_cutoffs(raw, males = NULL, females = NULL)

detectionP.minfi(rgSet)
```

## Arguments

- raw:

  Output of
  [`read_idats`](https://hhhh5.github.io/ewastools/reference/read_idats.md).
  Must include component `detP` for `mask()` and `eval_detP_cutoffs()`.

- threshold:

  p-value threshold above which observations are set to `NA`. Must be on
  the same scale as `raw$detP`, i.e., on the log10 scale if p-values
  were computed with `detectionP.neg()`.

- males, females:

  Column indices of male and female subjects.

- rgSet:

  A `minfi` `RGChannelSet` object.

## Value

For `detectionP()` and `detectionP.neg()`, a modified `raw` list object
with a `detP` component (a matrix of detection p-values) added.
`detectionP()` computes p-values on the linear scale, whereas
`detectionP.neg()` returns p-values on the log10 scale. If the data was
imported with `on_disk = TRUE`, `detP` is an `ff` object.

For `detectionP.minfi()`, a matrix of detection p-values with probe IDs
as row names.

For `mask()`, a modified `raw` list object with undetected observations
set to `NA`.

For `eval_detP_cutoffs()`, Plot of the 90% quantile (males) and 10%
quantile (females) of the number of undetected chrY probes for each
p-value cut-off.

## References

Heiss JA, Just AC. Improved filtering of DNA methylation microarray data
by detection p-values and its impact on downstream analyses. Clinical
Epigenetics (2019) 11:15

## Author

Jonathan Heiss
