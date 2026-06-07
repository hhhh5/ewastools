# Detection p-values

Compute detection p-values. p-values are based on the distribution of
the intensities of the negative control probes or the U (M) intensities
observed for completely methylated (unmethylated) probes, respectively
(Heiss and Just, 2019). `detP_threshold` generates a plot\\ showing the
number of undetected Y chromosome probes among male and female subjects
for various p-value thresholds, in order to empirically choose a
threshold. Finally, `mask` is masking an implementation for
`RGChannelSet` objects as used in the `minfi` package.

## Usage

``` r
summits(beta)

detectionP(raw)

detectionP.neg(raw)

mask(raw, threshold)

eval_detP_cutoffs(raw, males = NULL, females = NULL)

detectionP.minfi(rgSet)
```

## Arguments

- raw:

  Output of calling
  [`read_idats`](https://hhhh5.github.io/ewastools/reference/read_idats.md),
  must include component `detP` for `mask` and `detP_threshold`.

- threshold:

  p-value threshold (arithmetic scale) above which oberservations are
  set to NA.

- rgSet:

  minfi rgSet object

- male/female:

  Indices of male and female subjects

## Value

For `detectionP` and `detectionP.neg` a modified `raw` object with a
`detP` component, a matrix of detection p-values, added. `detectionP`
computes p-values on the linear scale, whereas `detectionP.neg` returns
p-values on the log10 scale.

For `detectionP.minfi` a matrix of detection p-values.

For `mask`, a modified `raw` object, with undetected probes set to `NA`.

Vector of length two

## References

Heiss JA, Just AC. Improved filtering of DNA methylation microarray data
by detection p-values and its impact on downstream analyses. Clinical
Epigenetics (2019) 11:15 Return peaks of beta-value distribution

Return the locations of the two peaks (completly methylated and
unmethylated CpG sites) of the beta-value distribution.

## Author

Jonathan A. Heiss
