# Quality control metrics

Calculate the QC metrics as described in the 'BeadArray Controls
Reporter Software Guide' from Illumina.

## Usage

``` r
control_metrics(raw)

sample_failure(metrics)
```

## Arguments

- raw:

  Output of calling
  [`read_idats`](https://hhhh5.github.io/ewastools/reference/read_idats.md)

- metrics:

  Output of callling `control_metrics`

## Value

For `control_metrics`, a list of 17 control metrics, each one a numeric
vector equal in length to the sample size

For `sample_failure`, a logical vector, `TRUE` if sample at
corresponding index failed on any of the 17 control metrics.

## Note

You can download the Software Guide at
https://support.illumina.com/downloads/beadarray-controls-reporter-software-guide-1000000004009.html
