# Quality control metrics

`control_metrics()` calculates the metrics described in the 'BeadArray
Controls Reporter Software Guide' from Illumina. `sample_failure()`
compares these metrics against the thresholds recommended by Illumina.

## Usage

``` r
control_metrics(raw, FFPE = FALSE)

sample_failure(metrics)
```

## Arguments

- raw:

  Output of
  [`read_idats`](https://hhhh5.github.io/ewastools/reference/read_idats.md).

- FFPE:

  Set to `TRUE` if the FFPE DNA restoration kit was used.

- metrics:

  Output of `control_metrics()`.

## Value

For `control_metrics`, a list of 21 control metrics, each one a numeric
vector equal in length to the number of samples. Each vector has an
attribute `"threshold"` holding the cut-off recommended by Illumina
(samples below the threshold fail).

- Restoration:

  Success of the FFPE DNA restoration. Green signal of the restoration
  control relative to background. Only evaluated if `FFPE = TRUE`.

- Staining Green:

  Biotin staining efficiency. Ratio of Biotin (High) to background green
  signals.

- Staining Red:

  DNP staining efficiency. Ratio of DNP (High) to background red
  signals.

- Extension Green:

  Extension efficiency of C/G base pairs. Ratio of green intensities.

- Extension Red:

  Extension efficiency of A/T base pairs. Ratio of red intensities.

- Hybridization High/Medium:

  Hybridization performance. Ratio of high-to-medium green signal
  intensities.

- Hybridization Medium/Low:

  Hybridization performance. Ratio of medium-to-low green signal
  intensities.

- Target Removal 1:

  Target removal efficiency 1. Ratio of green intensities (+3000
  background offset).

- Target Removal 2:

  Target removal efficiency 2. Ratio of green intensities (+3000
  background offset).

- Bisulfite Conversion I Green:

  Bisulfite conversion efficiency for Type I green probes. Ratio of
  green intensities.

- Bisulfite Conversion I Green (Bkg):

  Bisulfite conversion green background ratio (+3000 background offset).

- Bisulfite Conversion I Red:

  Bisulfite conversion efficiency for Type I red probes. Ratio of red
  intensities.

- Bisulfite Conversion I Red (Bkg):

  Bisulfite conversion red background ratio (+3000 background offset).

- Bisulfite Conversion II:

  Bisulfite conversion efficiency for Type II probes. Ratio of
  red-to-green intensities.

- Bisulfite Conversion II (Bkg):

  Bisulfite conversion Type II green background ratio (+3000 background
  offset).

- Specificity I Green:

  Specificity of hybridization for Type I green probes. Ratio of green
  intensities.

- Specificity I Red:

  Specificity of hybridization for Type I red probes. Ratio of red
  intensities.

- Specificity II:

  Specificity of hybridization for Type II probes. Ratio of red-to-green
  intensities.

- Specificity II (Bkg):

  Specificity of hybridization green background ratio (+3000 background
  offset).

- Non-polymorphic Green:

  Non-polymorphic query performance in the green channel.

- Non-polymorphic Red:

  Non-polymorphic query performance in the red channel.

For `sample_failure`, a logical vector, `TRUE` if the sample at the
corresponding index falls below the threshold for any of the control
metrics.

## Note

You can download the Software Guide at
<https://support.illumina.com/downloads/beadarray-controls-reporter-software-guide-1000000004009.html>
