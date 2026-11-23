# Infer the sex of sample donors

`check_sex()` computes the normalized average total intensities of
probes targeting the X and Y chromosomes. `predict_sex()` uses these
values and the recorded sex of the samples to predict the sex of each
sample, which allows to identify mislabeled samples.

## Usage

``` r
check_sex(raw)

predict_sex(X, Y, male, female)
```

## Arguments

- raw:

  Output of
  [`read_idats`](https://hhhh5.github.io/ewastools/reference/read_idats.md).

- X, Y:

  Normalized X and Y chromosome intensities as returned by
  `check_sex()`.

- male, female:

  Indices of male and female samples.

## Value

`check_sex` returns a list containing two numeric vectors:

- X:

  Normalized average total intensities of probes targeting the X
  chromosome.

- Y:

  Normalized average total intensities of probes targeting the Y
  chromosome.

`predict_sex` returns a factor with levels "m" (male) and "f" (female).
Samples that cannot be determined conclusively are assigned `NA`.

## Details

`check_sex()` should be applied to dye-bias corrected intensities (see
[`correct_dye_bias`](https://hhhh5.github.io/ewastools/reference/Preprocessing.md))
but **before** masking undetected probes with
[`mask`](https://hhhh5.github.io/ewastools/reference/detectionP.md), as
masking would remove many of the Y chromosome probes.

`predict_sex()` determines cutoffs for the X and Y chromosome
intensities as the Hodges-Lehmann estimators of the midpoint between
recorded male and female samples. Samples with low X and high Y
intensities are predicted as male, samples with high X and low Y
intensities are predicted as female.
