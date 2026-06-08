# Infer the sex of sample donor

Return the normalized average total intensities of probes targeting the
X and Y chromosomes.

## Usage

``` r
check_sex(raw)

predict_sex(X, Y, male, female)
```

## Arguments

- raw:

  Output of calling
  [`read_idats()`](https://hhhh5.github.io/ewastools/reference/read_idats.md)

- X, Y:

  Forwarded from `check_sex()`

- male, female:

  Indices of male and female samples

## Value

`check_sex` returns the normalized average total intensities of probes
targeting the X and Y chromosomes. These are forwarded to `predict_sex`
which returns a factor with levels "f","m" (and `NA` if the sex cannot
be determined conclusively).
