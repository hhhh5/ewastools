# Drop samples

Drop samples from the object returned by
[`read_idats()`](https://hhhh5.github.io/ewastools/reference/read_idats.md).
Used for removing samples that failed quality control before computing
beta-values.

## Usage

``` r
drop_samples(raw, j = NULL)
```

## Arguments

- raw:

  Output of calling
  [`read_idats`](https://hhhh5.github.io/ewastools/reference/read_idats.md)

- j:

  Indices of the samples to drop

## Value

A modified `raw` object
