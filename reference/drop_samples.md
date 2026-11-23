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

  Output of
  [`read_idats`](https://hhhh5.github.io/ewastools/reference/read_idats.md).

- j:

  Indices (column numbers) of the samples to drop.

## Value

A modified `raw` object without the dropped samples.
