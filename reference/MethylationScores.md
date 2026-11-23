# Various methylation scores

Calculate various epigenetic scores (e.g. maternal smoking, gestational
age, Horvath clock) based on beta-values of DNA methylation.

## Usage

``` r
methylation_score(beta, model)
```

## Arguments

- beta:

  Matrix of beta-values with probe IDs as row names and samples as
  columns.

- model:

  Name of the model. Currently supported models are:

  - `"cord_blood:maternal_smoking"`: Score for in utero exposure to
    maternal smoking based on cord blood \[1\].

  - `"cord_blood:gestational_age"`: Gestational age (in weeks) based on
    cord blood \[2\].

  - `"placenta:gestational_age"`: Gestational age (in weeks) based on
    placenta \[3\].

  - `"horvath_clock"`: Multi-tissue epigenetic age (in years) \[4\].

## Value

A numeric vector of calculated methylation scores, one for each sample
(column) of `beta`. Returns `NULL` if `model` is not one of the
supported models.

## References

\[1\] Reese SE, et al. DNA methylation score as a biomarker in newborns
for sustained maternal smoking during pregnancy. Environmental Health
Perspectives. 2017;125(4):760-6. DOI:10.1289/EHP333

\[2\] Knight AK, et al. An epigenetic clock for gestational age at birth
based on blood methylation data. Genome Biology. 2016;17:206.
DOI:10.1186/s13059-016-1068-z

\[3\] Lee Y, et al. Placental epigenetic clocks: estimating gestational
age using placental DNA methylation levels. Aging. 2019;11(12):4238-53.
DOI:10.18632/aging.102049

\[4\] Horvath S. DNA methylation age of human tissues and cell types.
Genome Biology. 2013;14:R115. DOI:10.1186/gb-2013-14-10-r115
