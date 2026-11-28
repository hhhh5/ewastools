# ewastools 1.8.2

## New
* Add support for Infinium Methylation Screening Array (MSA) BeadChip
* Add a benchmark vignette using a dataset with technical replicates (GSE250556). Vignette covers sample donor identification and cell type deconvolution.

## Fixes
* Replace `ff::is.ff()` calls with `inherits(..., "ff")` so in-memory workflows no longer fail if the optional package `ff` is not installed.
* Fix `call_genotypes(..., learn = TRUE)`, which previously failed in `eBeta()` due to checking an undefined `weights` variable
* Anonymized IDATs with missing/empty `RunInfo` previously errored on `nrow(red$RunInfo) > 1`
* Zero-bead out-of-band intensities (`oobG`/`oobR`) are assigned `NA` instead being left as `0`
* Extracting `sample_id` via `strsplit(..., "/")` failed on Windows backslash paths. It now uses `basename()`.
* Fix order of calling pre-processing functions in `exemplary_ewas`: `check_sex()` should be called before `mask()`.
* Use `log10(0.01)` threshold instead of `0.01` in `detectionP.neg()` 
* Fix `methylation_score()` in case of missing probes: previously, `NA`s caused errors for on-disk (`ff`) matrices, and samples with no observed model probes silently returned the bare model intercept instead of `NA`.
* Fix single-sample / single-sex dimension dropping (`drop = FALSE`) in `eval_detP_cutoffs()`, `detectionP.minfi()`, and `correct_dye_bias2()`
* Fix `drop_samples()` when passing an empty index (`integer(0)`, e.g. when no samples fail QC). This previously dropped all samples instead of keeping them

## Improvements
* Reduce memory usage by loading platform manifests lazily
* Pre-knit vignettes so that package/website builds do not re-download large GEO datasets
* Improved `read_idats()` performance by caching bead address matching and manifest channel subsets across samples with identical `MidBlock` ordering
* Standardize argument validation across exported functions using `checkmate`
* Expand unit tests covering QC, preprocessing, genotyping, deconvolution, and methylation score pipeline.

# ewastools 1.8.1
* Added experimental support for Mouse Methylation BeadChip (feedback welcome)
* Added `on_disk` option to `read_idats()`. Uses `ff` package to store methylation data on disk instead
  memory. While slower, should enable analyses of large datasets on small machines.
* Update vignette

# ewastools 1.8.0
* Removed deconvolution models based on combinations of two datasets
* Improved support for EPICv2 chip, adding a function to handle new probe IDs
* Refactored parsing of manifests files

# ewastools 1.7.2
* Added support for EPICv2 chip

# ewastools 1.7.1
* Added HRS cell composition model

# ewastools 1.7
* Saliva cell composition prediction
* Test `detectionP.minfi` for compatibility
* Drop `normalize`

# ewastools 1.6
* `read_idats`: extract Sentrix barcode and position from .idat files, stored in `meta` data.table (92eb8df289a473fa0c06471be197db1953e39efd0)
* `call_genotypes`: use pre-specified parameters for SNP mixture model by default, i.e., `learn=FALSE` (507abceeb8db572d72089e423678f3cd5a7f6769)
* Update recommedation of pre-processing workflow: call `mask()` after `check_sex()` (eb54d1b0d74d2b0231db6e03a8c60710abd42e09)
* Add some unit tests
* Fix: some bisulfite conversion control probes were ignored due to inconsistent naming (24233ffcbea4337b9a7d4d3060b03b484d5ba01c)