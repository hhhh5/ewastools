# Read .idat files

Import and annotate the data from .idat files.

## Usage

``` r
read_idats(idat_files, quiet = FALSE, on_disk = FALSE)
```

## Arguments

- idat_files:

  Character vector of relative or absolute filepaths without
  "\_Grn.idat" and "\_Red.idat" suffixes. IDATs for red and green
  channel must have the same prefix and be stored in the same folder.
  Gzipped files (".idat.gz") are supported as well.

  E.g., a sample with the files "file/path/sample1_Red.idat" and
  "file/path/sample1_Grn.idat" would be passed as "file/path/sample1".

  All .idat files are assumed to be from the same platform.

- quiet:

  If `TRUE`, suppresses the progress bar (useful for RMarkdown scripts).

- on_disk:

  If `TRUE`, data will be stored on disk using the `ff` package. This
  will be slower but enables processing datasets larger than memory.

## Value

A list containing

- platform:

  Name of the BeadChip (27K/450K/EPICv1/EPICv2/MSA/MOUSE)

- manifest:

  A data.table describing the probes

- M:

  Intensities of probes targeting methylated sequences

- U:

  Intensities of probes targeting unmethylated sequences

- S:

  Standard deviations of the bead intensities underlying `M`

- T:

  Standard deviations of the bead intensities underlying `U`

- N:

  Number of beads from which average intensities in `M` were derived

- V:

  Number of beads from which average intensities in `U` were derived

- controls:

  A data.table describing the control probes

- ctrlG:

  Intensities of the control probes in the green color channel

- ctrlR:

  Intensities of the control probes in the red color channel

- ctrlN:

  Number of beads underlying the control probe intensities

- oobG:

  List of matrices `M` and `U` with out-of-band intensities in the green
  color channel (of Infinium type I probes measured in the red color
  channel)

- oobR:

  List of matrices `M` and `U` with out-of-band intensities in the red
  color channel (of Infinium type I probes measured in the green color
  channel)

- meta:

  A tibble with sample metadata: `sample_id` (file name without path),
  scan `date` and `time`, Sentrix `barcode` and `position`

The intensity matrices have one row per probe (in the order of
`manifest`) and one column per sample. If `on_disk = TRUE`, they are
`ff` objects (except for control probes).

## Examples

``` r
if (FALSE) { # \dontrun{
read_idats('9976861004_R01C01')
} # }
```
