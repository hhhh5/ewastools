# Read .idat files

Read .idat files

## Usage

``` r
read_idats(idat_files, quiet = FALSE, on_disk = FALSE)
```

## Arguments

- idat_files:

  Character vector of relative or absolute filepaths, but without the
  suffixes \_Grn.idat' and '\_Red.idat'. IDATs for red and green channel
  must have the same prefix and be stored in the same folder. E.g., a
  sample with the idats 200607110235_R01C01_Red.idat and
  200607110235_R01C01_Grn.idat would be passed to `read_idats` as
  "200607110235_R01C01".

- quiet:

  If `TRUE`, suppresses the progress bar (useful for RMarkdown scripts).

- on_disk:

  If `TRUE`, data will be stored on disk using the ff package. This will
  be slower but enables processing datasets larger than memory.

## Value

A list containing

- manifest:

  A data.table describing the probes

- M:

  intensities of targeting methylated sequences

- U:

  intensities of targeting unmethylated sequences

- N:

  number of beads from which average intensities in M were derived

- V:

  number of beads from which average intensities in U were derived

- ctrlG:

  Intensities of the control probes in the green color channel

- ctrlR:

  Intensities of the control probes in the red color channel

- meta:

  A data.table containing unique sample_ids and metadata

## Examples

``` r
if (FALSE) { # \dontrun{
read_idats('9976861004_R01C01')
} # }
```
