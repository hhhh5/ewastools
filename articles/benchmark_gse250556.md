# GSE250556 Benchmarks: Sample Tracking and Cell Type Deconvolution

This vignette demonstrates two key quality control benchmarks using the
public dataset **GSE250556**, which contains repeated technical
replicates across multiple subjects:

1.  **Sample Tracking and Mix-up Detection:** Inferring donor identity
    from SNP probes to verify concordance between metadata annotations
    and genetic fingerprints
    ([`call_genotypes()`](https://hhhh5.github.io/ewastools/reference/call_genotypes.md),
    [`enumerate_sample_donors()`](https://hhhh5.github.io/ewastools/reference/check_snp_agreement.md)).
2.  **Cell Type Deconvolution Reproducibility:** Assessing the
    consistency of leukocyte composition estimation
    ([`estimateLC()`](https://hhhh5.github.io/ewastools/reference/cell_composition.md))
    across technical replicates of the same individual.

``` r

library(ewastools)
library(dplyr)
library(tidyr)
library(ggplot2)
```

## Data Preparation

We first locate the raw IDAT files and retrieve sample metadata from
GEO.

``` r

# Locate directory containing GSE250556 raw IDAT files
data_dir = if (dir.exists("../../ewastools-data/GSE250556_RAW")) {
  "../../ewastools-data/GSE250556_RAW"
} else if (dir.exists("../ewastools-data/GSE250556_RAW")) {
  "../ewastools-data/GSE250556_RAW"
}
```

``` r

# 1. Fetch metadata from GEO
pheno_url = "https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE250556&targ=gsm&form=text&view=brief"
pheno = readLines(pheno_url)
pheno = split(pheno, cumsum(stringr::str_detect(pheno, "^\\^SAMPLE = GSM")))
names(pheno) = purrr::map_chr(pheno, 1) |> stringr::str_extract("GSM\\d+")

pheno_df = pheno |>
  purrr::map(tibble::tibble) |>
  dplyr::bind_rows(.id = "gsm") |>
  tidyr::separate_wider_delim(2, delim = " = ", names = c("var", "val"), too_many = "merge", too_few = "align_start") |>
  dplyr::filter(var == "!Sample_title") |>
  dplyr::mutate(subject_id = stringr::str_extract(val, "subject[A-Z]"))

# 2. Find local IDAT files
idats = list.files(data_dir, pattern = "_Grn\\.idat(\\.gz)?$", full.names = TRUE)
files_df = tibble::tibble(
  prefix = stringr::str_remove(idats, "_Grn\\.idat(\\.gz)?$"),
  gsm = stringr::str_extract(basename(prefix), "^GSM\\d+"))

# Merge metadata with local files
pheno_df = dplyr::inner_join(pheno_df, files_df, by = "gsm", unmatched = "error")
```

``` r

# 3. Read IDAT files
meth = read_idats(pheno_df$prefix, quiet = TRUE, on_disk = TRUE)
```

## Benchmark 1: Sample Tracking and Mix-up Detection

Illumina DNA methylation microarrays contain a panel of SNP (single
nucleotide polymorphism) probes (prefixed with `rs`). These probes can
be used to call genotypes and verify that samples collected from the
same individual indeed share the same genetic profile.

``` r

# Identify SNP probes in manifest
snps = meth$manifest$index[meth$manifest$probe_type == "rs"]

# Calculate beta values for SNPs
M_snps = meth$M[snps, , drop = FALSE]
U_snps = meth$U[snps, , drop = FALSE]
beta_snps = M_snps / (M_snps + U_snps + 100)

# Call genotypes and cluster samples by shared donor profile
genotypes = call_genotypes(beta_snps, learn = FALSE)
pheno_df$donor_id = enumerate_sample_donors(genotypes)

# Verify mapping between metadata subjects and methylation-inferred donors
pheno_df |>
  dplyr::count(subject_id, donor_id) |>
  print()
```

As shown in the table above, all replicates for each subject group
together into a unique inferred donor ID, confirming that there were no
sample mislabelings or mix-ups.

## Benchmark 2: Reproducibility of Cell Type Deconvolution

Whole blood consists of heterogeneous cell populations. Cell type
deconvolution estimates the proportion of major leukocyte cell types
from DNA methylation profiles. Across technical replicates of the same
biological sample, cell type estimates should show high consistency and
minimal variance.

### Pre-processing

Before cell type estimation, we apply detection p-value masking and
dye-bias correction:

``` r

meth = detectionP(meth)
meth = mask(meth, 0.01)
meth = correct_dye_bias(meth)
beta_full = beta_values(meth)
```

### Estimating Leukocyte Composition

We estimate cell proportions using the reference dataset for blood
(`ref = "HRS"`):

``` r

LC = estimateLC(beta_full, ref = "HRS")
pheno_df = cbind(pheno_df, as.data.frame(LC))
cell_types = colnames(LC)

head(LC)
```

### Visualizing Deconvolution Reproducibility

We reshape the data to plot the estimated proportions of each cell type
across technical replicates, grouped by subject:

``` r

plot_df = pheno_df |>
  tidyr::pivot_longer(
    cols = dplyr::all_of(cell_types),
    names_to = "Cell Type",
    values_to = "Proportion"
  )

ggplot(plot_df, aes(x = subject_id, y = Proportion, color = subject_id)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.5) +
  geom_jitter(width = 0.15, size = 1.8, alpha = 0.8) +
  scale_y_continuous(trans = "log10") +
  facet_wrap(vars(`Cell Type`), scales = "free_y") +
  theme_minimal() +
  labs(
    title = "Cell Type Deconvolution across Technical Replicates",
    subtitle = "GSE250556 leukocyte composition estimates grouped by subject",
    x = "Subject ID",
    y = "Estimated Proportion (log10 scale)",
    color = "Subject"
  ) +
  theme(
    legend.position = "none",
    strip.text = element_text(face = "bold", size = 10),
    axis.text.x = element_text(angle = 45, hjust = 1)
  )
```

The plot illustrates high consistency across technical replicates within
each subject, demonstrating robust reproducibility for
[`estimateLC()`](https://hhhh5.github.io/ewastools/reference/cell_composition.md).
