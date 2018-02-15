# ewas_tools: a quality control toolset for the Illumina Infinium DNA methylation platforms

The following functionality is offered for the 450K and EPIC chips:

- Screen for problematic samples, e.g. failed assays, mislabeled or contaminated samples (`control_metrics, check_sex, snp_outliers`).
- Compute detection p-values and mask respective data points (`detectionP, mask`).
- Preprocess the data (`correct_dye_bias, normalize, dont_normalize`).
- Estimate leukocyte composition in case of blood samples (`estimateLC`)

A vignette demonstrating the application of quality control checks on a public dataset is provided.

**For installation please use the precompiled binaries or run `devtools::install_github("hhhh5/ewastools")`**
