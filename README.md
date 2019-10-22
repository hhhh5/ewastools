# ewas_tools: a quality control toolset for the Illumina Infinium DNA methylation platforms

The following functionality is offered for the 450K and EPIC chips:

- Screen for problematic samples, e.g. failed assays, mislabeled or contaminated samples (`control_metrics, check_sex, snp_outliers`).
- Compute detection p-values and mask respective data points (`detectionP, mask`).
- Preprocess the data (`correct_dye_bias, normalize, dont_normalize`).
- Estimate leukocyte composition in case of blood samples (`estimateLC`)

**Read the vignette on the recommended QC workflow at <https://hhhh5.github.io/ewastools/articles/exemplary_ewas.html>.**

Package documentation is available at <https://hhhh5.github.io/ewastools/>.


An open access paper describing the quality checks implemented in this package in detail is available at <http://doi.org/10.1186/s13148-018-0504-1>. A vignette demonstrating the application of quality control checks on a public dataset is provided as well.

A second paper (open access) describing the improved detection p-values implemented in `ewastools` can be found at <https://doi.org/10.1186/s13148-019-0615-3>.

**For installation please run `devtools::install_github("hhhh5/ewastools")`**
