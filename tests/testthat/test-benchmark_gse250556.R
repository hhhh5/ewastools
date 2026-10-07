testthat::test_that("GSE250556 benchmarks: sample tracking and cell type deconvolution", {
  # Locate the external data directory
  data_dir = testthat::test_path("../../../ewastools-data", "GSE250556_RAW")
  
  testthat::skip_if_not(dir.exists(data_dir), "GSE250556_RAW data is not available.")
  
  # 1. Scraping metadata from GEO...
  pheno_url = "https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE250556&targ=gsm&form=text&view=brief"
  pheno = readLines(pheno_url)
  pheno = split(pheno, cumsum(stringr::str_detect(pheno, "^\\^SAMPLE = GSM")))
  names(pheno) = purrr::map_chr(pheno, 1) |> stringr::str_extract("GSM\\d+")
  
  pheno_df = pheno |>
    purrr::map(tibble::tibble) |>
    dplyr::bind_rows(.id = "gsm") |>
    tidyr::separate_wider_delim(2, delim = " = ", names = c("var", "val"), too_many="merge", too_few="align_start") |>
    dplyr::filter(var == "!Sample_title") |>
    dplyr::mutate(subject_id = stringr::str_extract(val, "subject[A-Z]"))
    
  # 2. Finding local IDAT files...
  idats = list.files(data_dir, pattern = "_Grn\\.idat(\\.gz)?$", full.names = TRUE)
  files_df = tibble::tibble(
    prefix = stringr::str_remove(idats, "_Grn\\.idat(\\.gz)?$"),
    gsm = stringr::str_extract(basename(prefix), "^GSM\\d+"))
  
  # Merge metadata with local files. unmatched = "error" ensures no files or metadata are missing.
  pheno_df = dplyr::inner_join(pheno_df, files_df, by = "gsm", unmatched = "error")
  
  prefixes = pheno_df$prefix
  
  # 3. Reading IDAT files...
  meth = read_idats(prefixes, quiet = TRUE, on_disk = TRUE)
  
  # --- BENCHMARK 1: Sample Tracking and Mix-up Detection ---
  # meth$manifest is a data.table, but we can subset it using base R
  snps = meth$manifest$index[meth$manifest$probe_type == "rs"]
  
  # For genotype calling we need beta values for SNPs
  M_snps = meth$M[snps, , drop = FALSE]
  U_snps = meth$U[snps, , drop = FALSE]
  beta_snps = M_snps / (M_snps + U_snps + 100)
  
  genotypes = call_genotypes(beta_snps, learn = FALSE)
  
  # Group samples that have matching genotypes
  pheno_df$donor_id = enumerate_sample_donors(genotypes)
  
  # Check for perfect 1-to-1 mapping between metadata subjects and methylation donors
  
  # 1. Each metadata subject should have exactly 1 donor ID
  donors_per_subject = pheno_df |>
    dplyr::group_by(subject_id) |>
    dplyr::summarise(n = dplyr::n_distinct(donor_id))
    
  testthat::expect_true(all(donors_per_subject$n == 1))
  
  # 2. Each methylation donor ID should correspond to exactly 1 metadata subject
  subjects_per_donor = pheno_df |>
    dplyr::group_by(donor_id) |>
    dplyr::summarise(n = dplyr::n_distinct(subject_id))
    
  testthat::expect_true(all(subjects_per_donor$n == 1))
  
  # --- BENCHMARK 2: Reproducibility of Cell Type Deconvolution ---
  meth = detectionP(meth)
  meth = mask(meth, 0.01)
  meth = correct_dye_bias(meth)
  beta_full = beta_values(meth)
  
  # Estimate leukocyte composition (using HRS reference for blood)
  LC = estimateLC(beta_full, ref="HRS")
  pheno_df = cbind(pheno_df, as.data.frame(LC))
  
  # Calculate variance across technical replicates for each subject
  cell_types = colnames(LC)
  variance_summary = pheno_df |>
    dplyr::group_by(subject_id) |>
    dplyr::summarise(dplyr::across(dplyr::all_of(cell_types), ~ var(.x, na.rm = TRUE)))
    
  mean_vars = colMeans(variance_summary[, cell_types, drop = FALSE], na.rm = TRUE)
  
  # We expect the mean intra-subject variance for every cell type to be extremely low
  testthat::expect_true(all(mean_vars < 0.01))
  
  # --- VISUAL REGRESSION TESTING ---
  # Create a long format dataframe for plotting
  plot_df = pheno_df |>
    tidyr::pivot_longer(cols = dplyr::all_of(cell_types), names_to = "Cell_Type", values_to = "Proportion")
    
  p = ggplot2::ggplot(plot_df, ggplot2::aes(x = subject_id, y = Proportion, color = subject_id)) +
    ggplot2::geom_boxplot(outlier.shape = NA) +
    ggplot2::geom_jitter(width = 0.2, alpha = 0.7) +
    ggplot2::facet_wrap(~ Cell_Type, scales = "free_y") +
    ggplot2::theme_minimal() +
    ggplot2::labs(
      title = "Cell Type Deconvolution across Technical Replicates",
      x = "Subject ID",
      y = "Estimated Proportion"
    ) +
    ggplot2::theme(legend.position = "none", axis.text.x = ggplot2::element_text(angle = 45, hjust = 1))
    
  # visual regression test: this generates/checks tests/testthat/_snaps/leukocyte_variance_plot.svg
  vdiffr::expect_doppelganger("leukocyte_variance_plot", p)
})
