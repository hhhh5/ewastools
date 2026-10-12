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
  
  # Check for perfect 1-to-1 mapping between metadata subjects and methylation donors
  genotypes = call_genotypes(beta_snps, learn = FALSE)
  pheno_df$donor_id = enumerate_sample_donors(genotypes)
  
  pheno_df |>
    dplyr::count(subject_id, donor_id) |>
    testthat::expect_shape(nrow = 4)
  
  # --- BENCHMARK 2: Reproducibility of Cell Type Deconvolution ---
  meth = detectionP(meth)
  meth = mask(meth, 0.01)
  meth = correct_dye_bias(meth)
  beta_full = beta_values(meth)
  
  # Estimate leukocyte composition (using HRS reference for blood)
  LC = estimateLC(beta_full, ref="HRS")
  pheno_df = cbind(pheno_df, as.data.frame(LC))
  
  cell_types = colnames(LC)
  
  pheno_df |>
    tidyr::pivot_longer(cols = dplyr::all_of(cell_types),
      names_to = "Cell Type", values_to = "Proportion") |>
  ggplot2::ggplot(ggplot2::aes(x = subject_id, y = Proportion, color = subject_id)) +
  ggplot2::geom_boxplot() +
  ggplot2::scale_y_continuous(trans = "log10") +
  ggplot2::facet_wrap(ggplot2::vars(`Cell Type`), labeller = ggplot2::label_both)
})
