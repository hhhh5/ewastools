test_that("sample_failure handles single-sample metrics correctly", {
  # Mock metrics list for 1 sample
  metrics = list(
    m1 = structure(c(0.5), threshold = 1.0),
    m2 = structure(c(2.0), threshold = 1.0)
  )

  # Sample should fail because m1 (0.5) < threshold (1.0)
  sf = sample_failure(metrics)
  expect_type(sf, "logical")
  expect_length(sf, 1)
  expect_true(sf)

  # Passing sample
  metrics_pass = list(
    m1 = structure(c(1.5), threshold = 1.0),
    m2 = structure(c(2.0), threshold = 1.0)
  )
  sf_pass = sample_failure(metrics_pass)
  expect_length(sf_pass, 1)
  expect_false(sf_pass)
})

test_that("methylation_score handles single-sample matrix and warns on missing probes", {
  # Synthetic 1-sample beta matrix
  probes = c("cg02256631", "cg02482603", "cg04103532")
  beta = matrix(c(0.1, 0.5, 0.8), nrow = 3, ncol = 1,
                dimnames = list(probes, "sample1"))

  # Model expects 26 probes; with only 3 present, it should issue a warning but calculate anyway
  expect_warning(
    score <- methylation_score(beta, "cord_blood:maternal_smoking"),
    "model probes were not found")
  expect_length(score, 1)
  expect_false(is.na(score))
})

test_that("read_manifest lazy loads and caches manifest and controls", {
  m27k = ewastools:::read_manifest("27K")
  expect_type(m27k, "list")
  expect_named(m27k, c("manifest", "controls"))
  expect_s3_class(m27k$manifest, "data.frame")
  expect_s3_class(m27k$controls, "data.frame")
})

test_that("read_idats does not mutate cached manifests or controls in place", {
  k450 = testthat::test_path("../../data/5640269011_R01C01_Grn.idat")
  if (file.exists(k450)) {
    idat_base = sub("_Grn\\.idat$", "", k450)
    raw1 = read_idats(idat_base, quiet = TRUE)
    # Cached 450K manifest should still have chr37 and not index or OOBi
    m450 = ewastools:::read_manifest("450K")
    expect_true("chr37" %in% names(m450$manifest))
    expect_false("OOBi" %in% names(m450$manifest))
    expect_false("index" %in% names(m450$manifest))
    expect_false("index" %in% names(m450$controls))

    # A second read_idats call on the same platform should succeed cleanly
    expect_no_error(raw2 <- read_idats(idat_base, quiet = TRUE))
  }
})
