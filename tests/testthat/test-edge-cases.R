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
