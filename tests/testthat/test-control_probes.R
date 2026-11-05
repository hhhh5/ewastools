testthat::test_that("control probe metrics return correct results", {
  k450 = testthat::test_path("../../data/5640269011_R01C01")
  epic = testthat::test_path("../../data/200144450018_R04C01")
  ref  = testthat::test_path("metrics.rds")

  stopifnot(
    k450 = file.exists(paste0(k450, "_Grn.idat")),
    epic = file.exists(paste0(epic, "_Grn.idat")),
    ref  = file.exists(ref))

  k450 = read_idats(k450, quiet = TRUE)
  k450_metrics = control_metrics(k450)

  epic = read_idats(epic, quiet = TRUE)
  epic_metrics = control_metrics(epic)

  metrics = cbind(k450 = k450_metrics, epic = epic_metrics)

  ref = readRDS(testthat::test_path("metrics.rds"))
  testthat::expect_equal(metrics, ref)
})

testthat::test_that("sample_failure flags samples below threshold and ignores NAs", {
  k450 = testthat::test_path("../../data/5640269011_R01C01")
  stopifnot(file.exists(paste0(k450, "_Grn.idat")))

  raw2 = read_idats(c(k450, k450), quiet = TRUE)
  ctrls = control_metrics(raw2)

  failed = sample_failure(ctrls)
  testthat::expect_type(failed, "logical")
  testthat::expect_length(failed, 2)
  testthat::expect_false(any(is.na(failed)))

  # Force sample 1 below threshold on one metric and sample 2 to NA
  ctrls_mod = ctrls
  thresh = attr(ctrls_mod$`Staining Green`, "threshold")
  ctrls_mod$`Staining Green`[1] = thresh - 0.1
  ctrls_mod$`Staining Green`[2] = NA_real_

  failed_mod = sample_failure(ctrls_mod)
  testthat::expect_true(failed_mod[1])
  testthat::expect_equal(failed_mod[2], failed[2])
})