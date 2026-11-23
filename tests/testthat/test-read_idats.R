data_dir <- testthat::test_path("../../../ewastools-data")

test_that("read_idats imports 27K platform", {
  skip_if_not(dir.exists(data_dir), "Local test data directory 'ewastools-data' not found.")
  idat_prefix <- file.path(data_dir, "GSE17769_HM27/GSM443811_HCC38_4308918033_G")
  
  if (!file.exists(paste0(idat_prefix, "_Grn.idat.gz")))
    skip("Test data for 27K not found")

  capture.output(res <- read_idats(idat_prefix, quiet = TRUE))
  
  expect_equal(res$platform, "27K")
  expect_true(is.matrix(res$M))
  expect_equal(nrow(res$M), 27594)
  expect_equal(ncol(res$M), 1)
  expect_s3_class(res$meta, "data.frame")
  expect_equal(nrow(res$meta), 1)
})

test_that("read_idats imports 450K platform", {
  skip_if_not(dir.exists(data_dir), "Local test data directory 'ewastools-data' not found.")
  idat_prefix <- file.path(data_dir, "reinius/5684819001_R01C01")
  
  if (!file.exists(paste0(idat_prefix, "_Grn.idat")))
    skip("Test data for 450K not found")

  capture.output(res <- read_idats(idat_prefix, quiet = TRUE))
  
  expect_equal(res$platform, "450K")
  expect_true(is.matrix(res$M))
  expect_equal(nrow(res$M), 485577)
  expect_equal(ncol(res$M), 1)
  expect_s3_class(res$meta, "data.frame")
  expect_equal(nrow(res$meta), 1)
})

test_that("read_idats imports EPICv1 platform", {
  skip_if_not(dir.exists(data_dir), "Local test data directory 'ewastools-data' not found.")
  idat_prefix <- file.path(data_dir, "GSE110554/GSM2998021")
  
  if (!file.exists(paste0(idat_prefix, "_Grn.idat.gz")))
    skip("Test data for EPICv1 not found")
  
  capture.output(res <- read_idats(idat_prefix, quiet = TRUE))
  
  expect_equal(res$platform, "EPICv1")
  expect_true(is.matrix(res$M))
  expect_equal(nrow(res$M), 865918)
  expect_equal(ncol(res$M), 1)
  expect_s3_class(res$meta, "data.frame")
  expect_equal(nrow(res$meta), 1)
})

test_that("read_idats imports EPICv2 platform", {
  skip_if_not(dir.exists(data_dir), "Local test data directory 'ewastools-data' not found.")
  idat_prefix <- file.path(data_dir, "DemoDataEPICv2_NextSeq550/221004_NB551015_0621_206891110004/206891110004/206891110004_R01C01")
  
  if (!file.exists(paste0(idat_prefix, "_Grn.idat")))
    skip("Test data for EPICv2 not found")
  
  capture.output(res <- read_idats(idat_prefix, quiet = TRUE))
  
  expect_equal(res$platform, "EPICv2")
  expect_true(is.matrix(res$M))
  expect_equal(nrow(res$M), 937055)
  expect_equal(ncol(res$M), 1)
  expect_s3_class(res$meta, "data.frame")
  expect_equal(nrow(res$meta), 1)
})

test_that("read_idats imports MSA platform", {
  skip_if_not(dir.exists(data_dir), "Local test data directory 'ewastools-data' not found.")
  idat_prefix <- file.path(data_dir, "MSA-demo/207760740041/207760740041_R03C01")
  
  if (!file.exists(paste0(idat_prefix, "_Grn.idat")))
    skip("Test data for MSA not found")
  
  capture.output(res <- read_idats(idat_prefix, quiet = TRUE))
  
  expect_equal(res$platform, "MSA")
  expect_true(is.matrix(res$M))
  expect_equal(nrow(res$M), 281798)
  expect_equal(ncol(res$M), 1)
  expect_s3_class(res$meta, "data.frame")
  expect_equal(nrow(res$meta), 1)
})

test_that("read_idats imports MOUSE platform", {
  skip_if_not(dir.exists(data_dir), "Local test data directory 'ewastools-data' not found.")
  idat_prefix <- file.path(data_dir, "MOUSE/204617710004_R01C01")
  
  if (!file.exists(paste0(idat_prefix, "_Grn.idat")))
    skip("Test data for MOUSE not found")
  
  capture.output(res <- read_idats(idat_prefix, quiet = TRUE))
  
  expect_equal(res$platform, "MOUSE")
  expect_true(is.matrix(res$M))
  expect_equal(nrow(res$M), 287050)
  expect_equal(ncol(res$M), 1)
  expect_s3_class(res$meta, "data.frame")
  expect_equal(nrow(res$meta), 1)
})

test_that("read_idats and drop_samples work on bundled 450K and EPICv1 files", {
  k450 = testthat::test_path("../../data/5640269011_R01C01")
  epic = testthat::test_path("../../data/200144450018_R04C01")
  stopifnot(
    file.exists(paste0(k450, "_Grn.idat")),
    file.exists(paste0(epic, "_Grn.idat"))
  )

  res_epic = read_idats(epic, quiet = TRUE)
  expect_equal(res_epic$platform, "EPICv1")
  expect_equal(dim(res_epic$M), c(865918, 1))
  expect_equal(res_epic$meta$sample_id, "200144450018_R04C01")

  res2 = read_idats(c(k450, k450), quiet = TRUE)
  expect_equal(res2$platform, "450K")
  expect_equal(dim(res2$M), c(485577, 2))
  expect_equal(res2$M[, 1], res2$M[, 2])

  # Empty `j` and `NULL` keep all samples
  expect_equal(ncol(drop_samples(res2, NULL)$M), 2)
  expect_equal(ncol(drop_samples(res2, integer(0))$M), 2)
  expect_equal(ncol(drop_samples(res2, c(FALSE, FALSE))$M), 2)

  # Integer and logical `j` drop the specified sample
  d_int = drop_samples(res2, 1L)
  d_log = drop_samples(res2, c(TRUE, FALSE))
  expect_equal(ncol(d_int$M), 1)
  expect_equal(nrow(d_int$meta), 1)
  expect_equal(d_int$M, d_log$M)
})

test_that("read_idats and drop_samples work with on_disk = TRUE", {
  skip_if_not_installed("ff")
  k450 = testthat::test_path("../../data/5640269011_R01C01")
  stopifnot(file.exists(paste0(k450, "_Grn.idat")))

  res_ff = read_idats(c(k450, k450), quiet = TRUE, on_disk = TRUE)
  expect_s3_class(res_ff$M, "ff")
  expect_true(is.matrix(res_ff$ctrlG))

  dropped_ff = drop_samples(res_ff, 2L)
  expect_s3_class(dropped_ff$M, "ff")
  expect_equal(ncol(dropped_ff$M), 1)
  expect_true(is.matrix(dropped_ff$ctrlG))
  expect_equal(ncol(dropped_ff$ctrlG), 1)
})
