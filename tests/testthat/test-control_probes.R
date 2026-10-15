testthat::test_that("control probe metrics return correct results", {

  k450 = testthat::test_path("../../data/5640269011_R01C01_Grn.idat")
  stopifnot(file.exists(k450))

  k450 = sub("_Grn\\.idat$", "", k450)
  k450 = read_idats(k450, quiet = TRUE)
  k450 = control_metrics(k450)

  epic = testthat::test_path("../../data/200144450018_R04C01_Grn.idat")
  stopifnot(file.exists(epic))

  epic = sub("_Grn\\.idat$", "", epic)
  epic = read_idats(epic, quiet = TRUE)
  epic = control_metrics(epic)

  metrics = cbind(k450, epic)

  testthat::expect_equal_to_reference(metrics, file = "metrics.rds") 
})