library(stringi)

test_that("control probe metrics return correct results", {

  k450 = testthat::test_path("../../data/5640269011_R01C01_Grn.idat")
  k450 = stri_sub(k450,1,-10)
  k450 = read_idats(k450)
  k450 = control_metrics(k450)

  epic = testthat::test_path("../../data/200144450018_R04C01_Grn.idat")
  epic = stri_sub(epic,1,-10)
  epic = read_idats(epic)
  epic = control_metrics(epic)

  metrics = cbind(k450,epic)

  expect_equal_to_reference(metrics, file = "metrics.rds") 
})