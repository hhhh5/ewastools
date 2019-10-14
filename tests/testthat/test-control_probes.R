library(stringi)

test_that("control probe metrics return correct results", {

  k450 = system.file(paste0("data/5640269011_R01C01_Grn.idat"),package="ewastools")
  k450 = stri_sub(k450,1,-10)
  k450 = read_idats(k450)
  k450 = control_metrics(k450)

  epic = system.file(paste0("data/200144450018_R04C01_Grn.idat"),package="ewastools")
  epic = stri_sub(epic,1,-10)
  epic = read_idats(epic)
  epic = control_metrics(epic)

  metrics = cbind(k450,epic)  

  expect_equal_to_reference(metrics,file="metrics.rds") 

  #                              epic       k450      
  # Restoration                  0.06311688 0.04903418
  # Staining Green               57.12603   295.8571  
  # Staining Red                 45.75246   362.5732  
  # Extension Green              49.15882   63.74247  
  # Extension Red                19.75773   18.2598   
  # Hybridization High/Medium    1.547735   1.791453  
  # Hybridization Medium/Low     2.095093   1.552534  
  # Target Removal 1             10.57692   16.825    
  # Target Removal 2             9.101655   11.06908  
  # Bisulfite Conversion I Green 13.02981   28.45172  
  # Bisulfite Conversion I Red   9.493972   9.404153  
  # Bisulfite Conversion II      5.32149    7.090308  
  # Specificity I Green          16.97826   16.55665  
  # Specificity I Red            6.834112   3.780076  
  # Specificity II               23.56429   24.39921  
  # Non-polymorphic Green        16.21673   9.844952  
  # Non-polymorphic Red          10.25316   23.657  

	

})