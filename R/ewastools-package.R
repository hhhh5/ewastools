#' @import data.table
#' @importFrom graphics axis boxplot curve hist legend plot points
#' @importFrom stats coef dbeta density mad median na.omit pnorm quantile
#' @importFrom utils head read.table setTxtProgressBar txtProgressBar
NULL

utils::globalVariables(c(
   "M", "Mi", "OOBi", "U", "Ui",
   "address", "addressM", "addressU", "agreement",
   "channel", "chr", "controls", "ctrlG", "ctrlR",
   "donor1", "donor2", "group", "i", "index",
   "manifest", "name", "par", "snps", "x"
))
