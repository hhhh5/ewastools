#' Infer the sex of sample donor
#' 
#' Return the normalized average total intensities of probes targeting the X and Y chromosomes. 
#' 
#' @rdname check_sex
#' @param raw Output of \code{\link{read_idats}}
#' @return \code{check_sex} returns a list containing two numeric vectors:
#' \item{X}{Normalized average total intensities of probes targeting the X chromosome.}
#' \item{Y}{Normalized average total intensities of probes targeting the Y chromosome.}
#' @export
check_sex = function(raw) {
   
   if(!all(c('M', 'U', 'manifest') %in% names(raw))) stop('Invalid argument')

   with(raw, {

      # Select allosomal (X and Y chromosome) and autosomal probes from manifest
      i_chrX = manifest[ forcats::fct_match(chr,   "chrX"),         index]
      i_chrY = manifest[ forcats::fct_match(chr,   "chrY"),         index]
      i_auto = manifest[!forcats::fct_match(chr, c("chrX","chrY")), index]

      J = ncol(M)
      chrX = numeric(J)
      chrY = numeric(J)
      auto = numeric(J)

      # Calculate mean total intensity (M + U) for X, Y, and autosomes per sample
      for (j in 1:J) {
         chrX[j] = mean(M[i_chrX, j] + U[i_chrX, j], na.rm = TRUE)
         chrY[j] = mean(M[i_chrY, j] + U[i_chrY, j], na.rm = TRUE)
         auto[j] = mean(M[i_auto, j] + U[i_auto, j], na.rm = TRUE)
      }

      # Normalize allosomal intensities
      chrX = chrX / auto
      chrY = chrY / auto

      return(list(X = chrX, Y = chrY))
   })
}

#' @rdname check_sex
#' @param X,Y Forwarded from \code{check_sex()}
#' @param male,female Indices of male and female samples
#' @return \code{predict_sex} returns a factor with levels "f" (female) and "m" (male). 
#'   Samples that cannot be determined conclusively are assigned \code{NA}.
#' @export
predict_sex = function(X, Y, male, female){

   # Compute the robust Hodges-Lehmann estimator to determine the midpoint of X chromosome intensities
   # between known male and female samples.
   cutX = outer(X[male], X[female], "+")
   cutX = median(cutX) / 2

   # Likewise, compute the robust Hodges-Lehmann estimator for Y chromosome intensities
   cutY = outer(Y[male], Y[female], "+")
   cutY = median(cutY) / 2

   # Prediction is based on the quadrant (cutX/cutY) in which a sample falls.
   # Samples in the upper right (high X, high Y) and lower left (low X, low Y) quadrants 
   # are assigned NA (e.g. could represent Klinefelter syndrome, Turner syndrome, or failed assays).
   prediction = rep(NA, times = length(X))
   prediction[X >= cutX & Y <= cutY] =  "f"
   prediction[X <= cutX & Y >= cutY] =  "m"
   factor(prediction, levels=c("m", "f"), labels=c("m", "f"))
}