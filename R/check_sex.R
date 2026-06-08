#' Infer the sex of sample donor
#' 
#' Return the normalized average total intensities of probes targeting the X and Y chromosomes. 
#' 
#' @rdname check_sex
#' @param raw Output of calling \code{\link{read_idats()}}
#' @export
check_sex = function(raw) {
   
   if(!all(c('M', 'U', 'manifest') %in% names(raw))) stop('Invalid argument')

   with(raw, {

      # select allosomal and autosomal probes
      i_chrX = manifest[ forcats::fct_match(chr,   "chrX"),         index]
      i_chrY = manifest[ forcats::fct_match(chr,   "chrY"),         index]
      i_auto = manifest[!forcats::fct_match(chr, c("chrX","chrY")), index]

      J = ncol(M)
      chrX = numeric(J)
      chrY = numeric(J)
      auto = numeric(J)

      for (j in 1:J) {
         chrX[j] = mean(M[i_chrX, j] + U[i_chrX, j], na.rm = TRUE)
         chrY[j] = mean(M[i_chrY, j] + U[i_chrY, j], na.rm = TRUE)
         auto[j] = mean(M[i_auto, j] + U[i_auto, j], na.rm = TRUE)
      }

      # normalize allosomal intensities
      chrX = chrX / auto
      chrY = chrY / auto

      return(list(X = chrX, Y = chrY))
   })
}

#' @rdname check_sex
#' @param X,Y Forwarded from \code{check_sex()}
#' @param male,female Indices of male and female samples
#' @return \code{check_sex} returns the normalized average total intensities of probes targeting the X and Y chromosomes. These are forwarded to \code{predict_sex} which returns a factor with levels "f","m" (and \code{NA} if the sex cannot be determined conclusively).
#' @export
predict_sex = function(X, Y, male, female){

   # compute the robust Hodges-Lehmann estimator for the total intensity for X chr probes
   cutX = outer(X[male], X[female], "+")
   cutX = median(cutX) / 2

   # ... likewise for Y chr probes
   cutY = outer(Y[male], Y[female], "+")
   cutY = median(cutY) / 2

   # Prediction is based the quadrant (cutX/cutY) in which a sample falls
   # Samples in the upper right and lower left quadrant are assigned NA
   # (though there could be Klinefelter samples or similar)
   prediction = rep(NA, times = length(X))
   prediction[X >= cutX & Y <= cutY] =  "f"
   prediction[X <= cutX & Y >= cutY] =  "m"
   factor(prediction, levels=c("m", "f"), labels=c("m", "f"))
}