#' @title Genotype calling
#' @description \code{call_genotypes()} calls genotypes from the SNP probes on the array and
#' detects outlying observations. A mixture model (three Beta distributions representing the
#' genotypes AA, AB and BB, and one uniform distribution representing outliers) is fitted to the
#' beta-values of all SNP probes. Either predefined model parameters are used (default), or the
#' parameters are learned from the data via the EM algorithm. Subsequently, the a-posteriori
#' probabilities of each genotype and of being an outlier are computed for every SNP observation.
#'
#' \code{snp_outliers()} summarizes the outlier probabilities per sample as the average log2 odds
#' of being an outlier across all SNP probes. Samples with a score greater than -4 should be flagged
#' as they are likely contaminated or failed assays.
#'
#' @author Jonathan Heiss
#' @rdname call_genotypes
#'
#' @param snpmatrix Matrix of beta-values for SNP probes. Provide SNP probes as rows and samples as
#' columns.
#' @param learn Logical. If \code{TRUE}, model parameters are learned from the dataset using the EM algorithm.
#' If \code{FALSE} (default), predefined model parameters are used. The predefined parameters are
#' recommended for small datasets or datasets containing many outliers.
#' @param maxiter Maximum number of iterations of the EM algorithm. Only used if \code{learn = TRUE}.
#' @param genotypes Output of \code{call_genotypes()}.
#'
#' @return For \code{call_genotypes()}, a list containing:
#' \item{snps}{The input \code{snpmatrix}.}
#' \item{outliers}{A matrix of a-posteriori probabilities of each SNP observation being an outlier.}
#' \item{gamma}{A list of three matrices containing the a-posteriori probabilities for each of the
#'    three genotypes (AA, AB, BB) respectively.}
#' \item{par}{Parameters of the mixture model (class priors \code{pi}, Beta distribution parameters \code{shapes1} and \code{shapes2}, and outlier prior \code{alpha}).}
#' \item{loglik}{Log-likelihood in each iteration of the EM algorithm (\code{NULL} if \code{learn = FALSE}).}
#'
#' @return For \code{snp_outliers()}, a numeric vector with one outlier score per sample
#' (average log2 odds of SNP observations being outliers). High values indicate contaminated or failed samples.
#'
#' @seealso \code{\link{check_snp_agreement}} and \code{\link{enumerate_sample_donors}} for using the
#' called genotypes to track samples.
#'
#' @export
#'
call_genotypes = function(snpmatrix, learn=FALSE, maxiter = 50)
{
   checkmate::assert_matrix(snpmatrix, mode = "numeric")
   checkmate::assert_flag(learn)
   checkmate::assert_count(maxiter, positive = TRUE)

   snps = as.vector(snpmatrix)

   # Drop NAs to be able to compute likelihoods, but keep
   # track of which indices are NA to re-insert them later
   NAs = which(is.na(snps))
   snps = na.omit(snps)

   # Expectation Step (E-step)
   # Computes posterior probabilities for genotype classes and outlier status, plus log-likelihood
   e_step = function(snps, par) {
      p_outlier = 1 # Uniform(0, 1) density representing outliers

      gamma = cbind(
          par$pi[1] * dbeta(snps, shape1 = par$shapes1[1], shape2 = par$shapes2[1])
         ,par$pi[2] * dbeta(snps, shape1 = par$shapes1[2], shape2 = par$shapes2[2])
         ,par$pi[3] * dbeta(snps, shape1 = par$shapes1[3], shape2 = par$shapes2[3])
      )

      # Scale by non-outlier prior probability.
      # Normalizing by `non_outlier_lik` yields P(genotype = k | not outlier), which sums to 1 across k.
      gamma = (1 - par$alpha) * gamma
      non_outlier_lik = rowSums(gamma)
      gamma = gamma / non_outlier_lik

      # A-posteriori probability of each observation being an outlier
      outliers = (par$alpha * p_outlier) / ((par$alpha * p_outlier) + non_outlier_lik)

      # Overall log-likelihood
      loglik = sum(log((par$alpha * p_outlier) + non_outlier_lik))

      list(gamma = gamma, outliers = outliers, loglik = loglik)
   }

   # Maximization Step (M-step)
   # Updates model parameters (class priors, Beta shapes, outlier prior) from current posteriors
   m_step = function(snps, post) {
      # Joint posterior P(genotype = k, not outlier | data) downweights outliers during parameter updates
      weights = post$gamma * (1 - post$outliers)

      # Method of moments estimation of Beta distribution parameters for each class
      s1 = eBeta(snps, weights[, 1])
      s2 = eBeta(snps, weights[, 2])
      s3 = eBeta(snps, weights[, 3])

      # MLE of class priors and outlier proportion
      class_priors = colSums(weights)
      class_priors = class_priors / sum(class_priors)

      list(
          pi      = class_priors
         ,shapes1 = c(s1$shape1, s2$shape1, s3$shape1)
         ,shapes2 = c(s1$shape2, s2$shape2, s3$shape2)
         ,alpha   = mean(post$outliers)
      )
   }

   if (!learn) {

      # Use predefined model parameters
      # (might work better if training set is small or contains many outliers)
      par = list(
          pi      = c(0.2818387, 0.4330363, 0.2851250) # homozygous AA, heterozygous AB, homozygous BB
         ,shapes1 = c(2.206479, 80.830012, 40.640821)
         ,shapes2 = c(38.043029, 84.411900, 3.315509)
         ,alpha   = 0.06646095 # Class probability for outliers
      )

      post = e_step(snps, par)
      loglik = NULL

   } else { # Learn dataset-specific model parameters using the EM algorithm

      par = list(
          pi      = c(1/3, 1/3, 1/3) # Initial class probabilities for the three genotypes
         ,shapes1 = c(10, 80, 80)    # Initial shape1 values
         ,shapes2 = c(80, 80, 10)    # Initial shape2 values
         ,alpha   = 1e-2             # Initial outlier prior probability
      )

      loglik = rep(NA_real_, times = maxiter)
      post = e_step(snps, par)
      loglik[1] = post$loglik

      i = 2
      gain = Inf

      while (i < maxiter && gain > 1e-4) { # stop if maxiter reached or improvement is below threshold
         par = m_step(snps, post)
         post = e_step(snps, par)
         loglik[i] = post$loglik
         gain = loglik[i] - loglik[i - 1]
         i = i + 1
      }

      loglik = loglik[seq_len(i - 1)]
   }

   ## Re-insert missing values (NAs) that were dropped initially
   outliers = post$outliers
   if (length(NAs) != 0) {
      tmp = rep(NA_real_, times = length(snpmatrix))
      tmp[-NAs] = outliers
      outliers = tmp
   }

   dim(outliers) = dim(snpmatrix)

   gamma = lapply(1:3, function(k) {

      if (length(NAs) != 0) {
         tmp = rep(NA_real_, times = length(snpmatrix))
         tmp[-NAs] = post$gamma[, k]
      } else {
         tmp = post$gamma[, k]
      }

      dim(tmp) = dim(snpmatrix)
      tmp
   })


   return(list(
       snps = snpmatrix
      ,outliers = outliers
      ,gamma = gamma
      ,par = par
      ,loglik = loglik
   ))
}

#' Plot a histogram of SNP beta-values with the density of the fitted mixture model overlaid
#'
#' @param genotypes Output of \code{call_genotypes()}.
#' @noRd
mxm_ = function(genotypes)
{
   # Visualizes the distribution of SNP beta-values and overlays the fitted mixture model density
   with(genotypes,{

   hist(snps,breaks=200,freq=FALSE,xlab="beta-value",main=NA)
      curve(  
         par$alpha+(1-par$alpha)*(
         par$pi[1]*dbeta(x,shape1=par$shapes1[1],shape2=par$shapes2[1])+
         par$pi[2]*dbeta(x,shape1=par$shapes1[2],shape2=par$shapes2[2])+
         par$pi[3]*dbeta(x,shape1=par$shapes1[3],shape2=par$shapes2[3]))
         ,from=0,to=1,col=2,add=TRUE,lwd=2)

      return(invisible(NULL))

   })
}

#' @rdname call_genotypes
#' @export
#'
snp_outliers = function(genotypes)
{
   checkmate::assert_list(genotypes)
   checkmate::assert_names(names(genotypes), must.include = "outliers")

   # Average log odds of beta-values being outliers across all SNP probes
   log_odds = genotypes$outliers / (1-genotypes$outliers)
   log_odds = colMeans(log2(log_odds), na.rm = TRUE) # NAs stem from missing observations
   log_odds
}

#' Method of moments estimation of Beta distribution parameters from weighted observations
#'
#' @param x Vector of observations (beta-values).
#' @param w Vector of weights (posterior probabilities for the class).
#' @return A list with estimated Beta distribution parameters:
#' \item{shape1}{Estimated shape1 parameter.}
#' \item{shape2}{Estimated shape2 parameter.}
#' @noRd
eBeta = function(x, w)
{
  stopifnot(all(x > 0 & x < 1), all(w >= 0), length(x) == length(w))
  
  # Normalize weights
  w = w / sum(w)
  
  # Weighted mean and variance
  x_bar <- sum(w * x)
  s2 <- sum(w * (x - x_bar)^2)
  
  # Condition check: sample variance must be strictly less than x_bar * (1 - x_bar)
  term <- (x_bar * (1 - x_bar) / s2) - 1
  if (term <= 0) {
    stop("Weighted variance is too large to fit a valid Beta distribution.")
  }
  
  alpha <- x_bar * term
  beta  <- (1 - x_bar) * term
  
  return(list(shape1 = alpha, shape2 = beta))
}