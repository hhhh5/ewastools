#' @title Genotype calling
#' @description Detect SNP probes which do not fit into one of the three categories (AA, AB, BB).
#' A mixture model (3 Beta distributions, 1 uniform distribution for outliers) is fitted to all SNP
#' probes. After learning the model parameters via the EM algorithm, the probability of being an outlier
#' is computed for each SNP.
#'
#' @author Jonathan Heiss
#' @rdname call_genotypes
#'
#' @param snpmatrix Matrix of beta-values for SNP probes. Provide SNP probes as rows and samples as
#' columns. 
#' @param learn Logical. If \code{TRUE}, model parameters are learned from the dataset using the EM algorithm.
#' If \code{FALSE} (default), predefined model parameters are used.
#' @param genotypes Output of \code{call_genotypes}
#' @param maxiter Maximal number of iterations of the Expectation-Maximization algorithm learning
#' the mixture model.
#' 
#' @return For \code{call_genotypes()}, a list containing:
#' \item{snps}{The input \code{snpmatrix}.}
#' \item{outliers}{A matrix of a-posteriori probabilities of each SNP observation being an outlier.}
#' \item{gamma}{A list of three matrices containing the a-posteriori probabilities for each of the three genotypes (AA, AB, BB) respectively.}
#' \item{par}{Parameters of the fitted mixture model (class priors \code{pi}, beta parameters \code{shapes1} and \code{shapes2}, and outlier prior \code{alpha}).}
#' \item{loglik}{Log-likelihood in each iteration of the EM algorithm.}
#'
#' @return For \code{snp_outliers()}, a metric assessing the outlierness of the SNP beta-values.
#' 	High values may indicate either contaminated or failed samples.
#'
#' @return For \code{mxm_()}, a histogram showing the distribution of beta-values for SNP probes with
#' the density function of the mixture model overlaid.
#'
#' @export
#'
call_genotypes = function(snpmatrix, learn=FALSE, maxiter = 50)
{
   snps = snpmatrix
   dim(snps) = NULL

   # Drop NAs to be able to compute likelihoods, but keep
   # track of which indices are NA to re-insert them later
   NAs = which(is.na(snps))
   snps = na.omit(snps)
   n = length(snps)

   if(learn==FALSE){

      # Use predefined model parameters
      # (might work better if training set is small or contains many outliers)

      # Class probability for outliers
      alpha = 0.06646095

      # Class probabilities for the three genotypes: homozygous AA, heterozygous AB, homozygous BB
      pi = c(0.2818387,0.4330363,0.2851250)

      # Beta distribution parameters
      shapes1 = c(2.206479,80.830012,40.640821)
      shapes2 = c(38.043029,84.411900,3.315509)

      # Uniform distribution representing outliers
      p = 1

      loglik = NULL

      # Expectation step: Compute posterior probabilities (gamma) for each of the three genotype classes
      gamma = cbind(
          pi[1] * dbeta(snps,shape1=shapes1[1],shape2=shapes2[1])
         ,pi[2] * dbeta(snps,shape1=shapes1[2],shape2=shapes2[2])
         ,pi[3] * dbeta(snps,shape1=shapes1[3],shape2=shapes2[3])
         )

      # Scale by non-outlier prior probability
      gamma = (1-alpha) * gamma
      tmp = rowSums(gamma)
      gamma = gamma/tmp

      # Compute a-posteriori probability of each observation being an outlier
      outliers = (alpha*p) / ((alpha*p) + tmp)

   } else { # Learn dataset-specific model parameters using the EM algorithm

      alpha = 1e-2
      outliers = rep(alpha,times=n)
      pi = c(1/3,1/3,1/3) # Initial class probabilities for the three genotypes
      shapes1 = c(10,80,80) # Initial shape1 values
      shapes2 = c(80,80,10) # Initial shape2 values
      p = 1 # Uniform distribution for outliers
      gamma = NA

      # Expectation Step (E-step)
      # Computes the posterior probabilities for genotype classes and outlier status
      e_step = function(){
         gamma = cbind(
             pi[1] * dbeta(snps,shape1=shapes1[1],shape2=shapes2[1])
            ,pi[2] * dbeta(snps,shape1=shapes1[2],shape2=shapes2[2])
            ,pi[3] * dbeta(snps,shape1=shapes1[3],shape2=shapes2[3])
            )

         gamma = (1-alpha) * gamma
         tmp = rowSums(gamma)
         gamma <<- gamma/tmp

         outliers <<- (alpha*p) / ((alpha*p) + tmp)

         # Calculate overall log-likelihood
         loglik = (alpha*p) + tmp
         loglik = sum(log(loglik))

         return(loglik)
      }

      # Maximization Step (M-step)
      # Estimates the model parameters (class priors, beta shapes) that maximize the expected log-likelihood
      m_step = function(){

         gamma = gamma * (1-outliers)

         # Maximum Likelihood Estimation (MLE) of Beta distribution parameters for each class
         s1 = eBeta(snps,gamma[,1])
         s2 = eBeta(snps,gamma[,2])
         s3 = eBeta(snps,gamma[,3])

         shapes1 <<- c(s1$shape1,s2$shape1,s3$shape1)
         shapes2 <<- c(s1$shape2,s2$shape2,s3$shape2)

         # MLE of class priors
         pi = apply(gamma,2,sum)
         pi <<- pi/sum(pi)
         alpha <<- sum(outliers)/n

         invisible(NULL)
      }

      loglik = rep(NA_real_,times=maxiter)
      loglik[1] = e_step()

      i = 2; gain=Inf;

      while (i < maxiter & gain > 1e-4) { # stop if maxiter reached or improvement is below threshold
         m_step()
         loglik[i] = e_step()
         gain = loglik[i]-loglik[i-1]
         i = i + 1
      }

      loglik=loglik[1:(i-1)]
   }

   ## Re-insert missing values (NAs) that were dropped initially
   if (length(NAs) != 0) {
      tmp = rep(NA_real_,times=length(snpmatrix))
      tmp[-NAs] = outliers
      outliers = tmp
   }

   dim(outliers) = dim(snpmatrix)

   gamma = lapply(1:3,function(k){

      if (length(NAs) != 0) {
         tmp = rep(NA_real_, times = length(snpmatrix))
         tmp[-NAs] = gamma[, k]
      } else {
         tmp = gamma[, k]
      }

      dim(tmp) = dim(snpmatrix)
      tmp
   })


   return(list(
       snps=snpmatrix
      ,outliers=outliers
      ,gamma=gamma
      ,par=list(pi=pi,shapes1=shapes1,shapes2=shapes2,alpha=alpha)
      ,loglik=loglik
   ))
}

#' @rdname call_genotypes
#'
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
   stopifnot("outliers" %in% names(genotypes))

   # Average log odds of beta-values being outliers across all SNP probes
   log_odds = genotypes$outliers / (1-genotypes$outliers)
   log_odds = colMeans(log2(log_odds), na.rm = TRUE) # NAs stem from missing observations
   log_odds
}

#' @rdname call_genotypes
#' @param x Vector of observations (beta-values).
#' @param w Vector of weights (posterior probabilities for the class).
#' @return A list with estimated Beta distribution parameters:
#' \item{shape1}{Estimated shape1 parameter.}
#' \item{shape2}{Estimated shape2 parameter.}
#'
eBeta = function(x,w)
{
   # Method of moments estimation of Beta distribution parameters from weighted observations
   n = length(w)
   w = n*w/sum(w)
   sample.mean =  mean(w*x)
   sample.var  = (mean(w*x^2)-sample.mean^2) * n/(n-1)
   v = sample.mean * (1-sample.mean)
    
   if (sample.var < v){
      shape1 = sample.mean * (v/sample.var - 1)
      shape2 = (1 - sample.mean) * (v/sample.var - 1)
   } else {
      shape2 = sample.mean * (v/sample.var - 1)
      shape1 = (1 - sample.mean) * (v/sample.var - 1)
   }
    
   list(shape1 = shape1, shape2 = shape2)
}