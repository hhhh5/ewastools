#' @title Genotype calling
#' @description Detect SNP probes which do not fit into on of the three categories (AA,AB,BB).
#' A mixture model (3 Beta distributions, 1 uniform distribution for outliers) is fitted to all SNP probes.
#' After learning the model parameters via EM algorithm, the probability of being an outlier is computed for each SNP.
#'
#' @author Jonathan A. Heiss
#' @rdname call_genotypes
#'
#' @param snpmatrix Matrix of beta-values for SNP probes. Provide SNPs probes as rows and samples as columns. 
#' @param genotypes Output of \code{call_genotypes}
#' @param maxiter Maximal number of iterations of the Expectation-Maximization algorithm learning the mixture model
#' 
#' @return  For \code{call_genotypes}, a list containing
#' \item{par}{Parameters of the mixture model}
#' \item{loglik}{Log-likelihood in each iteration of the EM algorithm}
#' \item{outliers}{A-posteriori probability of SNP being an outlier}
#' \item{gamma}{A-posteriori probabilities for each of the three genotypes}
#' @return For \code{snp_outliers}, a metric assessing the outlierness of the SNP beta-values. High values may indicate either contaminated or failed samples.
#' @return For \code{mxm_}, a histogram showing the distribution of beta-values for SNP probes with the density function of the mixture model overlaid.
#' @export
#'
call_genotypes <- function(snpmatrix,learn=FALSE,maxiter=50){

	snps = snpmatrix
	dim(snps) = NULL

	# Drop NAs to be able to compute likelihoods, but keep
	# score of which entries to reinsert them again later
	NAs = which(is.na(snps))
	snps = na.omit(snps)
	n = length(snps)

	if(learn==FALSE){
		
		# Use predefined model parameters
		# (might work better if training set is small or contains many outliers)

		# Class probability for outliers
		alpha = 0.06646095
		
		# Class probabilities for homozygous and heterozygous genotypes
		pi = c(0.2818387,0.4330363,0.2851250)

		# Beta distribution parameters
		shapes1 = c(2.206479,80.830012,40.640821)
		shapes2 = c(38.043029,84.411900,3.315509)

		# Uniform distribution representing outliers
		p = 1
		
		loglik = NULL

		gamma = cbind(
			 pi[1] * dbeta(snps,shape1=shapes1[1],shape2=shapes2[1])
			,pi[2] * dbeta(snps,shape1=shapes1[2],shape2=shapes2[2])
			,pi[3] * dbeta(snps,shape1=shapes1[3],shape2=shapes2[3])
			)

		gamma = (1-alpha) * gamma
		tmp = rowSums(gamma)
		gamma = gamma/tmp
	 
		outliers = (alpha*p) / ((alpha*p) + tmp)

	}else{ # Learn dataset-specific model parameters using the EM algorithm

		alpha = 1e-2
		outliers = rep(alpha,times=n)
		pi = c(1/3,1/3,1/3) # Class probabilities
		shapes1 = c(10,80,80) 
		shapes2 = c(80,80,10)
		p = 1 # Uniform distribution
		gamma = NA
		
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

			loglik = (alpha*p) + tmp
			loglik = sum(log(loglik))

			return(loglik)
		}

		m_step = function(){

			gamma = gamma * (1-outliers)

			# MLE
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

		while(i<maxiter & gain>1e-4){ # stop if maxiter reached or improvement is below threshold
			m_step()
			loglik[i] = e_step()
			gain = loglik[i]-loglik[i-1]
			i=i+1
		}

		loglik=loglik[1:(i-1)]

	}
	
	## Re-insert missing values
	if(length(NAs)!=0){
		tmp = rep(NA_real_,times=length(snpmatrix))
		tmp[-NAs] = outliers
		outliers = tmp
	}
	
	dim(outliers) = dim(snpmatrix)

	gamma = lapply(1:3,function(k){

		if(length(NAs)!=0){
			tmp = rep(NA_real_,times=length(snpmatrix))
			tmp[-NAs] = gamma[,k]
		}else{
			tmp = gamma[,k]
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
mxm_ = function(genotypes){
	
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
snp_outliers = function(genotypes){

 	if(!"outliers"%in%names(genotypes)) stop('Invalid argument')

 	# Average log odds of beta-values being outliers across all SNP probes
 	log_odds = genotypes$outliers / (1-genotypes$outliers)
 	log_odds = colMeans(log2(log_odds),na.rm=TRUE)
 	log_odds
}

#' @rdname call_genotypes
#'
eBeta = function(x,w){
	
	# Beta distribution parameter estimation
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