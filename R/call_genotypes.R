#' @title Genotype calling
#' @description Detect SNP probes which do not fit into on of the three categories (AA,AB,BB).
#' A mixture model (3 normal distribution, 1 uniform distribution for outliers) is fitted to all SNP probes.
#' In order to achieve a more normal distribution beta-values are transformed to M-values first.
#' After learning the model parameters via EM algorithm, the probability of being an outlier is computed for each SNP.
#'
#' @author Jonathan A. Heiss
#' @rdname call_genotypes
#'
#' @param snpmatrix Matrix of beta-values for SNP probes. Provide SNPs probes as rows and samples as columns. 
#' @param mxm Mixture model, output of \code{call_genotypes}
#' @param maxiter Maximal number of iterations of the Expectation-Maximization algorithm learning the mixture model
#' 
#' @return  For \code{call_genotypes}, a list containing
#' \item{par}{parameters of the mixture model}
#' \item{loglik}{Log-likelihood in each iteration of the EM algorithm}
#' \item{outliers}{A-posteriori probability of SNP being an outlier}
#' \item{gamma}{A-posteriori probabilities for each of the three genotypes}
#' @return For \code{mxm_}, a histogram showing the distribution of beta-values for SNP probes with the density function of the mixture model overlaid.
#' @export
#'
call_genotypes <- function(snpmatrix,maxiter=50){

	snps = snpmatrix
	dim(snps) = NULL
	NAs = which(is.na(snps))
	snps = na.omit(snps)
	n = length(snps)
	alpha = 1e-2
	outliers = rep(alpha,n)
	pi = c(1/3,1/3,1/3)
	shapes1 = c(10,80,80) 
	shapes2 = c(80,80,10)
	p = 1
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

	loglik = rep(NA_real_,maxiter)
	loglik[1] = e_step()

	i = 2; gain=Inf;

	while(i<maxiter & gain>1e-4){
		m_step()
		loglik[i] = e_step()
		gain = loglik[i]-loglik[i-1]
		i=i+1
	}

	if(length(NAs)!=0){
		tmp = rep(NA,length(snpmatrix))
		tmp[-NAs] = outliers
	}else{
		tmp = outliers
	}
	dim(tmp) = dim(snpmatrix)

	tmp = rep(NA,length(snpmatrix))
	if(length(NAs)!=0){ tmp[-NAs] = outliers } else tmp = outliers
	dim(tmp) = dim(snpmatrix)

	gamma = lapply(1:3,function(k){

		if(length(NAs)!=0){
			tmp = rep(NA,length(snpmatrix))
			tmp[-NAs] = gamma[,k]
		}else{
			tmp = gamma[,k]
		}

		dim(tmp) = dim(snpmatrix)
		tmp
	})

	return(list(
		snps=snpmatrix
		,outliers=tmp
		,gamma=gamma
		,par=list(pi=pi,shapes1=shapes1,shapes2=shapes2,alpha=alpha)
		,loglik=loglik[1:(i-1)]
		))

}

#' @rdname call_genotypes
#' @export
#'
mxm_ = function(mxm){
	
	with(mxm,{

		hist(snps,breaks=200,freq=FALSE,xlab='β-value',main=NA)
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
snp_outliers = function(gt){

 	if(!"outliers"%in%names(gt)) stop('Invalid argument')

 	log_odds = gt$outliers / (1-gt$outliers)
 	log_odds = colMeans(log2(log_odds),na.rm=TRUE)
}

#' @rdname call_genotypes
#'
eBeta = function(x,w){
	
	n = sum(w)
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
