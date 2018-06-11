#' Compute relative agreement scores of SNP probes. Beta-values from the probes interrogating SNPs (59 for EPIC, 65 for 450K) are binned into three categories. Subsequently, conflicts are listed, either samples from supposedly same donor with a low agreement between SNP probes, or samples from supposedly different donors with a suspiciously high agreement.
#' @author Jonathan A. Heiss
#'
#' @title  Agreement of SNPs
#' @export
#'
#' @rdname check_snp_agreement
#' @param genotypes A list of length 3, each element a matrix with probablilities 
#' @param weights A matrix of weights (1-probability of being an outlier) for each SNP
#' @param donor_ids Vector of donor IDs, should include duplicates
#' @param sample_ids Vector of unique sample IDs
#' 
#' @return \code{check_snp_agreement} returns all connnected components of the conflict graph as list of data.tables. \code{agreement_} computes the same scores but returns a plot instead.
#' 
check_snp_agreement = function(genotypes,weights,donor_ids,sample_ids){

	stopifnot(anyDuplicated(sample_ids)==0)
	stopifnot(anyDuplicated(donor_ids)!=0)

	sample_ids = as.character(sample_ids)

	conflicts = cbind(CJ(donor_ids,donor_ids,sorted=FALSE),CJ(sample_ids,sample_ids,sorted=FALSE))
	setcolorder(conflicts,c(1,3,2,4))
	setnames(conflicts,1:4,c('donor1','sample1','donor2','sample2'))

	J = length(sample_ids)

	d = matrix(NA_real_,nrow=J,ncol=J)
	for(j in 1:J){
		tmp = 
		(weights * genotypes[[1]]) * (weights[,j] * genotypes[[1]][,j]) +
		(weights * genotypes[[2]]) * (weights[,j] * genotypes[[2]][,j]) +
		(weights * genotypes[[3]]) * (weights[,j] * genotypes[[3]][,j])

		d[,j] = colSums(tmp,na.rm=TRUE) / colSums(weights*weights[,j],na.rm=TRUE)

	}

	# first drop the duplicate entries and the diagonal
	d[upper.tri(d,diag=TRUE)] = NA_real_
	conflicts$similarity = as.numeric(d)
	rm(tmp,d)
	conflicts = conflicts[!is.na(similarity)]
	
	# drop cases of no interest
	conflicts = conflicts[!(donor1!=donor2 & similarity<0.90)]
	conflicts = conflicts[!(donor1==donor2 & similarity>0.90)]

	if(nrow(conflicts)==0) return(NULL)

	### find the components of the graph (consider list of conflicts as edges in a graph)
	
	# union/find data structure
	sets = unique(c(conflicts$sample1,conflicts$sample2))
	sets = hash::hash(keys=sets,values=sets)

	# union
	for(i in 1:nrow(conflicts)){ sets[[ conflicts[i]$sample1 ]] = sets [[ conflicts[i]$sample2 ]] }
	# find roots
	sets = sapply(hash::keys(sets),function(key){ while(key!=sets[[key]]){ key = sets[[key]] }; return(key) })
	sets = hash::hash(keys=names(sets),values=sets)
	conflicts[,group:=hash::values(sets,keys=sample1)]
	rm(sets)

	conflicts = split(conflicts,by='group',keep.by=FALSE)
	return(conflicts)
}

#' @rdname check_snp_agreement
#'
agreement_ = function(genotypes,weights,donor_ids,sample_ids,...){

	stopifnot(anyDuplicated(sample_ids)==0)
	stopifnot(anyDuplicated(donor_ids)!=0)

	conflicts = cbind(CJ(donor_ids,donor_ids,sorted=FALSE),CJ(sample_ids,sample_ids,sorted=FALSE))
	setcolorder(conflicts,c(1,3,2,4))
	setnames(conflicts,1:4,c('donor1','sample1','donor2','sample2'))

	J = length(sample_ids)

	d = matrix(NA_real_,nrow=J,ncol=J)
	for(j in 1:J){
		tmp = 
		(weights * genotypes[[1]]) * (weights[,j] * genotypes[[1]][,j]) +
		(weights * genotypes[[2]]) * (weights[,j] * genotypes[[2]][,j]) +
		(weights * genotypes[[3]]) * (weights[,j] * genotypes[[3]][,j])

		d[,j] = colSums(tmp,na.rm=TRUE) / colSums(weights*weights[,j],na.rm=TRUE)

	}

	# first drop the duplicate entries and the diagonal
	d[upper.tri(d,diag=TRUE)] = NA_real_
	conflicts$similarity = as.numeric(d)
	rm(tmp,d)
	conflicts = conflicts[!is.na(similarity)]

	boxplot(similarity ~ I(donor1==donor2),data=conflicts,horizontal=TRUE,ylim=c(0,1),xlab='Agreement',...)
	
	return(invisible(NULL))
}