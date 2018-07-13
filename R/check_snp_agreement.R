#' These function check the agreement of genetic fingerprints between samples (after genotype calling with \code{call_genotypes}). In case of \code{check_snp_agreement}, the user supplies for each sample the assumed donor ID and the function returns all conflicts (i.d., either samples that are supposed to come from the same donor but have distinct genetic fingerprints or samples from supposedly different donors with the same fingerprints). \code{check_snp_agreement} on the other hand is inferring donor IDs in case they are unknown (useful for example to detect technical replicates in a public dataset).
#' 
#' @author Jonathan A. Heiss
#' 
#' @title  Agreement of SNPs
#' @export
#'
#' @rdname check_snp_agreement
#' @param genotypes Output of \code{call_genotypes}
#' @param donor_ids Vector of donor IDs, should include duplicates
#' @param sample_ids Vector of unique sample IDs
#' 
#' @return \code{check_snp_agreement} returns all conflicts found between user-supplied donor IDs and those derived from the genetic fingerprint. Either \code{NULL} if not conflicts were found, or a list of \code{data.table}s, each representing one connected component of the conflict graph.
#' @return \code{enumerate_sample_donors} returns for each sample the donor IDs.
#' 
check_snp_agreement = function(genotypes,donor_ids,sample_ids){

	J = ncol(genotypes$snps)
	weights = 1-genotypes$outliers
	gamma = genotypes$gamma

	stopifnot(anyDuplicated(sample_ids)==0)
	stopifnot(length(sample_ids) == J)

	sample_ids = as.character(sample_ids)

	conflicts = cbind(CJ(donor_ids,donor_ids,sorted=FALSE),CJ(sample_ids,sample_ids,sorted=FALSE))
	setcolorder(conflicts,c(1,3,2,4))
	setnames(conflicts,1:4,c('donor1','sample1','donor2','sample2'))

	d = matrix(NA_real_,nrow=J,ncol=J)
	for(j in 1:J){
		tmp = 
		(weights * gamma[[1]]) * (weights[,j] * gamma[[1]][,j]) +
		(weights * gamma[[2]]) * (weights[,j] * gamma[[2]][,j]) +
		(weights * gamma[[3]]) * (weights[,j] * gamma[[3]][,j])

		d[,j] = colSums(tmp,na.rm=TRUE) / colSums(weights*weights[,j],na.rm=TRUE)

	}

	# first drop the duplicate entries and the diagonal
	d[upper.tri(d,diag=TRUE)] = NA_real_
	conflicts$agreement = as.numeric(d)
	rm(tmp,d)
	conflicts = conflicts[!is.na(agreement)]
	
	# drop cases of no interest
	conflicts = conflicts[!(donor1!=donor2 & agreement<0.90)]
	conflicts = conflicts[!(donor1==donor2 & agreement>0.90)]

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
agreement_ = function(genotypes,donor_ids,sample_ids,...){

	J = ncol(genotypes$snps)
	weights = 1-genotypes$outliers
	gamma = genotypes$gamma
	
	stopifnot(anyDuplicated(sample_ids)==0)
	stopifnot(length(sample_ids) == J)

	conflicts = cbind(CJ(donor_ids,donor_ids,sorted=FALSE),CJ(sample_ids,sample_ids,sorted=FALSE))
	setcolorder(conflicts,c(1,3,2,4))
	setnames(conflicts,1:4,c('donor1','sample1','donor2','sample2'))


	d = matrix(NA_real_,nrow=J,ncol=J)
	for(j in 1:J){
		tmp = 
		(weights * gamma[[1]]) * (weights[,j] * gamma[[1]][,j]) +
		(weights * gamma[[2]]) * (weights[,j] * gamma[[2]][,j]) +
		(weights * gamma[[3]]) * (weights[,j] * gamma[[3]][,j])

		d[,j] = colSums(tmp,na.rm=TRUE) / colSums(weights*weights[,j],na.rm=TRUE)

	}

	# first drop the duplicate entries and the diagonal
	d[upper.tri(d,diag=TRUE)] = NA_real_
	conflicts$agreement = as.numeric(d)
	rm(tmp,d)
	conflicts = conflicts[!is.na(agreement)]

	boxplot(agreement ~ I(donor1==donor2),data=conflicts,horizontal=TRUE,ylim=c(0,1),xlab='Agreement',...)
	
	return(invisible(NULL))
}

#' @rdname check_snp_agreement
#' @export
#' 
enumerate_sample_donors = function(genotypes){

	J = ncol(genotypes$snps)
	weights = 1-genotypes$outliers
	gamma = genotypes$gamma

	samesame = CJ(1:J,1:J,sorted=FALSE)
	setnames(samesame,1:2,c('sample1','sample2'))

	d = matrix(NA_real_,nrow=J,ncol=J)
	for(j in 1:J){
		tmp = 
		(weights * gamma[[1]]) * (weights[,j] * gamma[[1]][,j]) +
		(weights * gamma[[2]]) * (weights[,j] * gamma[[2]][,j]) +
		(weights * gamma[[3]]) * (weights[,j] * gamma[[3]][,j])

		d[,j] = colSums(tmp,na.rm=TRUE) / colSums(weights*weights[,j],na.rm=TRUE)

	}

	# first drop the duplicate entries and the diagonal
	d[upper.tri(d,diag=TRUE)] = NA_real_
	samesame$agreement = as.numeric(d)
	rm(tmp,d)
	samesame = samesame[!is.na(agreement)]
	
	# drop cases of no interest
	samesame = samesame[agreement>0.85]

	if(nrow(samesame)==0) return(1:J)

	### find the strongly connected components of the graph (consider samesame as edges in a graph)
	e = rep(NA,times=2*nrow(samesame))
	e[c(TRUE,FALSE)] = samesame$sample1
	e[c(FALSE,TRUE)] = samesame$sample2

	g = igraph::make_graph(edges=e,n=J,directed=FALSE)
	return(igraph::components(g,mode="strong")$membership)
}
