#' Check agreement of genetic fingerprint between samples
#'
#' These functions check the agreement of genetic fingerprints between samples
#' (after genotype calling with \code{\link{call_genotypes}}). In case of \code{check_snp_agreement},
#' the user supplies for each sample the assumed donor ID and the function returns all conflicts,
#' meaning either samples that are supposed to come from the same donor but have distinct genetic
#' fingerprints, or samples from supposedly different donors with the same fingerprint.
#' \code{enumerate_sample_donors} on the other hand infers donor IDs in case they are unknown
#' (useful, for example, to detect technical replicates in a public dataset).
#'
#' Two samples are considered to have the same genetic fingerprint if their agreement score
#' (the probability of sharing the same genotype, averaged across SNP probes and weighted by the
#' probability of the observations not being outliers) exceeds 0.9.
#'
#' @author Jonathan Heiss
#'
#' @export
#'
#' @rdname check_snp_agreement
#' @param genotypes Output of \code{\link{call_genotypes}}.
#' @param donor_ids Vector of donor IDs corresponding to the samples (columns of the SNP matrix).
#' Contains duplicates if several samples come from the same donor.
#' @param sample_ids Vector of unique sample IDs corresponding to the samples.
#'
#' @return \code{check_snp_agreement} returns \code{NULL} if no conflicts were found. Otherwise a list
#' of \code{data.table}s, each representing one connected component of the conflict graph (i.e.,
#' a group of related conflicts), with columns \code{donor1}, \code{sample1}, \code{donor2},
#' \code{sample2}, and \code{agreement}.
#' @return \code{enumerate_sample_donors} returns an integer vector of inferred donor IDs, one for
#' each sample. Samples with the same ID share the same genetic fingerprint.
#'
check_snp_agreement = function(genotypes, donor_ids, sample_ids) {

   J = ncol(genotypes$snps)

   stopifnot(anyDuplicated(sample_ids) == 0) # Stop if sample IDs are not unique
   stopifnot(length(sample_ids) == J)

   ## Generate all pairwise combinations of samples/donors
   conflicts = cbind(
      CJ(donor_ids,  donor_ids,  sorted = FALSE),
      CJ(sample_ids, sample_ids, sorted = FALSE))
   setcolorder(conflicts, c(1, 3, 2, 4)) # Rearrange columns
   setnames(conflicts, 1:4, c('donor1','sample1','donor2','sample2'))

   ## Compute the agreement between all sample pairs
   d = snp_agreement_matrix_(genotypes)
   conflicts$agreement = as.numeric(d)
   rm(d)
   conflicts = conflicts[!is.na(agreement)]
   
   ## Drop expected/non-conflict cases:
   # - Samples from different donors expected to have low agreement (< 0.90)
   # - Samples from the same donor expected to have high agreement (> 0.90)
   conflicts = conflicts[!(donor1!=donor2 & agreement<0.90)]
   conflicts = conflicts[!(donor1==donor2 & agreement>0.90)]

   if(nrow(conflicts)==0) return(NULL)
   
   ## Find the weakly connected components of the graph
   # (consider conflicts as edges in a graph to group related issues together)
   e = rep(NA_character_, times = 2 * nrow(conflicts))
   e[c(TRUE,FALSE)] = as.character(conflicts$sample1)
   e[c(FALSE,TRUE)] = as.character(conflicts$sample2)

   g = igraph::make_graph(edges=e,directed=FALSE)
   g = igraph::components(g,mode="weak")$membership

   conflicts$group = g[as.character(conflicts$sample1)]
   conflicts = split(conflicts,by='group',keep.by=FALSE)
   return(conflicts)
}

#' Boxplot of agreement scores for pairs of samples from the same vs. different donors
#'
#' @param ... Additional arguments passed to \code{\link[graphics]{boxplot}}.
#' @noRd
agreement_ = function(genotypes,donor_ids,sample_ids,...){

   J = ncol(genotypes$snps)
   
   stopifnot(anyDuplicated(sample_ids)==0) # Stop if sample IDs are not unique
   stopifnot(length(sample_ids) == J)

   conflicts = cbind(
      CJ(donor_ids,  donor_ids,  sorted = FALSE),
      CJ(sample_ids, sample_ids, sorted = FALSE))
   setcolorder(conflicts, c(1,3,2,4))
   setnames(conflicts, 1:4, c('donor1','sample1','donor2','sample2'))

   ## Compute the agreement between all sample pairs
   d = snp_agreement_matrix_(genotypes)
   conflicts$agreement = as.numeric(d)
   rm(d)
   conflicts = conflicts[!is.na(agreement)]

   # Plot a boxplot comparing agreement distribution for same vs. different donors
   boxplot(agreement ~ I(donor1==donor2),data=conflicts,horizontal=TRUE,ylim=c(0,1),xlab='Agreement',...)
   
   return(invisible(NULL))
}

#' @rdname check_snp_agreement
#' @export
#' 
enumerate_sample_donors = function(genotypes){

   J = ncol(genotypes$snps)

   samesame = CJ(1:J,1:J,sorted=FALSE)
   setnames(samesame,1:2,c('sample1','sample2'))

   ## Compute the agreement between all sample pairs
   d = snp_agreement_matrix_(genotypes)
   samesame$agreement = as.numeric(d)
   rm(d)
   samesame = samesame[!is.na(agreement)]
   
   ## Select sample pairs with the same fingerprint/high agreement (> 0.90)
   samesame = samesame[agreement>0.90]

   # If there are no highly agreeing pairs, every sample belongs to a unique donor
   if(nrow(samesame)==0) return(1:J)

   ## Find the strongly connected components of the graph (consider samesame pairs as edges in a graph)
   e = rep(NA,times=2*nrow(samesame))
   e[c(TRUE,FALSE)] = samesame$sample1
   e[c(FALSE,TRUE)] = samesame$sample2

   g = igraph::make_graph(edges=e,n=J,directed=FALSE)
   return(igraph::components(g,mode="strong")$membership)
}

#' Compute lower-triangular matrix of pairwise SNP agreement scores
#' @noRd
snp_agreement_matrix_ = function(genotypes) {
   weights = 1 - genotypes$outliers # we don't want outliers to be counted in the genetic fingerprint
   weights[is.na(weights)] = 0

   w_gamma = lapply(genotypes$gamma, function(g) {
      wg = weights * g
      wg[is.na(wg)] = 0
      wg
   })

   # Probability of sharing the same genotype (sum of product of genotype posteriors)
   # weighted by the probability of neither observation being an outlier,
   # normalized by the sum of joint outlier weights.
   num = crossprod(w_gamma[[1]]) + crossprod(w_gamma[[2]]) + crossprod(w_gamma[[3]])
   den = crossprod(weights)
   d = num / den

   ## Drop duplicate entries (since matrix is symmetric) and the diagonal
   d[upper.tri(d, diag = TRUE)] = NA_real_
   d
}
