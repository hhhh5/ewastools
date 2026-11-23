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

   checkmate::assert_list(genotypes)
   checkmate::assert_names(names(genotypes), must.include = c("snps", "outliers", "gamma"))

   J = ncol(genotypes$snps)

   checkmate::assert_atomic(donor_ids, len = J)
   checkmate::assert_atomic(sample_ids, len = J, unique = TRUE)

   AGREEMENT_THRESHOLD = 0.90

   pairs = pairwise_snp_agreement_(genotypes, donor_ids = donor_ids, sample_ids = sample_ids)

   # Keep pairs that conflict with recorded donor labels:
   # - Different recorded donors with high genetic agreement (>= 0.90)
   # - Same recorded donor with low genetic agreement (<= 0.90)
   conflicts = pairs[
      (donor1 != donor2 & agreement >= AGREEMENT_THRESHOLD) |
      (donor1 == donor2 & agreement <= AGREEMENT_THRESHOLD)
   ]

   if (nrow(conflicts) == 0) return(NULL)
   
   # Group related conflicts into weakly connected components of the conflict graph
   edges = data.frame(
      from = as.character(conflicts$sample1),
      to   = as.character(conflicts$sample2)
   )
   conflict_graph = igraph::graph_from_data_frame(edges, directed = FALSE)
   membership = igraph::components(conflict_graph, mode = "weak")$membership

   conflicts$group = membership[as.character(conflicts$sample1)]
   conflict_groups = split(conflicts, by = "group", keep.by = FALSE)
   return(conflict_groups)
}

#' Boxplot of agreement scores for pairs of samples from the same vs. different donors
#'
#' @param ... Additional arguments passed to \code{\link[graphics]{boxplot}}.
#' @noRd
agreement_ = function(genotypes, donor_ids, sample_ids, ...) {

   checkmate::assert_list(genotypes)
   checkmate::assert_names(names(genotypes), must.include = c("snps", "outliers", "gamma"))

   J = ncol(genotypes$snps)
   
   checkmate::assert_atomic(donor_ids, len = J)
   checkmate::assert_atomic(sample_ids, len = J, unique = TRUE)

   pairs = pairwise_snp_agreement_(genotypes, donor_ids = donor_ids, sample_ids = sample_ids)

   # Plot a boxplot comparing agreement distribution for same vs. different donors
   boxplot(agreement ~ I(donor1 == donor2), data = pairs, horizontal = TRUE, ylim = c(0, 1), xlab = "Agreement", ...)
   
   return(invisible(NULL))
}

#' @rdname check_snp_agreement
#' @export
#' 
enumerate_sample_donors = function(genotypes) {

   checkmate::assert_list(genotypes)
   checkmate::assert_names(names(genotypes), must.include = c("snps", "outliers", "gamma"))

   AGREEMENT_THRESHOLD = 0.90
   J = ncol(genotypes$snps)

   pairs = pairwise_snp_agreement_(genotypes, sample_ids = seq_len(J))

   # Select sample pairs sharing the same genetic fingerprint (high agreement)
   matching_pairs = pairs[agreement > AGREEMENT_THRESHOLD]

   # If there are no highly agreeing pairs, every sample belongs to a unique donor
   if (nrow(matching_pairs) == 0) return(seq_len(J))

   # Connected components of the matching-pairs graph define the inferred donors
   edges = data.frame(from = matching_pairs$sample1, to = matching_pairs$sample2)
   donor_graph = igraph::graph_from_data_frame(
      d        = edges,
      vertices = data.frame(name = seq_len(J)),
      directed = FALSE
   )
   unname(igraph::components(donor_graph, mode = "strong")$membership)
}

#' Build a data.table of unique sample pairs with their pairwise SNP agreement scores
#' @noRd
pairwise_snp_agreement_ = function(genotypes, donor_ids = NULL, sample_ids = seq_len(ncol(genotypes$snps))) {
   if (is.null(donor_ids)) {
      pairs = CJ(sample_ids, sample_ids, sorted = FALSE)
      setnames(pairs, 1:2, c("sample1", "sample2"))
   } else {
      pairs = cbind(
         CJ(donor_ids,  donor_ids,  sorted = FALSE),
         CJ(sample_ids, sample_ids, sorted = FALSE)
      )
      setcolorder(pairs, c(1, 3, 2, 4))
      setnames(pairs, 1:4, c("donor1", "sample1", "donor2", "sample2"))
   }

   pairs$agreement = as.numeric(snp_agreement_matrix_(genotypes))
   pairs[!is.na(agreement)]
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
