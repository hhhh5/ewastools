#' Dye-bias correction and beta-values
#'
#' @description \code{correct_dye_bias()} corrects the dye bias between the red and green color
#' channel. It implements an improved version of RELIC using a robust Theil-Sen estimator
#' to regress the intensities of the red on those of the green normalization control probes.
#'
#' \code{beta_values()} converts the methylated (M) and unmethylated (U) intensities into
#' beta-values (M / (M + U)). Intensities below 1 are set to 1 before the conversion.
#'
#' @author Jonathan Heiss
#' @name Preprocessing
#' @rdname Preprocessing
#' @param raw Output of \code{\link{read_idats}}.
#'
#' @return \code{correct_dye_bias} returns a modified \code{raw} list object with dye-bias corrected
#'   intensities in \code{M}, \code{U}, and the out-of-band signals.
#' @return \code{beta_values} returns a matrix of beta-values with probe IDs as row names and
#'   sample IDs as column names. If the data was imported with \code{on_disk = TRUE}, an \code{ff}
#'   object is returned instead. In this case, it is recommended to transpose the returned ff object
#'   for fast row-based access.
#'
#' @references{Xu Z, Langie SA, De Boever P, Taylor JA, Niu L. RELIC: a novel dye-bias correction method for Illumina Methylation BeadChip. BMC Genomics. 2017 Jan 3;18(1):4.}
#'
NULL

#' @rdname Preprocessing
#' @export
#' 
correct_dye_bias <- function(raw)
{
   checkmate::assert_list(raw)
   checkmate::assert_names(names(raw), must.include = c("manifest", "M", "U", "controls", "ctrlG", "ctrlR"))

   i1g = raw$manifest[channel == "Grn" ,]
   i2  = raw$manifest[channel == "Both",]

   ## Normalization control probes.
   # These normalization controls are paired: G is paired with A, C is paired with T.
   # TODO: Check that the indices of the paired probes match.
   # Green channel intensities of NORM_G/NORM_C are corrected
   #  using the Red channel intensities of NORM_A/NORM_T as reference.
   Ai = raw$controls[group == "NORM_A"][order(name)]$index
   Ti = raw$controls[group == "NORM_T"][order(name)]$index
   Ci = raw$controls[group == "NORM_C"][order(name)]$index
   Gi = raw$controls[group == "NORM_G"][order(name)]$index

   J = ncol(raw$M)

   for(j in seq_len(J)){
      # Regress intensities in the red color channel on those in the green color channel.
      # In the next step, the green intensities are plugged into the learned function
      #  and adjusted towards the red intensities.
      # The relation is assumed to be linear on the log scale.
      x = log(raw$ctrlG[c(Gi, Ci), j])
      y = log(raw$ctrlR[c(Ai, Ti), j])

      # Filter out missing or infinite values
      keep = !is.na(y) & !is.na(x) & is.finite(x) & is.finite(y)
      x = x[keep]; y = y[keep]

      # Theil-Sen robust linear regression on log-transformed intensities
      m = mblm::mblm(y~x, repeated = FALSE)

      rescale = function(intensities) exp(stats::coef(m)[1] + log(intensities) * stats::coef(m)[2])

      # Adjust only green-channel intensities toward the red channel
      i = i2$index
      raw$M[i,j] = rescale(raw$M[i,j])

      i = i1g$index
      raw$M[i,j] = rescale(raw$M[i,j])
      raw$U[i,j] = rescale(raw$U[i,j])

      # Also adjust the out-of-band green signals (from Type I Red probes)
      raw$oobG$U[,j] = rescale(raw$oobG$U[,j])
      raw$oobG$M[,j] = rescale(raw$oobG$M[,j])
   }

   return(raw)
}

#' Experimental dye-bias correction using separate regressions for A~G and T~C
#'
#' @param raw Output of \code{read_idats()}.
#' @noRd
correct_dye_bias2 = function (raw) 
{
   # Experimental version of dye-bias correction. I found that red~green differs for 
   # ... G~A and C~T. They are however not independent. I have not understood this relation yet.
   checkmate::assert_list(raw)
   checkmate::assert_names(names(raw), must.include = c("manifest", "M", "U", "controls", "ctrlG", "ctrlR"))

   i1g = raw$manifest[channel == "Grn", ]
   i2  = raw$manifest[channel == "Both",]

   Ai = raw$controls[group == "NORM_A"][order(name)]$index; Ai = raw$ctrlR[Ai, , drop = FALSE]
   Gi = raw$controls[group == "NORM_G"][order(name)]$index; Gi = raw$ctrlG[Gi, , drop = FALSE]
   Ti = raw$controls[group == "NORM_T"][order(name)]$index; Ti = raw$ctrlR[Ti, , drop = FALSE]
   Ci = raw$controls[group == "NORM_C"][order(name)]$index; Ci = raw$ctrlG[Ci, , drop = FALSE]

   J = ncol(raw$M)

   for (j in seq_len(J)) {
      # Regress A on G
      x = log(Gi[, j])
      y = log(Ai[, j])
      keep = !is.na(y) & !is.na(x) & is.finite(x) & is.finite(y)
      x = x[keep]
      y = y[keep]
      m = mblm::mblm(y ~ x, repeated = FALSE)

      i = i2$index
      raw$M[i, j] = exp(stats::coef(m)[1] + log(raw$M[i, j]) * stats::coef(m)[2])

      # Regress T on C
      x = log(Ci[, j])
      y = log(Ti[, j])
      keep = !is.na(y) & !is.na(x) & is.finite(x) & is.finite(y)
      x = x[keep]
      y = y[keep]
      m = mblm::mblm(y ~ x, repeated = FALSE)

      i = i1g$index
      raw$U[i, j] = exp(stats::coef(m)[1] + log(raw$U[i, j]) * stats::coef(m)[2])
      raw$M[i, j] = exp(stats::coef(m)[1] + log(raw$M[i, j]) * stats::coef(m)[2])
   }

   return(raw)
}

#' @rdname Preprocessing
#' @export
#'
beta_values = function(raw)
{
   checkmate::assert_list(raw)
   checkmate::assert_names(names(raw), must.include = c("manifest", "M", "U", "meta"))

   I = nrow(raw$M)
   J = ncol(raw$M)

   if (inherits(raw$M, "ff")) {
      meth = ff::ff(NA_real_, dim = c(I, J))
   } else {
      meth = matrix(NA_real_, nrow = I, ncol = J)
   }

   # Calculate beta = M / (M + U).
   # Floor intensity values at 1 to prevent division by zero
   for (j in seq_len(J)) {
      Uj = raw$U[, j]; Uj[Uj < 1] = 1
      Mj = raw$M[, j]; Mj[Mj < 1] = 1
      meth[, j] = Mj / (Uj + Mj)
   }
     
   rownames(meth) = as.character(raw$manifest$ilmn_id)
   colnames(meth) = raw$meta$sample_id

   return(meth)
}

#' Match probe IDs between a query list and matrix row names
#' 
#' @description This function parses legacy and EPICv2-style probe ID formats to find matching
#' indices in the dataset. It supports mapping legacy probe IDs to EPICv2/MSA/MOUSE data.
#' 
#' @param query Character vector of probe IDs to look up.
#' @param rownames Character vector of row names (probe IDs) of the target matrix.
#' 
#' @return An integer vector of the matching indices of \code{query} in \code{rownames}.
#' @noRd
find_matching_rows = function(query, rownames)
{
   # Regex for legacy Illumina probe formats (cg/rs/ch followed by numbers)
   legacy_regex = stringr::regex("^cg\\d{8,8}$|^rs\\d+$|^ch\\.\\w+\\.\\d+[FR]$")
   # Regex for EPICv2/MSA/MOUSE probe formats (containing strand and replicate tags like _BC11)
   epicv2_regex = stringr::regex("^(?:cg|ch|rs|nv|mu|rp).*_[TB][CO]\\d+$")

   if (all(stringr::str_detect(query, pattern = legacy_regex))) {
      query_type = "legacy"
   } else if (all(stringr::str_detect(query, pattern = epicv2_regex))) {
      query_type = "epicv2"
   } else {
      stop("Query is (partly) invalid!")
   }

   if (all(stringr::str_detect(rownames, pattern = legacy_regex))) {
      row_type = "legacy"
   } else if (all(stringr::str_detect(rownames, pattern = epicv2_regex))) {
      row_type = "epicv2"
   } else {
      stop("(Some) row names are invalid!")
   }

   if (query_type == row_type) {
      # Directly match if formats are identical
      return(match(query, rownames))
   } else if (query_type == "legacy" & row_type == "epicv2") {
      # Strip strand/design/replicate suffix (e.g. `_BC11`) from row names to match legacy `probe_id`
      # (`match()` returns the first matching replicate in manifest order)
      base_rownames = sub("_[TB][CO]\\d+$", "", rownames)
      return(match(query, base_rownames))
   } else if (query_type == "epicv2" & row_type == "legacy") {
      stop("Query contains EPICv2 probe IDs but dataset is of legacy type")
   } else {
      stop("This should be a dead branch")
   }
}
