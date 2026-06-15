#' @import data.table
NULL

#' Dye-bias correction and beta-values
#'
#' @description Preprocessing functions to correct dye-bias using the RELIC method and convert 
#' intensities to beta-values.
#'
#' @author Jonathan Heiss
#' @name Preprocessing
#' @rdname Preprocessing
#' @param raw Output of \code{\link{read_idats}}.
#'
#' @return \code{correct_dye_bias} and \code{correct_dye_bias2} return a modified \code{raw} list 
#'   object with dye-bias corrected intensities in \code{M}, \code{U}, and out-of-band signals.
#' @return \code{beta_values} returns a matrix of beta-values. If IDATs were imported with 
#'   \code{on_disk = TRUE}, this returns an \code{ff} object.
#' 
#' @references{Xu Z, Langie SA, De Boever P, Taylor JA, Niu L. RELIC: a novel dye-bias correction method for Illumina Methylation BeadChip. BMC Genomics. 2017 Jan 3;18(1):4.}
#' @references{Heiss JA, Brenner H. Between-array normalization for 450K data. Frontiers in Genetics. 2015;6.}
#' 
NULL

#' @rdname Preprocessing
#' @export
#' 
correct_dye_bias <- function(raw)
{
   stopifnot(c("manifest", "M", "U", "controls", "ctrlG", "ctrlR") %in% names(raw))

   i1g = raw$manifest[channel == "Grn" ,]
   i2  = raw$manifest[channel == "Both",]

   ## Normalization control probes.
   # These normalization controls are paired: G is paird with A, C is paired with T.
   # TODO: Check that the indices of the paired probes match.
   # Green channel intensities of NORM_G/NORM_C are regressed on Red channel intensities of NORM_A/NORM_T.
   Ai = raw$controls[group == "NORM_A"][order(name)]$index
   Ti = raw$controls[group == "NORM_T"][order(name)]$index
   Ci = raw$controls[group == "NORM_C"][order(name)]$index
   Gi = raw$controls[group == "NORM_G"][order(name)]$index

   J = ncol(raw$M)

   for(j in 1:J){
      # Regress intensities in the red color channel on those in the green color channel.
      # The relation is assumed to be linear on the log scale.
      x = log(raw$ctrlG[c(Gi, Ci), j])
      y = log(raw$ctrlR[c(Ai, Ti), j])

      # Filter out missing or infinite values
      keep = !is.na(y) & !is.na(x) & is.finite(x) & is.finite(y)
      x = x[keep]; y = y[keep]

      # Theil-Sen robust linear regression to estimate the slope (intercept is ignored)
      m = mblm::mblm(y~x, repeated = FALSE)

      # Adjust the intensities in the Green channel using the regression coefficients
      i = i2$index
      raw$M[i,j] = exp(stats::coef(m)[1] + log(raw$M[i,j]) * stats::coef(m)[2])

      i = i1g$index
      raw$M[i,j] = exp(stats::coef(m)[1] + log(raw$M[i,j]) * stats::coef(m)[2])
      raw$U[i,j] = exp(stats::coef(m)[1] + log(raw$U[i,j]) * stats::coef(m)[2])

      # Also adjust the out-of-band green signals
      raw$oobG$U[,j] = exp(stats::coef(m)[1] + log(raw$oobG$U[,j]) * stats::coef(m)[2])
      raw$oobG$M[,j] = exp(stats::coef(m)[1] + log(raw$oobG$M[,j]) * stats::coef(m)[2])
   }

   return(raw)
}

#' @rdname Preprocessing
#' 
correct_dye_bias2 = function (raw) 
{
   # Experimental version of dye-bias correction. I found that red~green differs for 
   # ... G~A and C~T. They are however not independent. I have not understood this relation yet.
   stopifnot(c("manifest", "M", "U", "controls", "ctrlG", "ctrlR") %in% names(raw))

   i1g = raw$manifest[channel == "Grn", ]
   i2  = raw$manifest[channel == "Both",]

   Ai = raw$controls[group == "NORM_A"][order(name)]$index; Ai = raw$ctrlR[Ai,]
   Gi = raw$controls[group == "NORM_G"][order(name)]$index; Gi = raw$ctrlG[Gi,]
   Ti = raw$controls[group == "NORM_T"][order(name)]$index; Ti = raw$ctrlR[Ti,]
   Ci = raw$controls[group == "NORM_C"][order(name)]$index; Ci = raw$ctrlG[Ci,]

   J = ncol(raw$M)

   for (j in 1:J) {
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

#' Convert fluorescence intensities to beta-values
#' @rdname Preprocessing
#' @return A matrix of beta-values.
#'   If IDATs were imported with option \code{on_disk = TRUE}, a \code{ff} object.
#'   If the data was imported using \code{on_disk = TRUE}, it is recommended to transpose
#'   the returned ff object for fast row-based access.
#' @export
#'
beta_values = function(raw)
{
   stopifnot(c("manifest", "M", "U", "meta") %in% names(raw))

   I = nrow(raw$M)
   J = ncol(raw$M)

   if (ff::is.ff(raw$M)) {
      meth = ff::ff(NA_real_, dim = c(I, J))
   } else {
      meth = matrix(NA_real_, nrow = I, ncol = J)
   }

   # Calculate beta = M / (M + U).
   # Floor intensity values at 1 to prevent division by zero
   for (j in 1:J) {
      Uj = raw$U[, j]; Uj[Uj < 1] = 1
      Mj = raw$M[, j]; Mj[Mj < 1] = 1
      meth[, j] = Mj / (Uj + Mj)
   }
     
   rownames(meth) = paste0(raw$manifest$ilmn_id)
   colnames(meth) = raw$meta$sample_id

   return(meth)
}

#' Match probe IDs between a query list and matrix row names
#' 
#' @description This function parses legacy and EPICv2 probe ID formats to find matching
#' indices in the dataset. It supports mapping legacy probe IDs to EPICv2 data.
#' 
#' @param query Character vector of probe IDs to look up.
#' @param rownames Character vector of row names (probe IDs) of the target matrix.
#' 
#' @return An integer vector of the matching indices of \code{query} in \code{rownames}.
#'
find_matching_rows = function(query, rownames)
{
   # Regex for legacy Illumina probe formats (cg/rs/ch followed by numbers)
   legacy_regex = stringr::regex("^cg\\d{8,8}$|^rs\\d+$|^ch\\.\\w+\\.\\d+[FR]$")
   # Regex for EPICv2 probe formats (containing strand and replicate tags like _BC11)
   epicv2_regex = stringr::regex("^[cg|ch|rs|nv].*_[TB][CO]\\d+$")

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
      # Find for each user-provided `ilmn_id` the `probe_id` with `probe_rep` == 1
      df_map = data.table(ewastools:::MANIFESTS$EPICv2[,c("ilmn_id", "probe_id")])
      setkey(df_map, probe_id)
      query = df_map[query, nomatch = NA, mult = "first"]$ilmn_id
      return(match(query, rownames))
   } else if (query_type == "epicv2" & row_type == "legacy") {
      stop("Query contains EPICv2 probe IDs but dataset is of legacy type")
   } else {
      stop("This should be a dead branch")
   }
}
