#' Various methylation scores
#'
#' @description Calculate various epigenetic scores (e.g. maternal smoking, gestational age, Horvath clock)
#' based on beta-values of DNA methylation.
#'
#' @rdname MethylationScores
#' @param beta Matrix of beta-values with probe IDs as row names and samples as columns.
#' @param model Name of the model. Currently supported models are:
#' \itemize{
#'   \item \code{"cord_blood:maternal_smoking"}: Score for in utero exposure to maternal smoking based on cord blood [1].
#'   \item \code{"cord_blood:gestational_age"}: Gestational age (in weeks) based on cord blood [2].
#'   \item \code{"placenta:gestational_age"}: Gestational age (in weeks) based on placenta [3].
#'   \item \code{"horvath_clock"}: Multi-tissue epigenetic age (in years) [4].
#' }
#' @return A numeric vector of calculated methylation scores, one for each sample (column) of \code{beta}.
#' Returns \code{NULL} if \code{model} is not one of the supported models.
#'
#' @references{[1] Reese SE, et al. DNA methylation score as a biomarker in newborns for sustained maternal smoking during pregnancy. Environmental Health Perspectives. 2017;125(4):760-6. DOI:10.1289/EHP333}
#' @references{[2] Knight AK, et al. An epigenetic clock for gestational age at birth based on blood methylation data. Genome Biology. 2016;17:206. DOI:10.1186/s13059-016-1068-z}
#' @references{[3] Lee Y, et al. Placental epigenetic clocks: estimating gestational age using placental DNA methylation levels. Aging. 2019;11(12):4238-53. DOI:10.18632/aging.102049}
#' @references{[4] Horvath S. DNA methylation age of human tissues and cell types. Genome Biology. 2013;14:R115. DOI:10.1186/gb-2013-14-10-r115}
#' @export
#'
methylation_score = function(beta, model) {

   if (model == "cord_blood:maternal_smoking") {

      # Model based on Table S5 from Reese et al. (2017) available at doi.org/10.1289/EHP333
      # Predicts newborn maternal smoking exposure from cord blood DNA methylation.
      markers = data.table(
          probe_id = c("cg02256631","cg02482603","cg04103532","cg04180046","cg04506190","cg05549655","cg05575921","cg08698721","cg09743950","cg10799846","cg12186702","cg13834112","cg13893782","cg14179389","cg14351425","cg14633298","cg14743346","cg17397069","cg19381766","cg22154659","cg22802102","cg23304605","cg25189904","cg25949550","cg26764244","cg27291468")
         ,coef = c(-0.191,2.706,1.786,14.027,2.318,6.210,-10.909,1.142,-6.330,-4.963,3.847,1.514,-0.963,-6.304,6.361,5.050,2.286,-2.912,5.245,-0.773,-0.254,0.011,-3.903,-46.991,-0.246,0.836)
         )
      intercept = 0

   } else if (model == "cord_blood:gestational_age"){

      # Knight, A.K., Craig, J.M., Theda, C. et al.
      # An epigenetic clock for gestational age at birth based on blood methylation data.
      # Genome Biol 17, 206 (2016).
      # Reference: doi.org/10.1186/s13059-016-1068-z
      markers = system.file("data/13059_2016_1068_MOESM3_ESM.csv", package="ewastools")
      markers = fread(markers)
      setnames(markers, 1:2, c("probe_id","coef"))
      markers = markers[-1]
      intercept = 41.72579759

   } else if (model == "placenta:gestational_age"){

      # Placental epigenetic clocks: estimating gestational age using placental DNA methylation levels
      # Reference: https://doi.org/10.18632/aging.102049
      markers = system.file("data/102049-SupFile1.csv", package="ewastools")
      markers = fread(markers)
      setnames(markers, 1:2, c("probe_id", "coef"))
      markers = markers[-1]
      intercept = 24.99772133

   } else if (model == "horvath_clock"){

      # Multi-tissue epigenetic age clock from Horvath (2013).
      # Reference: https://doi.org/10.1186/gb-2013-14-10-r115
      markers = system.file("data/horvath_clock.csv", package="ewastools")
      markers = fread(markers, header=TRUE)
      intercept = 0.695507258

   } else {
      return(NULL)
   }

   i = find_matching_rows(markers$probe_id, rownames(beta))
   if (any(is.na(i))) {
      warning(glue::glue("{sum(is.na(i))} of {length(i)} model probes were not found in the beta matrix."))
      # Drop missing probes before indexing `beta` (`ff` matrices do not allow `NA` subscripts)
      markers = markers[!is.na(i)]
      i = i[!is.na(i)]
   }

   # If none of the model probes are present in `beta`, return `NA` for every sample:
   # `ff` matrices cannot be indexed with `integer(0)`, and `colSums(..., na.rm = TRUE)`
   # on a 0-row matrix would otherwise evaluate to `0` and return the bare `intercept`.
   if (length(i) == 0L) {
      y = rep(NA_real_, ncol(beta))
      names(y) = colnames(beta)
      return(y)
   }

   sub_beta = beta[i, , drop = FALSE]
   all_missing = apply(is.na(sub_beta), 2, all)

   y = sub_beta * markers$coef
   y = colSums(y, na.rm = TRUE) + intercept

   # Samples with only `NA` values across matched probes should receive `NA`, not `intercept`
   y[all_missing] = NA_real_

   if (model == "horvath_clock") {
      # Epigenetic age predictor uses a transformation function to handle young ages.
      # Below is the inverse of the transformation function applied to linear predictions.
      inverse_transformation = function(x, adult.age = 20) {
         ifelse(x<0, (1+adult.age)*exp(x)-1, (1+adult.age)*x+adult.age) }

      y = inverse_transformation(y)
   }

   return(y)
}
