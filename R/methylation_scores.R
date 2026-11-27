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

   checkmate::assert_string(model)

   spec = load_score_model(model)
   if (is.null(spec)) return(NULL)

   markers   = spec$markers
   intercept = spec$intercept

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

   y = spec$transform(y)
   return(y)
}

#' Load probe weights, intercept, and post-transformation for a methylation score model
#' @noRd
load_score_model = function(model) {
   switch(
      model,

      "cord_blood:maternal_smoking" = {
         # Model based on Table S5 from Reese et al. (2017) available at doi.org/10.1289/EHP333
         # Predicts newborn maternal smoking exposure from cord blood DNA methylation.
         csv_path = system.file("data/reese_maternal_smoking.csv", package = "ewastools")
         markers = fread(csv_path, header = TRUE)
         list(markers = markers, intercept = 0, transform = identity)
      },

      "cord_blood:gestational_age" = {
         # Knight, A.K., Craig, J.M., Theda, C. et al.
         # An epigenetic clock for gestational age at birth based on blood methylation data.
         # Genome Biol 17, 206 (2016). Reference: doi.org/10.1186/s13059-016-1068-z
         csv_path = system.file("data/13059_2016_1068_MOESM3_ESM.csv", package = "ewastools")
         markers = fread(csv_path)
         setnames(markers, 1:2, c("probe_id", "coef"))
         markers = markers[-1] # Row 1 holds the intercept
         list(markers = markers, intercept = 41.72579759, transform = identity)
      },

      "placenta:gestational_age" = {
         # Placental epigenetic clocks: estimating gestational age using placental DNA methylation levels
         # Reference: https://doi.org/10.18632/aging.102049
         csv_path = system.file("data/102049-SupFile1.csv", package = "ewastools")
         markers = fread(csv_path)
         setnames(markers, 1:2, c("probe_id", "coef"))
         markers = markers[-1] # Row 1 holds the intercept
         list(markers = markers, intercept = 24.99772133, transform = identity)
      },

      "horvath_clock" = {
         # Multi-tissue epigenetic age clock from Horvath (2013).
         # Reference: https://doi.org/10.1186/gb-2013-14-10-r115
         csv_path = system.file("data/horvath_clock.csv", package = "ewastools")
         markers = fread(csv_path, header = TRUE)

         # Inverse of Horvath's age transformation (handles young ages < adult.age)
         inverse_transformation = function(x, adult.age = 20) {
            ifelse(x < 0, (1 + adult.age) * exp(x) - 1, (1 + adult.age) * x + adult.age)
         }

         list(markers = markers, intercept = 0.695507258, transform = inverse_transformation)
      },

      NULL
   )
}
