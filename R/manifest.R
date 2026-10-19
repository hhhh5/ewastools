.manifest_cache = new.env(parent = emptyenv())

#' Load Illumina methylation array manifest and control probes
#'
#' Lazy-loads the manifest and controls for a specified platform from \code{inst/extdata/}
#' and caches it in memory.
#'
#' @param platform Platform name: "27K", "450K", "EPICv1", "EPICv2", "MSA", or "MOUSE"
#' @return A list with elements \code{manifest} and \code{controls}.
#' @noRd
read_manifest = function(platform) {
   if (!exists(platform, envir = .manifest_cache, inherits = FALSE)) {
      rds_file = system.file("extdata", paste0(platform, ".rds"), package = "ewastools")
      if (!nzchar(rds_file) || !file.exists(rds_file))
         stop(sprintf("Manifest for platform '%s' not found.", platform), call. = FALSE)
      assign(platform, readRDS(rds_file), envir = .manifest_cache)
   }
   get(platform, envir = .manifest_cache, inherits = FALSE)
}
