#' @import data.table
NULL

#' Read .idat files
#' 
#' @description Import and annotate the data from .idat files.
#' 
#' @param idat_files Character vector of relative or absolute filepaths without
#' "_Grn.idat" and "_Red.idat" suffixes.
#' IDATs for red and green channel must have the same prefix and be stored in the same folder.
#'
#' E.g., a sample with the files "file/path/sample1_Red.idat" and 
#' "file/path/sample1_Grn.idat" would be passed as "file/path/sample1".
#'
#' All .idat files are assumed to be from the same platform.
#'
#' @param quiet If \code{TRUE}, suppresses the progress bar (useful for RMarkdown scripts).
#' @param on_disk If \code{TRUE}, data will be stored on disk using the \code{ff} package.
#'    This will be slower but enables processing datasets larger than memory.
#' 
#' @return A list containing
#' \item{platform}{Name of the BeadChip (450K/EPICv1/EPICv2/MOUSE)}
#' \item{manifest}{A data.table describing the probes}
#' \item{M}{intensities of targeting methylated sequences}
#' \item{U}{intensities of targeting unmethylated sequences}
#' \item{N}{number of beads from which average intensities in M were derived}
#' \item{V}{number of beads from which average intensities in U were derived}
#' \item{ctrlG}{Intensities of the control probes in the green color channel}
#' \item{ctrlR}{Intensities of the control probes in the red color channel}
#' \item{meta}{A data.table containing unique sample_ids and metadata}
#' 
#' @examples
#' \dontrun{
#' read_idats('9976861004_R01C01')
#' }
#'
#' @export
#'
read_idats = function(idat_files, quiet = FALSE, on_disk = FALSE)
{
   if (on_disk) {
      if(!requireNamespace("ff", quietly = TRUE))
         stop("Package 'ff' is required but not installed.")
   }

   ## illuminaio can handle gzipped .idats
   df_idats = tibble::tibble(basename = idat_files)

   zipped = !file.exists(paste0(idat_files, "_Grn.idat"))
   suffix = ifelse(zipped, ".idat.gz", ".idat")

   ex = file.exists(paste0(idat_files,"_Grn",suffix)) & file.exists(paste0(idat_files,"_Red",suffix))
   if (!all(ex)) stop("Some .idat files are missing")

   ## How many different features/bead types are there?
   # All idats are assumed to be from the same platform.
   # (Consider that Type I probes have each two different beads)
   P = illuminaio::readIDAT(paste0(idat_files[1], "_Grn", suffix[1]))$nSNPsRead
   print(P) # Print this in case someone encounters a new variant of a chip

   ## Pick appropriate manifest
   # (numbers are possible #features I have encountered in the wild so far)
   # We assume that all .idats use the same platform they can differ in #features, though

   if (P == 55300) {
      platform = "27K"
      chr      = "chr36"
      mapinfo  = "mapinfo36"
   } else if (P == 622399) { 
      platform = "450K"
      chr      = "chr37"
      mapinfo  = "mapinfo37"
   } else if (P %in% c(1051815, 1051943, 1052641)) {
      platform = "EPICv1"
      chr      = "chr38"
      mapinfo  = "mapinfo38"
   } else if (P == 1105209) {
      platform = "EPICv2"
      chr      = "chr38"
      mapinfo  = "mapinfo38"
   } else if (P == 361821) {
      platform = "MOUSE"
      chr      = "chr"
      mapinfo  = "mapinfo"
   } else {
      stop("Unknown platform")
   }

   manifest = ewastools:::MANIFESTS[[ platform ]]
   controls = ewastools:::CONTROLS [[ platform ]]

   setDT(manifest)
   setDT(controls)

   # Harmonize chromosome and mapping coordinate column names
   setnames(manifest, c(chr, mapinfo), c("chr", "mapinfo"))

   # (Consider moving this to "manifest.R")
   manifest[, index := 1L:.N]
   controls[, index := 1L:.N]

   manifest[channel == "Grn", OOBi := 1:.N]
   manifest[channel == "Red", OOBi := 1:.N]

   I = nrow(manifest)
   J = length(idat_files)

   # Allocate the matrices that will hold the data
   if (on_disk) {

      # Methylated (M) and unmethylated (U) signal intensities
      # ("double" because these intensities will be adjusted)
      M = ff::ff(vmode = "double",  dim = c(I, J))
      U = ff::ff(vmode = "double",  dim = c(I, J))
      # Standard deviations for the methylated (S) and unmethylated (T) intensities
      S = ff::ff(vmode = "integer", dim = c(I, J))
      T = ff::ff(vmode = "integer", dim = c(I, J))
      # Number of beads underlying methylated (N) and unmethylated (V) intensities
      N = ff::ff(vmode = "integer", dim = c(I, J))
      V = ff::ff(vmode = "integer", dim = c(I, J))

      oobG = list(
         M = ff::ff(vmode = "double", dim = c(manifest[channel == "Red", .N], J)),
         U = ff::ff(vmode = "double", dim = c(manifest[channel == "Red", .N], J)))

      oobR = list(
         M = ff::ff(vmode = "double", dim = c(manifest[channel == "Grn", .N], J)),
         U = ff::ff(vmode = "double", dim = c(manifest[channel == "Grn", .N], J)))

   } else {

      # Methylated (M) and unmethylated (U) signal intensities
      # ("double" because these intensities will be adjusted)
      M = U = matrix(NA_real_   ,nrow = I, ncol = J)
      # Standard deviations for the methylated (S) and unmethylated (T) intensities
      S = T = matrix(NA_integer_,nrow = I, ncol = J)
      # Number of beads underlying methylated (N) and unmethylated (V) signal intensities
      N = V = matrix(NA_integer_,nrow = I, ncol = J)

      oobG = list(
         M = matrix(NA_real_, nrow = manifest[channel == "Red", .N], ncol = J),
         U = matrix(NA_real_, nrow = manifest[channel == "Red", .N], ncol = J))

      oobR = list(
         M = matrix(NA_real_, nrow = manifest[channel == "Grn", .N], ncol = J),
         U = matrix(NA_real_, nrow = manifest[channel == "Grn", .N], ncol = J))
   }

   # Signal intensities of control probes
   # (as their number is low, these will be always stored in memory, even for `on_disk == TRUE`)
   ctrlG = ctrlR = matrix(NA_real_, nrow = nrow(controls), ncol = J)
   ctrlN = matrix(NA_integer_, nrow = nrow(controls), ncol = J)

   if (!quiet) pb = txtProgressBar(min = 0, max = J, style = 3)

   barcodes  = rep(NA_character_, J)
   positions = rep(NA_character_, J)
   dates     = rep(NA_character_, J)

   for (j in 1:J) {

      red = illuminaio::readIDAT(paste0(idat_files[j], "_Red", suffix[j]))
      grn = illuminaio::readIDAT(paste0(idat_files[j], "_Grn", suffix[j]))

      idat_order = red$MidBlock
      if(!identical(idat_order, grn$MidBlock)) stop("Red and green .idat files do not agree!")

      barcodes [j] = red$Barcode
      positions[j] = red$Unknowns$MostlyA

      # This information is sometimes not recorded
      if (nrow(red$RunInfo) >1 ) dates[j] = red$RunInfo[2,1]

      manifest[,Ui := match(addressU, idat_order)]
      manifest[,Mi := match(addressM, idat_order)]
      controls[, i := match(address , idat_order)]

      setindexv(manifest,"channel")
      manifest["Both", Mi := Ui, on = "channel"]

      ## Type I Red probes
      i = manifest["Red", on = "channel"]

      U[ i$index,j ] = red$Quants[ i$Ui, 1] # Mean
      T[ i$index,j ] = red$Quants[ i$Ui, 2] # SD
      V[ i$index,j ] = red$Quants[ i$Ui, 3] # NBeads
     
      M[ i$index,j ] = red$Quants[ i$Mi, 1] # Mean
      S[ i$index,j ] = red$Quants[ i$Mi, 2] # SD
      N[ i$index,j ] = red$Quants[ i$Mi, 3] # NBeads

      oobG$U[ i$OOBi,j ] = grn$Quants[ i$Ui,1 ]
      oobG$M[ i$OOBi,j ] = grn$Quants[ i$Mi,1 ]

      ## Type I Green probes
      i = manifest["Grn", on = "channel"]

      U[ i$index,j ] = grn$Quants[ i$Ui,1 ] # Mean
      T[ i$index,j ] = grn$Quants[ i$Ui,2 ] # SD
      V[ i$index,j ] = grn$Quants[ i$Ui,3 ] # NBeads

      M[ i$index,j ] = grn$Quants[ i$Mi,1 ] # Mean
      S[ i$index,j ] = grn$Quants[ i$Mi,2 ] # SD
      N[ i$index,j ] = grn$Quants[ i$Mi,3 ] # NBeads

      oobR$U[ i$OOBi,j ] = red$Quants[ i$Ui,1 ]
      oobR$M[ i$OOBi,j ] = red$Quants[ i$Mi,1 ]

      ## Type II probes (measured in both channels: Red for Unmethylated, Green for Methylated)
      i = manifest["Both", on = "channel"]

      U[ i$index,j ] = red$Quants[ i$Ui,1 ] # Mean
      T[ i$index,j ] = red$Quants[ i$Ui,2 ] # SD
      V[ i$index,j ] = red$Quants[ i$Ui,3 ] # NBeads

      M[ i$index,j ] = grn$Quants[ i$Mi,1 ] # Mean
      S[ i$index,j ] = grn$Quants[ i$Mi,2 ] # SD
      N[ i$index,j ] = grn$Quants[ i$Mi,3 ] # NBeads

      ## Control probes
      # Not keeping the SD for control probes ATM
      ctrlR[ controls$index,j ] = red$Quants[ controls$i,1 ]
      ctrlG[ controls$index,j ] = grn$Quants[ controls$i,1 ]
      ctrlN[ controls$index,j ] = red$Quants[ controls$i,3 ]

      if(!quiet) setTxtProgressBar(pb, j)

   }

   if(!quiet) close(pb)

   # Probes with zero beads on the chip (result of random chip assembly)
   # are set to intensity zero in the .idat files. Mark them as missing (NA).
   for (j in 1:J) {
      U[V[, j] == 0, j] = NA
      M[N[, j] == 0, j] = NA
      
      # Need at least two bead observations to calculate standard deviation
      T[V[, j] < 2, j] = NA
      S[N[, j] < 2, j] = NA
   }

   ctrlG[ctrlN == 0] = NA
   ctrlR[ctrlN == 0] = NA

   # Extract sample IDs from the file names (using the last portion of the path)
   sample_ids = strsplit(x = idat_files, split = "/")
   sample_ids = sapply(sample_ids, tail, n = 1L)

   meta = tibble::tibble(
       sample_id = sample_ids
      ,date = as.IDate(dates, "%m/%d/%Y %r")
      ,time = as.ITime(dates, "%m/%d/%Y %r")
      ,barcode = barcodes
      ,position = positions
      )

   raw = list(
       platform = platform
      ,manifest = manifest
      ,U = U, T = T, V = V
      ,M = M, S = S, N = N
      ,controls = controls
      ,ctrlG = ctrlG, ctrlR = ctrlR, ctrlN = ctrlN
      ,oobG = oobG, oobR = oobR
      ,meta = meta
   )

   return(raw)
}

#' Drop samples
#' @description Drop samples from the object returned by \code{read_idats()}.
#'    Used for removing samples that failed quality control before computing
#'    beta-values.
#'
#' @param raw Output of \code{\link{read_idats}}
#' @param j Indices of the samples to drop
#'
#' @return A modified \code{raw} object
#'
#' @export
#'
drop_samples = function(raw, j = NULL){

   if(is.null(j)) return(raw)


   if (ff::is.ff(raw$M)) {

      drop_cols = function(x, j) { # x must be an ff matrix
         cols_to_keep = setdiff(1:ncol(x), j)
         res = ff::ff(vmode = ff::vmode(x), dim = c(nrow(x), length(cols_to_keep)))
         if (length(cols_to_keep) > 0) {
            for (k in seq_along(cols_to_keep))
               res[, k] = x[, cols_to_keep[k]]
         }
         rownames(res) = rownames(x)
         colnames(res) = colnames(x)[cols_to_keep]
         return(res)
      }

      raw$U = drop_cols(raw$U, j)
      raw$M = drop_cols(raw$M, j)

      raw$T = drop_cols(raw$T, j)
      raw$S = drop_cols(raw$S, j)

      raw$V = drop_cols(raw$V, j)
      raw$N = drop_cols(raw$N, j)

      raw$ctrlG = drop_cols(raw$ctrlG, j)
      raw$ctrlR = drop_cols(raw$ctrlR, j)

      raw$oobG$M = drop_cols(raw$oobG$M, j)
      raw$oobG$U = drop_cols(raw$oobG$U, j)
      raw$oobR$M = drop_cols(raw$oobR$M, j)
      raw$oobR$U = drop_cols(raw$oobR$U, j)

      if("detP" %in% names(raw))
         raw$detP = drop_cols(raw$detP, j)

   } else {
      # Standard in-memory matrix column removal
      raw$U = raw$U[, -j, drop = FALSE]
      raw$M = raw$M[, -j, drop = FALSE]

      raw$T = raw$T[, -j, drop = FALSE]
      raw$S = raw$S[, -j, drop = FALSE]

      raw$V = raw$V[, -j, drop = FALSE]
      raw$N = raw$N[, -j, drop = FALSE]

      raw$ctrlG = raw$ctrlG[, -j, drop = FALSE]
      raw$ctrlR = raw$ctrlR[, -j, drop = FALSE]

      raw$oobG$M = raw$oobG$M[, -j, drop = FALSE]
      raw$oobG$U = raw$oobG$U[, -j, drop = FALSE]
      raw$oobR$M = raw$oobR$M[, -j, drop = FALSE]
      raw$oobR$U = raw$oobR$U[, -j, drop = FALSE]

      if("detP" %in% names(raw))
        raw$detP = raw$detP[, -j, drop = FALSE]
   }

   raw$meta = raw$meta[-j,, drop = FALSE]
   
   return(raw)
}
