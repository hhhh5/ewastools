#' Quality control metrics
#'
#' \code{control_metrics()} calculates the metrics described in the
#' 'BeadArray Controls Reporter Software Guide' from Illumina.
#' \code{sample_failure()} compares these metrics against the thresholds recommended by Illumina.
#'
#' @export
#' @note You can download the Software Guide at
#' \url{https://support.illumina.com/downloads/beadarray-controls-reporter-software-guide-1000000004009.html}
#'
#' @param raw Output of \code{\link{read_idats}}.
#' @param FFPE Set to \code{TRUE} if the FFPE DNA restoration kit was used.
#' @param metrics Output of \code{control_metrics()}.
#'
#' @return For \code{control_metrics}, a list of 21 control metrics, each one a numeric vector equal in length to the number of samples.
#' Each vector has an attribute \code{"threshold"} holding the cut-off recommended by Illumina (samples below the threshold fail).
#' \item{Restoration}{Success of the FFPE DNA restoration. Green signal of the restoration control relative to background.
#'   Only evaluated if \code{FFPE = TRUE}.}
#' \item{Staining Green}{Biotin staining efficiency. Ratio of Biotin (High) to background green signals.}
#' \item{Staining Red}{DNP staining efficiency. Ratio of DNP (High) to background red signals.}
#' \item{Extension Green}{Extension efficiency of C/G base pairs. Ratio of green intensities.}
#' \item{Extension Red}{Extension efficiency of A/T base pairs. Ratio of red intensities.}
#' \item{Hybridization High/Medium}{Hybridization performance. Ratio of high-to-medium green signal intensities.}
#' \item{Hybridization Medium/Low}{Hybridization performance. Ratio of medium-to-low green signal intensities.}
#' \item{Target Removal 1}{Target removal efficiency 1. Ratio of green intensities (+3000 background offset).}
#' \item{Target Removal 2}{Target removal efficiency 2. Ratio of green intensities (+3000 background offset).}
#' \item{Bisulfite Conversion I Green}{Bisulfite conversion efficiency for Type I green probes. Ratio of green intensities.}
#' \item{Bisulfite Conversion I Green (Bkg)}{Bisulfite conversion green background ratio (+3000 background offset).}
#' \item{Bisulfite Conversion I Red}{Bisulfite conversion efficiency for Type I red probes. Ratio of red intensities.}
#' \item{Bisulfite Conversion I Red (Bkg)}{Bisulfite conversion red background ratio (+3000 background offset).}
#' \item{Bisulfite Conversion II}{Bisulfite conversion efficiency for Type II probes. Ratio of red-to-green intensities.}
#' \item{Bisulfite Conversion II (Bkg)}{Bisulfite conversion Type II green background ratio (+3000 background offset).}
#' \item{Specificity I Green}{Specificity of hybridization for Type I green probes. Ratio of green intensities.}
#' \item{Specificity I Red}{Specificity of hybridization for Type I red probes. Ratio of red intensities.}
#' \item{Specificity II}{Specificity of hybridization for Type II probes. Ratio of red-to-green intensities.}
#' \item{Specificity II (Bkg)}{Specificity of hybridization green background ratio (+3000 background offset).}
#' \item{Non-polymorphic Green}{Non-polymorphic query performance in the green channel.}
#' \item{Non-polymorphic Red}{Non-polymorphic query performance in the red channel.}
#' 
#' @return For \code{sample_failure}, a logical vector, \code{TRUE} if the sample at the corresponding index falls below
#' the threshold for any of the control metrics.
#'
control_metrics = function(raw, FFPE = FALSE){

   checkmate::assert_list(raw)
   checkmate::assert_names(names(raw), must.include = c("platform", "controls", "ctrlG", "ctrlR"))
   checkmate::assert_choice(raw$platform, unique(CONTROL_METRICS$platform))
   checkmate::assert_flag(FFPE)

   controls = raw$controls

   # Build channel lookup helpers returning the submatrix of matching control probes
   select_probes = function(mat) {
      function(group, pattern) {
         keep = controls$group == group & grepl(pattern, controls$name)
         if (!any(keep)) {
            stop(sprintf("No control probes matched group '%s' and pattern '%s' for platform '%s'.",
                         group, pattern, raw$platform), call. = FALSE)
         }
         mat[controls$index[keep], , drop = FALSE]
      }
   }

   G = select_probes(raw$ctrlG)
   R = select_probes(raw$ctrlR)

   # Extension probes incorporate C/G (green) or A/T (red); their opposite-channel intensities
   # (A/T in green, C/G in red) plus a +3000 offset serve as channel-specific background references
   # per Illumina's BeadArray Controls Reporter (BACR) guide.
   bkgG = cmax(G("EXTENSION", "^Extension \\([AT]\\)$")) + 3000
   bkgR = cmax(R("EXTENSION", "^Extension \\([CG]\\)$")) + 3000

   specs = CONTROL_METRICS[CONTROL_METRICS$platform == raw$platform, ]
   # Threshold 0 disables failure flagging in `sample_failure()` for non-FFPE samples
   if (FFPE) specs$threshold[specs$name == "Restoration"] = 1

   J = ncol(raw$ctrlG)
   metrics = vector("list", nrow(specs))
   names(metrics) = specs$name

   for (k in seq_len(nrow(specs))) {
      val = as.numeric(specs$fn[[k]](G, R, bkgG, bkgR))
      if (length(val) != J) {
         stop(sprintf("Control metric '%s' produced %d values for %d samples.",
                      specs$name[k], length(val), J), call. = FALSE)
      }
      attr(val, "threshold") = specs$threshold[k]
      metrics[[k]] = val
   }

   return(metrics)
}

#' @rdname control_metrics
#' @export
#' 
sample_failure = function(metrics){
   checkmate::assert_list(metrics, min.len = 1)
   # Returns TRUE if any metric falls below its threshold (ignoring NA metrics so missing probes
   # neither mask failures in other metrics nor trigger false positives)
   Reduce(`|`, lapply(metrics, \(m) !is.na(m) & m < attr(m, "threshold")))
}

cmin = function(x) if (nrow(x) == 1L) x[1, ] else apply(x, 2, min)
cmax = function(x) if (nrow(x) == 1L) x[1, ] else apply(x, 2, max)

# Specification of Illumina's BeadArray Controls Reporter (BACR) QC metrics per platform.
# Every multi-probe metric aggregates across probes with cmin()/cmax() (order-invariant), and
# every single-probe metric matches exactly one probe (enforced via length check in control_metrics()).
#
# Assay biology notes:
# - Restoration: FFPE DNA Restore Kit efficiency (threshold raised to 1 when FFPE = TRUE).
# - Bisulfite conversion: queries [C/T] polymorphisms created by bisulfite conversion of non-CpG C.
#   (For BS Conversion I Green Bkg, the manual omits highest/lowest; `cmax` mirrors the red channel.)
# - Hybridization: synthetic targets at High/Medium/Low (5/1/0.2 pM) in IBX buffer (green channel).
# - Extension / Staining: sample-independent controls (hairpin extension; DNP/biotin beads).
# - Target removal: stripping efficiency after extension (low signal indicates efficient removal).
CONTROL_METRICS = tibble::tribble(
   ~platform, ~name,                                ~threshold, ~fn,
   # 27K (15 metrics; no Restoration, BS Conversion II, Target Removal 2, or Specificity II)
   "27K",     "Bisulfite Conversion I Green",       1,          \(G, R, bkgG, bkgR) cmin(G("BISULFITE CONVERSION I", "^BS Conversion I-C[123]$")) / cmax(G("BISULFITE CONVERSION I", "^BS Conversion I-U[123]$")),
   "27K",     "Bisulfite Conversion I Green (Bkg)", 1,          \(G, R, bkgG, bkgR) bkgG / cmax(G("BISULFITE CONVERSION I", "^BS Conversion I-U[123]$")),
   "27K",     "Bisulfite Conversion I Red",         1,          \(G, R, bkgG, bkgR) cmin(R("BISULFITE CONVERSION I", "^BS Conversion I-C[456]$")) / cmax(R("BISULFITE CONVERSION I", "^BS Conversion I-U[456]$")),
   "27K",     "Bisulfite Conversion I Red (Bkg)",   1,          \(G, R, bkgG, bkgR) bkgR / cmax(R("BISULFITE CONVERSION I", "^BS Conversion I-U[456]$")),
   "27K",     "Hybridization High/Medium",          1,          \(G, R, bkgG, bkgR) G("HYBRIDIZATION", "^Hyb \\(High\\)$") / G("HYBRIDIZATION", "^Hyb \\(Medium\\)$"),
   "27K",     "Hybridization Medium/Low",           1,          \(G, R, bkgG, bkgR) G("HYBRIDIZATION", "^Hyb \\(Medium\\)$") / G("HYBRIDIZATION", "^Hyb \\(Low\\)$"),
   "27K",     "Extension Green",                    5,          \(G, R, bkgG, bkgR) cmin(G("EXTENSION", "^Extension \\([CG]\\)$")) / cmax(G("EXTENSION", "^Extension \\([AT]\\)$")),
   "27K",     "Extension Red",                      5,          \(G, R, bkgG, bkgR) cmin(R("EXTENSION", "^Extension \\([AT]\\)$")) / cmax(R("EXTENSION", "^Extension \\([CG]\\)$")),
   "27K",     "Staining Green",                     5,          \(G, R, bkgG, bkgR) G("STAINING", "^Biotin \\(High\\)$") / G("STAINING", "^Biotin \\(Bkg\\)$"),
   "27K",     "Staining Red",                       5,          \(G, R, bkgG, bkgR) R("STAINING", "^DNP \\(High\\)$") / R("STAINING", "^DNP \\(Bkg\\)$"),
   "27K",     "Target Removal 1",                   1,          \(G, R, bkgG, bkgR) bkgG / G("TARGET REMOVAL", "^Target Removal 1$"),
   "27K",     "Specificity I Green",                1,          \(G, R, bkgG, bkgR) cmin(G("SPECIFICITY I", "^GT Mismatch [123] \\(PM\\)$")) / cmax(G("SPECIFICITY I", "^GT Mismatch [123] \\(MM\\)$")),
   "27K",     "Specificity I Red",                  1,          \(G, R, bkgG, bkgR) cmin(R("SPECIFICITY I", "^GT Mismatch 4 \\(PM\\)$")) / cmax(R("SPECIFICITY I", "^GT Mismatch 4 \\(MM\\)$")),
   "27K",     "Non-polymorphic Green",              5,          \(G, R, bkgG, bkgR) cmin(G("NON-POLYMORPHIC", "^NP \\([CG]\\)$")) / cmax(G("NON-POLYMORPHIC", "^NP \\([AT]\\)$")),
   "27K",     "Non-polymorphic Red",                5,          \(G, R, bkgG, bkgR) cmin(R("NON-POLYMORPHIC", "^NP \\([AT]\\)$")) / cmax(R("NON-POLYMORPHIC", "^NP \\([CG]\\)$")),

   # 450K (21 metrics; BS Conversion I Green = C1-C3/U1-U3, Red = C4-C6/U4-U6)
   "450K",    "Restoration",                        0,          \(G, R, bkgG, bkgR) G("RESTORATION", "^Restore$") / bkgG,
   "450K",    "Bisulfite Conversion I Green",       1,          \(G, R, bkgG, bkgR) cmin(G("BISULFITE CONVERSION I", "^BS Conversion I-C[123]$")) / cmax(G("BISULFITE CONVERSION I", "^BS Conversion I-U[123]$")),
   "450K",    "Bisulfite Conversion I Green (Bkg)", 1,          \(G, R, bkgG, bkgR) bkgG / cmax(G("BISULFITE CONVERSION I", "^BS Conversion I-U[123]$")),
   "450K",    "Bisulfite Conversion I Red",         1,          \(G, R, bkgG, bkgR) cmin(R("BISULFITE CONVERSION I", "^BS Conversion I-C[456]$")) / cmax(R("BISULFITE CONVERSION I", "^BS Conversion I-U[456]$")),
   "450K",    "Bisulfite Conversion I Red (Bkg)",   1,          \(G, R, bkgG, bkgR) bkgR / cmax(R("BISULFITE CONVERSION I", "^BS Conversion I-U[456]$")),
   "450K",    "Bisulfite Conversion II",            1,          \(G, R, bkgG, bkgR) cmin(R("BISULFITE CONVERSION II", "^BS Conversion II-[1-4]$")) / cmax(G("BISULFITE CONVERSION II", "^BS Conversion II-[1-4]$")),
   "450K",    "Bisulfite Conversion II (Bkg)",      1,          \(G, R, bkgG, bkgR) bkgG / cmax(G("BISULFITE CONVERSION II", "^BS Conversion II-[1-4]$")),
   "450K",    "Hybridization High/Medium",          1,          \(G, R, bkgG, bkgR) G("HYBRIDIZATION", "^Hyb \\(High\\)$") / G("HYBRIDIZATION", "^Hyb \\(Medium\\)$"),
   "450K",    "Hybridization Medium/Low",           1,          \(G, R, bkgG, bkgR) G("HYBRIDIZATION", "^Hyb \\(Medium\\)$") / G("HYBRIDIZATION", "^Hyb \\(Low\\)$"),
   "450K",    "Extension Green",                    5,          \(G, R, bkgG, bkgR) cmin(G("EXTENSION", "^Extension \\([CG]\\)$")) / cmax(G("EXTENSION", "^Extension \\([AT]\\)$")),
   "450K",    "Extension Red",                      5,          \(G, R, bkgG, bkgR) cmin(R("EXTENSION", "^Extension \\([AT]\\)$")) / cmax(R("EXTENSION", "^Extension \\([CG]\\)$")),
   "450K",    "Staining Green",                     5,          \(G, R, bkgG, bkgR) G("STAINING", "^Biotin \\(High\\)$") / G("STAINING", "^Biotin \\(Bkg\\)$"),
   "450K",    "Staining Red",                       5,          \(G, R, bkgG, bkgR) R("STAINING", "^DNP \\(High\\)$") / R("STAINING", "^DNP \\(Bkg\\)$"),
   "450K",    "Target Removal 1",                   1,          \(G, R, bkgG, bkgR) bkgG / G("TARGET REMOVAL", "^Target Removal 1$"),
   "450K",    "Target Removal 2",                   1,          \(G, R, bkgG, bkgR) bkgG / G("TARGET REMOVAL", "^Target Removal 2$"),
   "450K",    "Specificity I Green",                1,          \(G, R, bkgG, bkgR) cmin(G("SPECIFICITY I", "^GT Mismatch [123] \\(PM\\)$")) / cmax(G("SPECIFICITY I", "^GT Mismatch [123] \\(MM\\)$")),
   "450K",    "Specificity I Red",                  1,          \(G, R, bkgG, bkgR) cmin(R("SPECIFICITY I", "^GT Mismatch [456] \\(PM\\)$")) / cmax(R("SPECIFICITY I", "^GT Mismatch [456] \\(MM\\)$")),
   "450K",    "Specificity II",                     1,          \(G, R, bkgG, bkgR) cmin(R("SPECIFICITY II", "^Specificity [123]$")) / cmax(G("SPECIFICITY II", "^Specificity [123]$")),
   "450K",    "Specificity II (Bkg)",               1,          \(G, R, bkgG, bkgR) bkgG / cmax(G("SPECIFICITY II", "^Specificity [123]$")),
   "450K",    "Non-polymorphic Green",              5,          \(G, R, bkgG, bkgR) cmin(G("NON-POLYMORPHIC", "^NP \\([CG]\\)$")) / cmax(G("NON-POLYMORPHIC", "^NP \\([AT]\\)$")),
   "450K",    "Non-polymorphic Red",                5,          \(G, R, bkgG, bkgR) cmin(R("NON-POLYMORPHIC", "^NP \\([AT]\\)$")) / cmax(R("NON-POLYMORPHIC", "^NP \\([CG]\\)$")),

   # EPICv1 (21 metrics; BS Conversion I Green = C1-C2/U1-U2, Red = C3-C5/U3-U5 because Illumina
   # dropped 450K's C1 (address 22711390) and U1 (address 46651360) and shifted subsequent probe numbers down by 1)
   "EPICv1",  "Restoration",                        0,          \(G, R, bkgG, bkgR) G("RESTORATION", "^Restore$") / bkgG,
   "EPICv1",  "Bisulfite Conversion I Green",       1,          \(G, R, bkgG, bkgR) cmin(G("BISULFITE CONVERSION I", "^BS Conversion I-C[12]$")) / cmax(G("BISULFITE CONVERSION I", "^BS Conversion I-U[12]$")),
   "EPICv1",  "Bisulfite Conversion I Green (Bkg)", 1,          \(G, R, bkgG, bkgR) bkgG / cmax(G("BISULFITE CONVERSION I", "^BS Conversion I-U[12]$")),
   "EPICv1",  "Bisulfite Conversion I Red",         1,          \(G, R, bkgG, bkgR) cmin(R("BISULFITE CONVERSION I", "^BS Conversion I-C[345]$")) / cmax(R("BISULFITE CONVERSION I", "^BS Conversion I-U[345]$")),
   "EPICv1",  "Bisulfite Conversion I Red (Bkg)",   1,          \(G, R, bkgG, bkgR) bkgR / cmax(R("BISULFITE CONVERSION I", "^BS Conversion I-U[345]$")),
   "EPICv1",  "Bisulfite Conversion II",            1,          \(G, R, bkgG, bkgR) cmin(R("BISULFITE CONVERSION II", "^BS Conversion II-[1-4]$")) / cmax(G("BISULFITE CONVERSION II", "^BS Conversion II-[1-4]$")),
   "EPICv1",  "Bisulfite Conversion II (Bkg)",      1,          \(G, R, bkgG, bkgR) bkgG / cmax(G("BISULFITE CONVERSION II", "^BS Conversion II-[1-4]$")),
   "EPICv1",  "Hybridization High/Medium",          1,          \(G, R, bkgG, bkgR) G("HYBRIDIZATION", "^Hyb \\(High\\)$") / G("HYBRIDIZATION", "^Hyb \\(Medium\\)$"),
   "EPICv1",  "Hybridization Medium/Low",           1,          \(G, R, bkgG, bkgR) G("HYBRIDIZATION", "^Hyb \\(Medium\\)$") / G("HYBRIDIZATION", "^Hyb \\(Low\\)$"),
   "EPICv1",  "Extension Green",                    5,          \(G, R, bkgG, bkgR) cmin(G("EXTENSION", "^Extension \\([CG]\\)$")) / cmax(G("EXTENSION", "^Extension \\([AT]\\)$")),
   "EPICv1",  "Extension Red",                      5,          \(G, R, bkgG, bkgR) cmin(R("EXTENSION", "^Extension \\([AT]\\)$")) / cmax(R("EXTENSION", "^Extension \\([CG]\\)$")),
   "EPICv1",  "Staining Green",                     5,          \(G, R, bkgG, bkgR) G("STAINING", "^Biotin \\(High\\)$") / G("STAINING", "^Biotin \\(Bkg\\)$"),
   "EPICv1",  "Staining Red",                       5,          \(G, R, bkgG, bkgR) R("STAINING", "^DNP \\(High\\)$") / R("STAINING", "^DNP \\(Bkg\\)$"),
   "EPICv1",  "Target Removal 1",                   1,          \(G, R, bkgG, bkgR) bkgG / G("TARGET REMOVAL", "^Target Removal 1$"),
   "EPICv1",  "Target Removal 2",                   1,          \(G, R, bkgG, bkgR) bkgG / G("TARGET REMOVAL", "^Target Removal 2$"),
   "EPICv1",  "Specificity I Green",                1,          \(G, R, bkgG, bkgR) cmin(G("SPECIFICITY I", "^GT Mismatch [123] \\(PM\\)$")) / cmax(G("SPECIFICITY I", "^GT Mismatch [123] \\(MM\\)$")),
   "EPICv1",  "Specificity I Red",                  1,          \(G, R, bkgG, bkgR) cmin(R("SPECIFICITY I", "^GT Mismatch [456] \\(PM\\)$")) / cmax(R("SPECIFICITY I", "^GT Mismatch [456] \\(MM\\)$")),
   "EPICv1",  "Specificity II",                     1,          \(G, R, bkgG, bkgR) cmin(R("SPECIFICITY II", "^Specificity [123]$")) / cmax(G("SPECIFICITY II", "^Specificity [123]$")),
   "EPICv1",  "Specificity II (Bkg)",               1,          \(G, R, bkgG, bkgR) bkgG / cmax(G("SPECIFICITY II", "^Specificity [123]$")),
   "EPICv1",  "Non-polymorphic Green",              5,          \(G, R, bkgG, bkgR) cmin(G("NON-POLYMORPHIC", "^NP \\([CG]\\)( [1-5])?$")) / cmax(G("NON-POLYMORPHIC", "^NP \\([AT]\\)$")),
   "EPICv1",  "Non-polymorphic Red",                5,          \(G, R, bkgG, bkgR) cmin(R("NON-POLYMORPHIC", "^NP \\([AT]\\)$")) / cmax(R("NON-POLYMORPHIC", "^NP \\([CG]\\)( [1-5])?$")),

   # EPICv2 (21 metrics; same control probe structure as EPICv1)
   "EPICv2",  "Restoration",                        0,          \(G, R, bkgG, bkgR) G("RESTORATION", "^Restore$") / bkgG,
   "EPICv2",  "Bisulfite Conversion I Green",       1,          \(G, R, bkgG, bkgR) cmin(G("BISULFITE CONVERSION I", "^BS Conversion I-C[12]$")) / cmax(G("BISULFITE CONVERSION I", "^BS Conversion I-U[12]$")),
   "EPICv2",  "Bisulfite Conversion I Green (Bkg)", 1,          \(G, R, bkgG, bkgR) bkgG / cmax(G("BISULFITE CONVERSION I", "^BS Conversion I-U[12]$")),
   "EPICv2",  "Bisulfite Conversion I Red",         1,          \(G, R, bkgG, bkgR) cmin(R("BISULFITE CONVERSION I", "^BS Conversion I-C[345]$")) / cmax(R("BISULFITE CONVERSION I", "^BS Conversion I-U[345]$")),
   "EPICv2",  "Bisulfite Conversion I Red (Bkg)",   1,          \(G, R, bkgG, bkgR) bkgR / cmax(R("BISULFITE CONVERSION I", "^BS Conversion I-U[345]$")),
   "EPICv2",  "Bisulfite Conversion II",            1,          \(G, R, bkgG, bkgR) cmin(R("BISULFITE CONVERSION II", "^BS Conversion II-[1-4]$")) / cmax(G("BISULFITE CONVERSION II", "^BS Conversion II-[1-4]$")),
   "EPICv2",  "Bisulfite Conversion II (Bkg)",      1,          \(G, R, bkgG, bkgR) bkgG / cmax(G("BISULFITE CONVERSION II", "^BS Conversion II-[1-4]$")),
   "EPICv2",  "Hybridization High/Medium",          1,          \(G, R, bkgG, bkgR) G("HYBRIDIZATION", "^Hyb \\(High\\)$") / G("HYBRIDIZATION", "^Hyb \\(Medium\\)$"),
   "EPICv2",  "Hybridization Medium/Low",           1,          \(G, R, bkgG, bkgR) G("HYBRIDIZATION", "^Hyb \\(Medium\\)$") / G("HYBRIDIZATION", "^Hyb \\(Low\\)$"),
   "EPICv2",  "Extension Green",                    5,          \(G, R, bkgG, bkgR) cmin(G("EXTENSION", "^Extension \\([CG]\\)$")) / cmax(G("EXTENSION", "^Extension \\([AT]\\)$")),
   "EPICv2",  "Extension Red",                      5,          \(G, R, bkgG, bkgR) cmin(R("EXTENSION", "^Extension \\([AT]\\)$")) / cmax(R("EXTENSION", "^Extension \\([CG]\\)$")),
   "EPICv2",  "Staining Green",                     5,          \(G, R, bkgG, bkgR) G("STAINING", "^Biotin \\(High\\)$") / G("STAINING", "^Biotin \\(Bkg\\)$"),
   "EPICv2",  "Staining Red",                       5,          \(G, R, bkgG, bkgR) R("STAINING", "^DNP \\(High\\)$") / R("STAINING", "^DNP \\(Bkg\\)$"),
   "EPICv2",  "Target Removal 1",                   1,          \(G, R, bkgG, bkgR) bkgG / G("TARGET REMOVAL", "^Target Removal 1$"),
   "EPICv2",  "Target Removal 2",                   1,          \(G, R, bkgG, bkgR) bkgG / G("TARGET REMOVAL", "^Target Removal 2$"),
   "EPICv2",  "Specificity I Green",                1,          \(G, R, bkgG, bkgR) cmin(G("SPECIFICITY I", "^GT Mismatch [123] \\(PM\\)$")) / cmax(G("SPECIFICITY I", "^GT Mismatch [123] \\(MM\\)$")),
   "EPICv2",  "Specificity I Red",                  1,          \(G, R, bkgG, bkgR) cmin(R("SPECIFICITY I", "^GT Mismatch [456] \\(PM\\)$")) / cmax(R("SPECIFICITY I", "^GT Mismatch [456] \\(MM\\)$")),
   "EPICv2",  "Specificity II",                     1,          \(G, R, bkgG, bkgR) cmin(R("SPECIFICITY II", "^Specificity [123]$")) / cmax(G("SPECIFICITY II", "^Specificity [123]$")),
   "EPICv2",  "Specificity II (Bkg)",               1,          \(G, R, bkgG, bkgR) bkgG / cmax(G("SPECIFICITY II", "^Specificity [123]$")),
   "EPICv2",  "Non-polymorphic Green",              5,          \(G, R, bkgG, bkgR) cmin(G("NON-POLYMORPHIC", "^NP \\([CG]\\)( [1-5])?$")) / cmax(G("NON-POLYMORPHIC", "^NP \\([AT]\\)$")),
   "EPICv2",  "Non-polymorphic Red",                5,          \(G, R, bkgG, bkgR) cmin(R("NON-POLYMORPHIC", "^NP \\([AT]\\)$")) / cmax(R("NON-POLYMORPHIC", "^NP \\([CG]\\)( [1-5])?$")),

   # MSA (21 metrics; Non-polymorphic has only NP (A/T/C/G))
   "MSA",     "Restoration",                        0,          \(G, R, bkgG, bkgR) G("RESTORATION", "^Restore$") / bkgG,
   "MSA",     "Bisulfite Conversion I Green",       1,          \(G, R, bkgG, bkgR) cmin(G("BISULFITE CONVERSION I", "^BS Conversion I-C[12]$")) / cmax(G("BISULFITE CONVERSION I", "^BS Conversion I-U[12]$")),
   "MSA",     "Bisulfite Conversion I Green (Bkg)", 1,          \(G, R, bkgG, bkgR) bkgG / cmax(G("BISULFITE CONVERSION I", "^BS Conversion I-U[12]$")),
   "MSA",     "Bisulfite Conversion I Red",         1,          \(G, R, bkgG, bkgR) cmin(R("BISULFITE CONVERSION I", "^BS Conversion I-C[45]$")) / cmax(R("BISULFITE CONVERSION I", "^BS Conversion I-U[45]$")),
   "MSA",     "Bisulfite Conversion I Red (Bkg)",   1,          \(G, R, bkgG, bkgR) bkgR / cmax(R("BISULFITE CONVERSION I", "^BS Conversion I-U[45]$")),
   "MSA",     "Bisulfite Conversion II",            1,          \(G, R, bkgG, bkgR) cmin(R("BISULFITE CONVERSION II", "^BS Conversion II-[1-4]$")) / cmax(G("BISULFITE CONVERSION II", "^BS Conversion II-[1-4]$")),
   "MSA",     "Bisulfite Conversion II (Bkg)",      1,          \(G, R, bkgG, bkgR) bkgG / cmax(G("BISULFITE CONVERSION II", "^BS Conversion II-[1-4]$")),
   "MSA",     "Hybridization High/Medium",          1,          \(G, R, bkgG, bkgR) G("HYBRIDIZATION", "^Hyb \\(High\\)$") / G("HYBRIDIZATION", "^Hyb \\(Medium\\)$"),
   "MSA",     "Hybridization Medium/Low",           1,          \(G, R, bkgG, bkgR) G("HYBRIDIZATION", "^Hyb \\(Medium\\)$") / G("HYBRIDIZATION", "^Hyb \\(Low\\)$"),
   "MSA",     "Extension Green",                    5,          \(G, R, bkgG, bkgR) cmin(G("EXTENSION", "^Extension \\([CG]\\)$")) / cmax(G("EXTENSION", "^Extension \\([AT]\\)$")),
   "MSA",     "Extension Red",                      5,          \(G, R, bkgG, bkgR) cmin(R("EXTENSION", "^Extension \\([AT]\\)$")) / cmax(R("EXTENSION", "^Extension \\([CG]\\)$")),
   "MSA",     "Staining Green",                     5,          \(G, R, bkgG, bkgR) G("STAINING", "^Biotin \\(High\\)$") / G("STAINING", "^Biotin \\(Bkg\\)$"),
   "MSA",     "Staining Red",                       5,          \(G, R, bkgG, bkgR) R("STAINING", "^DNP \\(High\\)$") / R("STAINING", "^DNP \\(Bkg\\)$"),
   "MSA",     "Target Removal 1",                   1,          \(G, R, bkgG, bkgR) (cmax(R("EXTENSION", "^Extension \\([AT]\\)$")) + 3000) / R("TARGET REMOVAL", "^Target Removal 1$"),
   "MSA",     "Target Removal 2",                   1,          \(G, R, bkgG, bkgR) (cmax(R("EXTENSION", "^Extension \\([AT]\\)$")) + 3000) / R("TARGET REMOVAL", "^Target Removal 2$"),
   "MSA",     "Specificity I Green",                1,          \(G, R, bkgG, bkgR) cmin(G("SPECIFICITY I", "^GT Mismatch [123] \\(PM\\)$")) / cmax(G("SPECIFICITY I", "^GT Mismatch [123] \\(MM\\)$")),
   "MSA",     "Specificity I Red",                  1,          \(G, R, bkgG, bkgR) cmin(R("SPECIFICITY I", "^GT Mismatch [456] \\(PM\\)$")) / cmax(R("SPECIFICITY I", "^GT Mismatch [456] \\(MM\\)$")),
   "MSA",     "Specificity II",                     1,          \(G, R, bkgG, bkgR) cmin(R("SPECIFICITY II", "^Specificity [123]$")) / cmax(G("SPECIFICITY II", "^Specificity [123]$")),
   "MSA",     "Specificity II (Bkg)",               1,          \(G, R, bkgG, bkgR) bkgG / cmax(G("SPECIFICITY II", "^Specificity [123]$")),
   "MSA",     "Non-polymorphic Green",              5,          \(G, R, bkgG, bkgR) cmin(G("NON-POLYMORPHIC", "^NP \\([CG]\\)$")) / cmax(G("NON-POLYMORPHIC", "^NP \\([AT]\\)$")),
   "MSA",     "Non-polymorphic Red",                5,          \(G, R, bkgG, bkgR) cmin(R("NON-POLYMORPHIC", "^NP \\([AT]\\)$")) / cmax(R("NON-POLYMORPHIC", "^NP \\([CG]\\)$")),

   # MOUSE (21 metrics; BS Conversion II includes mouse-specific probes 330/505/649)
   "MOUSE",   "Restoration",                        0,          \(G, R, bkgG, bkgR) G("RESTORATION", "^Restore$") / bkgG,
   "MOUSE",   "Bisulfite Conversion I Green",       1,          \(G, R, bkgG, bkgR) cmin(G("BISULFITE CONVERSION I", "^BS Conversion I-C[12]$")) / cmax(G("BISULFITE CONVERSION I", "^BS Conversion I-U[12]$")),
   "MOUSE",   "Bisulfite Conversion I Green (Bkg)", 1,          \(G, R, bkgG, bkgR) bkgG / cmax(G("BISULFITE CONVERSION I", "^BS Conversion I-U[12]$")),
   "MOUSE",   "Bisulfite Conversion I Red",         1,          \(G, R, bkgG, bkgR) cmin(R("BISULFITE CONVERSION I", "^BS Conversion I-C[45]$")) / cmax(R("BISULFITE CONVERSION I", "^BS Conversion I-U[45]$")),
   "MOUSE",   "Bisulfite Conversion I Red (Bkg)",   1,          \(G, R, bkgG, bkgR) bkgR / cmax(R("BISULFITE CONVERSION I", "^BS Conversion I-U[45]$")),
   "MOUSE",   "Bisulfite Conversion II",            1,          \(G, R, bkgG, bkgR) cmin(R("BISULFITE CONVERSION II", "^BS Conversion II-\\d+$")) / cmax(G("BISULFITE CONVERSION II", "^BS Conversion II-\\d+$")),
   "MOUSE",   "Bisulfite Conversion II (Bkg)",      1,          \(G, R, bkgG, bkgR) bkgG / cmax(G("BISULFITE CONVERSION II", "^BS Conversion II-\\d+$")),
   "MOUSE",   "Hybridization High/Medium",          1,          \(G, R, bkgG, bkgR) G("HYBRIDIZATION", "^Hyb \\(High\\)$") / G("HYBRIDIZATION", "^Hyb \\(Medium\\)$"),
   "MOUSE",   "Hybridization Medium/Low",           1,          \(G, R, bkgG, bkgR) G("HYBRIDIZATION", "^Hyb \\(Medium\\)$") / G("HYBRIDIZATION", "^Hyb \\(Low\\)$"),
   "MOUSE",   "Extension Green",                    5,          \(G, R, bkgG, bkgR) cmin(G("EXTENSION", "^Extension \\([CG]\\)$")) / cmax(G("EXTENSION", "^Extension \\([AT]\\)$")),
   "MOUSE",   "Extension Red",                      5,          \(G, R, bkgG, bkgR) cmin(R("EXTENSION", "^Extension \\([AT]\\)$")) / cmax(R("EXTENSION", "^Extension \\([CG]\\)$")),
   "MOUSE",   "Staining Green",                     5,          \(G, R, bkgG, bkgR) G("STAINING", "^Biotin \\(High\\)$") / G("STAINING", "^Biotin \\(Bkg\\)$"),
   "MOUSE",   "Staining Red",                       5,          \(G, R, bkgG, bkgR) R("STAINING", "^DNP \\(High\\)$") / R("STAINING", "^DNP \\(Bkg\\)$"),
   "MOUSE",   "Target Removal 1",                   1,          \(G, R, bkgG, bkgR) (cmax(R("EXTENSION", "^Extension \\([AT]\\)$")) + 3000) / R("TARGET REMOVAL", "^Target Removal 1$"),
   "MOUSE",   "Target Removal 2",                   1,          \(G, R, bkgG, bkgR) (cmax(R("EXTENSION", "^Extension \\([AT]\\)$")) + 3000) / R("TARGET REMOVAL", "^Target Removal 2$"),
   "MOUSE",   "Specificity I Green",                1,          \(G, R, bkgG, bkgR) cmin(G("SPECIFICITY I", "^GT Mismatch [123] \\(PM\\)$")) / cmax(G("SPECIFICITY I", "^GT Mismatch [123] \\(MM\\)$")),
   "MOUSE",   "Specificity I Red",                  1,          \(G, R, bkgG, bkgR) cmin(R("SPECIFICITY I", "^GT Mismatch [456] \\(PM\\)$")) / cmax(R("SPECIFICITY I", "^GT Mismatch [456] \\(MM\\)$")),
   "MOUSE",   "Specificity II",                     1,          \(G, R, bkgG, bkgR) cmin(R("SPECIFICITY II", "^Specificity [123]$")) / cmax(G("SPECIFICITY II", "^Specificity [123]$")),
   "MOUSE",   "Specificity II (Bkg)",               1,          \(G, R, bkgG, bkgR) bkgG / cmax(G("SPECIFICITY II", "^Specificity [123]$")),
   "MOUSE",   "Non-polymorphic Green",              5,          \(G, R, bkgG, bkgR) cmin(G("NON-POLYMORPHIC", "^NP \\([CG]\\)( [1-5])?$")) / cmax(G("NON-POLYMORPHIC", "^NP \\([AT]\\)$")),
   "MOUSE",   "Non-polymorphic Red",                5,          \(G, R, bkgG, bkgR) cmin(R("NON-POLYMORPHIC", "^NP \\([AT]\\)$")) / cmax(R("NON-POLYMORPHIC", "^NP \\([CG]\\)( [1-5])?$"))
)
