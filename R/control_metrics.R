#' @import data.table
NULL

#' Calculate the QC metrics as described in the 'BeadArray Controls Reporter Software Guide' from Illumina.
#' 
#' @export
#' @note You can download the Software Guide at https://support.illumina.com/downloads/beadarray-controls-reporter-software-guide-1000000004009.html
#' @title Quality control metrics
#' 
#' @param raw Output of calling \code{\link{read_idats}}
#' @param metrics Output of callling \code{control_metrics}
#' 
#' @return For \code{control_metrics}, a list of 17 control metrics, each one a numeric vector equal in length to the sample size
#' @return For \code{sample_failure}, a logical vector, \code{TRUE} if sample at corresponding index failed on any of the 17 control metrics.
#' 
control_metrics = function(raw){

   if(!all(c("controls", "ctrlG", "ctrlR", "meta") %in% names(raw))) stop("Invalid argument")

   with(raw,{

      ### calculate the control probe metrics
      metrics = list()

      cg = controls[group == "EXTENSION" & name %like% "\\([CG]\\)|EXT-[CG]$", index]
      at = controls[group == "EXTENSION" & name %like% "\\([AT]\\)|EXT-[AT]$", index]

      # Restoration
      ii  = controls[group == 'RESTORATION', index]
      metrics$Restoration = as.numeric(ctrlG[ii,,drop=FALSE] / (apply(ctrlG[at,,drop=FALSE],2,max)+3000))
      attr(metrics$Restoration, "threshold") <- 0

      # Staining
      ii  = controls[name %ilike% "biotin.*high", index]
      bkg = controls[name %ilike% "biotin.*bkg" , index]
      metrics$`Staining Green` = ctrlG[ii,] / ctrlG[bkg,]
      attr(metrics$`Staining Green`, "threshold") <- 5

      ii  = controls[name %ilike% "DNP.*high", index]
      bkg = controls[name %ilike% "DNP.*bkg",  index]
      metrics$`Staining Red` = ctrlR[ii,] / ctrlR[bkg,]
      attr(metrics$`Staining Red`, "threshold") <- 5

      # Extension
      metrics$`Extension Green` = apply(ctrlG[cg,,drop=FALSE],2,min) / apply(ctrlG[at,,drop=FALSE],2,max) 
      metrics$`Extension Red`   = apply(ctrlR[at,,drop=FALSE],2,min) / apply(ctrlR[cg,,drop=FALSE],2,max) 
      attr(metrics$`Extension Green`,"threshold") <- 5
      attr(metrics$`Extension Red`,  "threshold") <- 5

      # Hybridization
      hyb_l = controls[name %ilike% "hyb.*low"   , index]
      hyb_m = controls[name %ilike% "hyb.*medium", index]
      hyb_h = controls[name %ilike% "hyb.*high"  , index]
      metrics$`Hybridization High/Medium` = ctrlG[hyb_h,] / ctrlG[hyb_m,]
      metrics$`Hybridization Medium/Low`  = ctrlG[hyb_m,] / ctrlG[hyb_l,]
      attr(metrics$`Hybridization High/Medium`, "threshold") <- 1
      attr(metrics$`Hybridization Medium/Low`,  "threshold") <- 1

      # Target removal
      ii = controls[name %in% c("Target Removal 1", "TRM-1"), index]
      metrics$`Target Removal 1` =  (apply(ctrlG[at,,drop=FALSE],2,max)+3000) / ctrlG[ii,]
      ii = controls[name %in% c("Target Removal 2", "TRM-2"), index]
      metrics$`Target Removal 2` =  (apply(ctrlG[at,,drop=FALSE],2,max)+3000) / ctrlG[ii,]
      attr(metrics$`Target Removal 1`, "threshold") <- 1
      attr(metrics$`Target Removal 2`, "threshold") <- 1

      # Bisulfite conversion I (no I-C6/I-U6 probes for EPIC but 450K)
      ii  = controls[group == "BISULFITE CONVERSION I" & name %like% "[I1].C[12]", index] # 450K: "I C1"; EPIC "I-C1", MOUSE "1-C1"
      bkg = controls[group == "BISULFITE CONVERSION I" & name %like% "[I1].U[12]", index] # 450K: "I U1"; EPIC "I-U1", MOUSE "1-U1"
      metrics$`Bisulfite Conversion I Green` = apply(ctrlG[ii,,drop=FALSE],2,min) / apply(ctrlG[bkg,,drop=FALSE],2,max)
      attr(metrics$`Bisulfite Conversion I Green`,'threshold') <- 1
      
      # Background/(U1, U2, or U3) [manual doesn't specify highest/lowest. I chose `max` mirroring the calculation for red]
      # Green channel-bkg = Extension Green highest AT
      ii  = controls[group == "BISULFITE CONVERSION I" & name %like% "[I1].U[12]", index]
      metrics$`Bisulfite Conversion I Green (Bkg)` = (apply(ctrlG[at,,drop=FALSE],2,max)+3000) / apply(ctrlG[ii,,drop=FALSE],2,max)
      attr(metrics$`Bisulfite Conversion I Green (Bkg)`,'threshold') <- 1

      # Lowest value of C4,5, or 6 / Highest value of U4, 5, or 6
      ii  = controls[group == "BISULFITE CONVERSION I" & name %like% "[I1].C[456]", index]
      bkg = controls[group == "BISULFITE CONVERSION I" & name %like% "[I1].U[456]", index] 
      metrics$`Bisulfite Conversion I Red` = apply(ctrlR[ii,,drop=FALSE],2,min) / apply(ctrlR[bkg,,drop=FALSE],2,max)
      attr(metrics$`Bisulfite Conversion I Red`,'threshold') <- 1

      # Background /(Highest value of U4, U5, or U6)
      # Red Channel-bkg = Extension Red highest CG
      ii  = controls[group == "BISULFITE CONVERSION I" & name %like% "[I1].U[456]", index]
      metrics$`Bisulfite Conversion I Red (Bkg)` = (apply(ctrlR[cg,,drop=FALSE],2,max)+3000) / apply(ctrlR[ii,,drop=FALSE],2,max)
      attr(metrics$`Bisulfite Conversion I Red (Bkg)`,'threshold') <- 1

      # Bisulfite conversion II
      # (Lowest of red C1, 2, 3, or 4) / (Highest of green C1, 2, 3, or 4)
      ii  = controls[group == "BISULFITE CONVERSION II", index]
      metrics$`Bisulfite Conversion II` = apply(ctrlR[ii,,drop=FALSE],2,min) / apply(ctrlG[ii,,drop=FALSE],2,max)
      attr(metrics$`Bisulfite Conversion II`, "threshold") <- 1

      # Bisulfite conversion II (Bkg)
      # Green channel-bkg = Extension Green highest AT
      # Background/(Highest C1, C2, C3, or C4 green).
      ii = controls[group == "BISULFITE CONVERSION II", index]
      metrics$`Bisulfite Conversion II (Bkg)` = (apply(ctrlG[at,,drop=FALSE],2,max)+3000) / apply(ctrlG[ii,,drop=FALSE],2,max)
      attr(metrics$`Bisulfite Conversion II (Bkg)`, "threshold") <- 1

      # Specificity I
      pm  = controls[group == "SPECIFICITY I" & name %like% "Mismatch [123] \\(PM\\)|-PM[123]", index]
      mm  = controls[group == "SPECIFICITY I" & name %like% "Mismatch [123] \\(MM\\)|-MM[123]", index]
      metrics$`Specificity I Green` = apply(ctrlG[pm,,drop=FALSE],2,min) / apply(ctrlG[mm,,drop=FALSE],2,max)
      pm  = controls[group == "SPECIFICITY I" & name %like% "Mismatch [456] \\(PM\\)|-PM[456]", index]
      mm  = controls[group == "SPECIFICITY I" & name %like% "Mismatch [456] \\(MM\\)|-MM[456]", index]
      metrics$`Specificity I Red` = apply(ctrlR[pm,,drop=FALSE],2,min) / apply(ctrlR[mm,,drop=FALSE],2,max)
      attr(metrics$`Specificity I Green`,'threshold') <- 1
      attr(metrics$`Specificity I Red`,'threshold') <- 1

      # Specificity II
      # (Lowest intensity of S1, S2, or S3 red) / (Highest intensity of S1, S2, or S3 green)
      ii  = controls[group == "SPECIFICITY II", index]
      metrics$`Specificity II` = apply(ctrlR[ii,,drop=FALSE],2,min) / apply(ctrlG[ii,,drop=FALSE],2,max)
      attr(metrics$`Specificity II`,'threshold') <- 1

      # Specificity II (Bkg)
      # Background/(Highest intensity S1, S2, S3, or S4 green)
      # bkg = Extension Green highest A or T intensity
      metrics$`Specificity II (Bkg)` = (apply(ctrlG[at,,drop=FALSE],2,max)+3000) / apply(ctrlG[ii,,drop=FALSE],2,max)
      attr(metrics$`Specificity II (Bkg)`,'threshold') <- 1

      # Non-polymorphic
      cg  = controls[group == "NON-POLYMORPHIC" & name %like% "NP \\([CG]\\)$|NPM-[CG]", index]
      at  = controls[group == "NON-POLYMORPHIC" & name %like% "NP \\([AT]\\)$|NPM-[AT]", index]
      metrics$`Non-polymorphic Green` = apply(ctrlG[cg,,drop=FALSE],2,min) / apply(ctrlG[at,,drop=FALSE],2,max)
      metrics$`Non-polymorphic Red`   = apply(ctrlR[at,,drop=FALSE],2,min) / apply(ctrlR[cg,,drop=FALSE],2,max)
      attr(metrics$`Non-polymorphic Green`,'threshold') <- 5
      attr(metrics$`Non-polymorphic Red`,  'threshold') <- 5

      return(metrics)
   })
}

#' @rdname control_metrics
#' @export
#' 
sample_failure = function(metrics){
   failed = sapply(metrics,function(metric){
      metric < attr(metric,"threshold") 
   })

   apply(failed,1,any,na.rm=TRUE)
}
