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

	# if(!all(c('controls','ctrlG','ctrlR','meta','platform')%in%names(raw))) stop('Invalid argument')
	if(!all(c('controls','ctrlG','ctrlR','meta')%in%names(raw))) stop('Invalid argument')

	with(raw,{

		### calculate the control probe metrics
		### Note: color channels are switched for the 450K platform compared to EPIC and the manual
		metrics = list()

		# Restoration
		ii  = controls[name=='Restore',index]
		bkg = controls[name%like%'Extension \\([AT]\\)',index]
		metrics$Restoration = ctrlG[ii,,drop=FALSE] / (apply(ctrlG[bkg,,drop=FALSE],2,max)+3000)
		attr(metrics$Restoration,'threshold') <- 0

		# Staining
		ii  = controls[name=='Biotin (High)',index]
		bkg = controls[name=='Biotin (Bkg)' ,index]
		metrics$`Staining Green` = ctrlG[ii,,drop=FALSE]/ctrlG[bkg,,drop=FALSE]
		attr(metrics$`Staining Green`,'threshold') <- 5

		ii  = controls[name=='DNP (High)',index]
		bkg = controls[name=='DNP (Bkg)' ,index]
		metrics$`Staining Red` = ctrlR[ii,,drop=FALSE]/ctrlR[bkg,,drop=FALSE]
		attr(metrics$`Staining Red`,'threshold') <- 5

		# Extension
		cg = controls[name%like%'Extension \\([CG]\\)',index]
		at = controls[name%like%'Extension \\([AT]\\)',index]
		metrics$`Extension Green` = apply(ctrlG[cg,,drop=FALSE],2,min) / apply(ctrlG[at,,drop=FALSE],2,max) 
		metrics$`Extension Red` = apply(ctrlR[at,,drop=FALSE],2,min) / apply(ctrlR[cg,,drop=FALSE],2,max) 
		attr(metrics$`Extension Green`,'threshold') <- 5
		attr(metrics$`Extension Red`,'threshold') <- 5

		# Hybridization
		hyb_l = controls[name=='Hyb (Low)'   ,index]
		hyb_m = controls[name=='Hyb (Medium)',index]
		hyb_h = controls[name=='Hyb (High)'  ,index]
		metrics$`Hybridization High/Medium` = ctrlG[hyb_h,,drop=FALSE]/ctrlG[hyb_m,,drop=FALSE] 
		metrics$`Hybridization Medium/Low` = ctrlG[hyb_m,,drop=FALSE]/ctrlG[hyb_l,,drop=FALSE]
		attr(metrics$`Hybridization High/Medium`,'threshold') <- 1
		attr(metrics$`Hybridization Medium/Low`,'threshold') <- 1

		# Target removal
		bkg = controls[name%like%'Extension \\([AT]\\)',index]
		ii = controls[name=='Target Removal 1',index]
		metrics$`Target Removal 1` =  (apply(ctrlG[bkg,,drop=FALSE],2,max)+3000) / ctrlG[ii,,drop=FALSE]
		ii = controls[name=='Target Removal 2',index]
		metrics$`Target Removal 2` =  (apply(ctrlG[bkg,,drop=FALSE],2,max)+3000) / ctrlG[ii,,drop=FALSE]
		attr(metrics$`Target Removal 1`,'threshold') <- 1
		attr(metrics$`Target Removal 2`,'threshold') <- 1

		# Bisulfite conversion I
		ii  = controls[name%like%'I-C[123]',index]
		bkg = controls[name%like%'I-U[123]',index]
		metrics$`Bisulfite Conversion I Green` = apply(ctrlG[ii,,drop=FALSE],2,min) / apply(ctrlG[bkg,,drop=FALSE],2,max)
		ii  = controls[name%like%'I-C[456]',index]
		bkg = controls[name%like%'I-U[456]',index] 
		metrics$`Bisulfite Conversion I Red` = apply(ctrlR[ii,,drop=FALSE],2,min) / apply(ctrlR[bkg,,drop=FALSE],2,max) 
		attr(metrics$`Bisulfite Conversion I Green`,'threshold') <- 1
		attr(metrics$`Bisulfite Conversion I Red`,'threshold') <- 1

		# Bisulfite conversion II
		ii  = controls[group=='BISULFITE CONVERSION II',index]
		metrics$`Bisulfite Conversion II` = apply(ctrlR[ii,,drop=FALSE],2,min) / apply(ctrlG[ii,,drop=FALSE],2,max)
		attr(metrics$`Bisulfite Conversion II`,'threshold') <- 1

		# Specificity I
		pm  = controls[name%like%'Mismatch [123] \\(PM\\)',index]
		mm  = controls[name%like%'Mismatch [123] \\(MM\\)',index]
		metrics$`Specificity I Green` = apply(ctrlG[pm,,drop=FALSE],2,min) / apply(ctrlG[mm,,drop=FALSE],2,max)
		pm  = controls[name%like%'Mismatch [456] \\(PM\\)',index]
		mm  = controls[name%like%'Mismatch [456] \\(MM\\)',index]
		metrics$`Specificity I Red` = apply(ctrlR[pm,,drop=FALSE],2,min) / apply(ctrlR[mm,,drop=FALSE],2,max)
		attr(metrics$`Specificity I Green`,'threshold') <- 1
		attr(metrics$`Specificity I Red`,'threshold') <- 1

		# Specificity II
		ii  = controls[group=='SPECIFICITY II',index]
		metrics$`Specificity II` = apply(ctrlR[ii,,drop=FALSE],2,min) / apply(ctrlG[ii,,drop=FALSE],2,max)
		attr(metrics$`Specificity II`,'threshold') <- 1

		# Non-polymorphic
		cg  = controls[name%like%'NP \\([CG]\\)$',index]
		at  = controls[name%like%'NP \\([AT]\\)$',index]
		metrics$`Non-polymorphic Green` = apply(ctrlG[cg,,drop=FALSE],2,min) / apply(ctrlG[at,,drop=FALSE],2,max)
		metrics$`Non-polymorphic Red` = apply(ctrlR[at,,drop=FALSE],2,min) / apply(ctrlR[cg,,drop=FALSE],2,max)
		attr(metrics$`Non-polymorphic Green`,'threshold') <- 5
		attr(metrics$`Non-polymorphic Red`,'threshold') <- 5

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
