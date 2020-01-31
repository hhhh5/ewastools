library(minfi)
library(ENmix)
library(IlluminaHumanMethylationEPICanno.ilm10b4.hg19)

test_that("compatibility with minfi", {

	idat1 = system.file(paste0("data/200144450018_R04C01_Grn.idat"),package="ewastools")
	# the horror: test can only work on Windows at the moment, but I don't want to add another bulky file.
	idat2 = system.file(paste0("data/200144450018_r04c01_Grn.idat"),package="ewastools")
	meth = read.metharray(c(idat1,idat2),extended=TRUE)
	detP = detectionP.minfi(meth)

	ctrls <- getProbeInfo(meth,type="Control")
	ctrls <- ctrls[ctrls$Address %in% featureNames(meth),]
	ctrl_r <- getRed(meth)[ctrls$Address,,drop=FALSE]
	ctrl_g <- getGreen(meth)[ctrls$Address,,drop=FALSE]
	CG.controls <- ctrls$Type %in% c("NORM_C","NORM_G")
	AT.controls <- ctrls$Type %in% c("NORM_A","NORM_T")
	cg_grn <- ctrl_g[CG.controls,,drop=FALSE]
	rownames(cg_grn) = ctrls$ExtendedType[CG.controls]
	at_red <- ctrl_r[AT.controls,,drop=FALSE]
	rownames(at_red) = ctrls$ExtendedType[AT.controls]
	meth <- preprocessRaw(meth)
	meth <- ENmix::relic(meth,at_red,cg_grn)

	meth = getBeta(meth)

	expect_equal(rownames(meth),rownames(detP))

})
