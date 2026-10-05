library(minfi)
library(ENmix)
library(IlluminaHumanMethylationEPICanno.ilm10b4.hg19)

data_dir <- normalizePath("~/repos/ewastools-data", mustWork = FALSE)

test_that("compatibility with minfi", {

	idat_prefix <- file.path(data_dir, "GSE110554/GSM2998021")
	meth = read.metharray(idat_prefix, extended = TRUE)
	detP = detectionP.minfi(meth)

	ctrls <- getProbeInfo(meth, type="Control")
	ctrls <- ctrls[ctrls$Address %in% featureNames(meth),]
	ctrl_r <- getRed(meth)[ctrls$Address,, drop=FALSE]
	ctrl_g <- getGreen(meth)[ctrls$Address,, drop=FALSE]
	CG.controls <- ctrls$Type %in% c("NORM_C", "NORM_G")
	AT.controls <- ctrls$Type %in% c("NORM_A", "NORM_T")
	cg_grn <- ctrl_g[CG.controls,, drop=FALSE]
	rownames(cg_grn) = ctrls$ExtendedType[CG.controls]
	at_red <- ctrl_r[AT.controls,, drop=FALSE]
	rownames(at_red) = ctrls$ExtendedType[AT.controls]
	meth <- preprocessRaw(meth)
	meth <- ENmix::relic(meth, at_red, cg_grn)

	meth = getBeta(meth)

	expect_equal(rownames(meth), rownames(detP))
})
