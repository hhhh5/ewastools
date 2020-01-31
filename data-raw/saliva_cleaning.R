library(data.table)
library(ewastools)
library(magrittr)
library(svd)
library(ggplot2)
library(purrr)

meta = readRDS("saliva/all_samples.rds")
setDT(meta)

meta = meta[,.(
	 donor = id
	,sex = factor(Male_label)
	,age
	,cell_type = celltype
	,barcode = Sentrix_ID
	,position = Sentrix_Position
	,rep = ifelse(NOTE=="technical duplicate",TRUE,FALSE)
	,White,Hispanic,Black,Asian,Other
	,Sick
	)]

meta[,file:=paste0("~/Jdrive/PM/Just_Lab/projects/saliva/data-raw/saliva/",barcode,"_",position)]
meta[,sample_id:=paste0(barcode,"_",position)]
meta[,j:=1:.N]

meth = read_idats(meta$file)

ctrl = control_metrics(meth)
meta$failed = sample_failure(ctrl)

meth %<>% detectionP
chrY = meth$manifest$chr == "Y"
meta$detP = colSums(meth$detP[-chrY,] > 0.05,na.rm=TRUE)

meta[,c("X","Y"):=check_sex(meth)]

plot(Y~X,meta,col=ifelse(meta$sex=="m",1,2),pch=ifelse(meta$sex=="m",3,1))

meta[,predicted_sex:=predict_sex(X,Y,which(sex=="m"),which(sex=="f"))]

# ---------------------- PCA

beta = meth %>% dont_normalize
set.seed(982278)
chrXY = meth$manifest[chr%in%c("X","Y") & probe_type!="rs",index]
j = meta[failed==FALSE & cell_type%in%c("large","CD45pos","whole")]$j
pcs = beta[-chrXY,j]
pcs = pcs - rowMeans(pcs)
pcs = na.omit(pcs)
pcs = t(pcs)
pcs = trlan.svd(pcs,neig=2) # compute the first two principal components
pcs = pcs$u

meta$pc1[j] = pcs[,1]
meta$pc2[j] = pcs[,2]

ggplot(meta[j]) + geom_text(aes(x=pc1,y=pc2,label=donor,color=cell_type),size=3)

meta$pca = FALSE
meta[sample_id %in% c("203287590243_R03C01","203287590025_R05C01","203287590228_R05C01","203293440009_R05C01"),pca:=TRUE]

# ----------------------

beta = meth %>% mask(0.05) %>% dont_normalize

snps = meth$manifest$probe_type == "rs"
snps = beta[snps,]

gt = call_genotypes(snps,learn=FALSE)

meta$snp_outlier = snp_outliers(gt)

stripchart(meta$snp_outlier,m="j")

j = meta$snp_outlier < -4
gt = call_genotypes(snps[,j],learn=FALSE)
check_snp_agreement(gt,meta$donor[j],meta$sample_id[j])

# $`1`
#    donor1             sample1 donor2             sample2 agreement
# 1:     18 203282450208_R05C01     22 203287590025_R05C01 0.9999988
# 2:     18 203282450209_R04C01     22 203287590025_R05C01 0.9999912
# 3:     22 203286230053_R06C01     22 203287590025_R05C01 0.3954723
# 4:     22 203286230103_R05C01     22 203287590025_R05C01 0.3970531
# 5:     18 203286230120_R03C01     22 203287590025_R05C01 0.9999996
# 6:     22 203287590025_R05C01     18 203287590228_R01C01 0.9999994
# 7:     22 203287590025_R05C01     18 203287590228_R08C01 0.9999974

# ----------------------

ex = meta[failed==TRUE | pca==TRUE | snp_outlier > -2 | sex!=predicted_sex | detP > 2e6]
meta = meta[-ex$j]

meta[cell_type=="CD45pos",cell_type:="Leukocytes" ]
meta[cell_type=="large"  ,cell_type:="Epithelial cells"]
meta[cell_type=="whole"  ,cell_type:="Whole"]

meta = meta[,.(
	 sample_id=paste0(barcode,"_",position)
	,donor
	,sex
	,cell_type
	,tissue="Saliva"
	,study="Bakulski"
	)]

meta[,file:=paste0("saliva/",sample_id)]

write.csv(meta,"saliva/clean_saliva_samples.csv",row.names=FALSE)
