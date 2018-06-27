library(data.table)
library(magrittr)
library(purrr)
library(stringi)
library(car)
library(ewastools)


samples = c("GSE77797","GSE99863","GSE42861","GSE85042")
samples = "https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=" %s+% samples %s+% "&targ=gsm&form=text&view=brief"
samples %<>% map(readLines) %>% unlist
samples %<>% split(.,cumsum(. %like% "^\\^SAMPLE = GSM") )

# Extract GSM accessions
names(samples) = map(samples,1) %>% stri_match_first(regex="GSM\\d+")

# Parse meta data
imap(samples,function(s,acc){
	s = strsplit(s,split=" = ",fixed=TRUE)	
	data.table(gsm=acc,variable=map_chr(s,1),value=map_chr(s,2))
}) -> samples

samples = rbindlist(samples)

# Keep only information on sample characteristics and supplementary files
samples = samples[variable %chin% c("!Sample_characteristics_ch1","!Sample_supplementary_file","!Sample_title","!Sample_series_id")]
samples[variable=="!Sample_title",variable:="name"]
samples[variable=="!Sample_series_id",variable:="gse"]

# Find the URLs pointing to the two .idat files
samples[variable=="!Sample_supplementary_file" & value %like% "_Red\\.idat",variable:="red"]
samples[variable=="!Sample_supplementary_file" & value %like% "_Grn\\.idat",variable:="grn"]

i = samples[variable == "!Sample_characteristics_ch1",which=TRUE]
ch = samples$value[i] %>% stri_split(fixed=": ")
samples$variable[i] = map_chr(ch,1)
samples$value   [i] = map_chr(ch,2)
rm(ch)

samples[,variable:= recode(variable,"'b cell (%)'='B';'cd4+ t cell (%)'='CD4';'cd8+ t cell (%)'='CD8';'granulocyte (%)'='GR';'monocyte (%)'='MO';'natural killer cell (%)'='NK'")]
samples[,variable:= recode(variable,"'gender'='sex';'Sex'='sex'")]

# Reshape data.table from long to wide format
samples = dcast(samples, gsm ~ variable)

map2(samples$red, samples$gsm %s+% "_Red.idat.gz", ~ download.file(.x,.y) ) %>% invisible
map2(samples$grn, samples$gsm %s+% "_Grn.idat.gz", ~ download.file(.x,.y) ) %>% invisible

samples = split(samples,by="gse")

# ---------------------------------------------------------------------------
# test set (GSE77797)  https://doi.org/10.1186/s12859-016-0943-7
gse77797 = samples$GSE77797[,list(gsm,B,CD4,CD8,GR,MO,NK,red,grn)]
gse77797$NK  %<>% as.numeric
gse77797$GR  %<>% as.numeric
gse77797$MO  %<>% as.numeric
gse77797$CD8 %<>% as.numeric
gse77797$CD4 %<>% as.numeric
gse77797$B   %<>% as.numeric

meth = gse77797 %$% gsm %>% read_idats %>% correct_dye_bias %>% dont_normalize

cell_types = c("B","CD4","CD8","GR","MO","NK")
training = c("reinius","gervin","gervin+reinius","goede","bakulski","bakulski+goede","bakulski+gervin","bakulski+reinius","goede+reinius")

R1 = sapply(training,function(ref){

	LC = estimateLC(meth,ref=ref)
	sapply(cell_types,function(ct){
		cor(gse77797[[ct]],LC[[ct]])
	})
})

R1 = t(R1)

R2 = sapply(training,function(ref){

	cat(ref)
	LC = estimateLC(meth,ref=ref)
	LC[,meth:=runif(.N)]
	m = lm(meth ~ B+CD4+CD8+GR+MO+NK,data=LC)
	mm = model.matrix(m)

	mean(apply(na.omit(meth),1,function(meth_i){
		m[1:8] = lm.fit(mm,meth_i)
		summary(m)$adj.r.squared
	}))
})

# ---------------------------------------------------------------------------
# 689 Rheumatoid arthritis GSE42861

gse42861 = samples$GSE42861[,list(gsm,sex,red,grn)]
gse42861[,meth:=runif(.N)]
gse42861[,sex:=factor(sex)]

meth = gse42861 %$% gsm %>% read_idats %>% correct_dye_bias %>% dont_normalize

R3 = sapply(training,function(ref){

	cat(ref)
	LC = estimateLC(meth,ref=ref)
	LC[,meth:=runif(.N)]
	m = lm(meth ~ B+CD4+CD8+GR+MO+NK,data=LC)
	mm = model.matrix(m)

	mean(apply(na.omit(meth),1,function(meth_i){
		m[1:8] = lm.fit(mm,meth_i)
		summary(m)$adj.r.squared
	}))
})

# ---------------------------------------------------------------------------
# 251 2y olds from Gambia GSE99863

gse99863 = samples$GSE99863[,list(gsm,sex,red,grn)]
gse99863[,meth:=runif(.N)]
gse99863[,sex:=factor(sex)]

meth = gse99863 %$% gsm %>% read_idats %>% correct_dye_bias %>% dont_normalize

R4 = sapply(training,function(ref){

	cat(ref)
	LC = estimateLC(meth,ref=ref)
	LC[,meth:=runif(.N)]
	m = lm(meth ~ B+CD4+CD8+GR+MO+NK,data=LC)
	mm = model.matrix(m)

	mean(apply(na.omit(meth),1,function(meth_i){
		m[1:8] = lm.fit(mm,meth_i)
		summary(m)$adj.r.squared
	}))
})


# ---------------------------------------------------------------------------
# 71 cord blood samples (GSE85042)
gse85042 = samples$GSE85042[,list(gsm,red,grn)]
gse85042[,meth:=runif(.N)]

meth = gse85042 %$% gsm %>% read_idats %>% correct_dye_bias %>% dont_normalize

R5 = sapply(training,function(ref){

	cat(ref)
	LC = estimateLC(meth,ref=ref)
	LC[,meth:=runif(.N)]
	m = lm(meth ~ B+CD4+CD8+GR+MO+NK,data=LC)
	mm = model.matrix(m)

	mean(apply(na.omit(meth),1,function(meth_i){
		m[1:8] = lm.fit(mm,meth_i)
		summary(m)$adj.r.squared
	}))
})




#                      B   CD4   CD8    GR    MO    NK        R2        R3        R5        R4
# reinius          0.995 0.936 0.986 0.998 0.981 0.830 0.2823696 0.1741940 0.1339679 0.2587702
# gervin           0.996 0.751 0.854 0.997 0.983 0.942 0.2762569 0.1582688 0.1108876 0.2387656
# gervin+reinius   0.996 0.796 0.883 0.998 0.982 0.945 0.2864695 0.1670512 0.1186851 0.2339372
# goede            0.994 0.845 0.972 0.998 0.980 0.933 0.2716647 0.2041458 0.1868533 0.4210794
# bakulski         0.997 0.799 0.926 0.998 0.983 0.596 0.2904914 0.2269359 0.2058096 0.4206478
# bakulski+goede   0.996 0.855 0.975 0.998 0.980 0.925 0.2813881 0.2409419 0.2025652 0.4231429 0.4529437 (+nRBC)
# bakulski+gervin  0.997 0.675 0.927 0.998 0.987 0.936 0.2971072 0.1653110 0.1563384 0.2376338
# bakulski+reinius 0.996 0.899 0.976 0.998 0.977 0.770 0.2865207 0.1719765 0.1682092 0.2305330
# goede+reinius    0.996 0.875 0.979 0.998 0.982 0.946 0.2828753 0.1651178 0.1298767 0.2294623 
