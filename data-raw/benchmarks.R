library(data.table)
library(magrittr)
library(purrr)
library(stringi)
library(car)
library(ewastools)
library(parallel)

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
samples[,file:="benchmarks/" %s+% gsm]
samples = split(samples,by="gse")

# ---------------------------------------------------------------------------

estimateLC = function (meth,ref) 
{
	J = ncol(meth)
	coefs = read.table(ref)
	coefs = as.matrix(coefs)
	n_celltypes = ncol(coefs)
	markers = match(rownames(coefs), rownames(meth))
	EST = sapply(1:J, function(j) {
		tmp = meth[markers, j]
		i = !is.na(tmp)
		quadprog::solve.QP(t(coefs[i, ]) %*% coefs[i, ], t(coefs[i, 
		 ]) %*% tmp[i], diag(n_celltypes), rep(0, n_celltypes))$sol
	})
	EST = t(EST)
	colnames(EST) = colnames(coefs)
	EST = data.table(EST)
	return(EST)
}

avg_r2 = function(ref){
	
	LC = estimateLC(meth,ref=ref)
	frml = paste(c("meth~1",names(LC)),collapse="+")
	LC[,meth:=runif(.N)]
	m = lm(frml,data=LC)
	mm = model.matrix(m)
	
	f = function(meth_i)
	{
		m[1:8] = lm.fit(mm,meth_i)
		summary(m)$adj.r.squared
	}
	f = possibly(f,otherwise=NA)
	
	rr = apply(meth,1,f)
	mean(rr,na.rm=TRUE)
	
}

training = list.files(path="../data",pattern="\\.txt$",full.names=TRUE)
training = setdiff(training,c("../data/blood_norm.txt","../data/hkcpgs.txt","../data/saliva.txt"))

# ---------------------------------------------------------------------------
# test set (GSE77797)  https://doi.org/10.1186/s12859-016-0943-7
gse77797 = samples$GSE77797[,list(gsm,file)]
meth = gse77797 %$% file %>% read_idats %>% correct_dye_bias %>% dont_normalize
R2  = mclapply(training,avg_r2,mc.cores=length(training)) %>% unlist

# ---------------------------------------------------------------------------
# 689 Rheumatoid arthritis GSE42861

gse42861 = samples$GSE42861[,list(gsm,file)]
meth = gse42861 %$% file %>% read_idats %>% correct_dye_bias %>% dont_normalize
R3  = mclapply(training,avg_r2,mc.cores=length(training)) %>% unlist

# ---------------------------------------------------------------------------
# 251 2y olds from Gambia GSE99863

gse99863 = samples$GSE99863[,list(gsm,file)]
meth = gse99863 %$% file %>% read_idats %>% correct_dye_bias %>% dont_normalize
R4  = mclapply(training,avg_r2,mc.cores=length(training)) %>% unlist

# ---------------------------------------------------------------------------
# 71 cord blood samples (GSE85042)
gse85042 = samples$GSE85042[,list(gsm,file)]
meth = gse85042 %$% file %>% read_idats %>% correct_dye_bias %>% dont_normalize
R5  = mclapply(training,avg_r2,mc.cores=length(training)) %>% unlist

round(cbind(R2,R3,R4,R5),3)

#                            R2    R3    R4    R5
# "Bakulski.txt"          0.313 0.287 0.444 0.236
# "Bakulski+deGoede.txt"  0.291 0.283 0.452 0.229
# "Bakulski+Gervin.txt"   0.301 0.291 0.446 0.225
# "Bakulski+Lin.txt"      0.303 0.288 0.450 0.232
# "Bakulski+Mill.txt"     0.284 0.301 0.462 0.226
# "Bakulski+Reinius.txt"  0.324 0.298 0.457 0.234
# "Bakulski+Salas.txt"    0.298 0.298 0.487 0.233
# "deGoede.txt"           0.270 0.280 0.478 0.229
# "deGoede+Lin.txt"       0.278 0.281 0.471 0.229
# "deGoede+Mill.txt"      0.260 0.280 0.473 0.225
# "deGoede+Reinius.txt"   0.285 0.288 0.465 0.226
# "deGoede+Salas.txt"     0.286 0.273 0.505 0.232
# "Gervin.txt"            0.277 0.175 0.244 0.125
# "Gervin+deGoede.txt"    0.281 0.286 0.461 0.225
# "Gervin+Lin.txt"        0.284 0.162 0.245 0.128
# "Gervin+Mill.txt"       0.273 0.176 0.288 0.131
# "Gervin+Reinius.txt"    0.289 0.168 0.235 0.127
# "Gervin+Salas.txt"      0.287 0.167 0.249 0.130
# "Lin.txt"               0.278 0.159 0.255 0.129
# "Lin+Mill.txt"          0.274 0.173 0.284 0.137
# "Lin+Reinius.txt"       0.288 0.165 0.230 0.125
# "Lin+Salas.txt"         0.284 0.166 0.242 0.126
# "Mill.txt"              0.257 0.147 0.216 0.122
# "Reinius.txt"           0.283 0.210 0.307 0.122
# "Reinius+Mill.txt"      0.261 0.187 0.294 0.125
# "Reinius+Salas.txt"     0.285 0.159 0.245 0.122
# "Salas.txt"             0.278 0.165 0.245 0.104
# "Salas+Mill.txt"        0.268 0.165 0.289 0.146
# "Lolipop.txt"           0.270 0.108 0.314 0.156


bakulski = fread("saliva/clean_saliva_samples.csv")
bakulski = bakulski[cell_type=="Whole"]

meth = bakulski %$% file %>% read_idats %>% correct_dye_bias %>% dont_normalize

estimateLC(meth,"../data/saliva.txt")
estimateLC(meth,"../data/Bakulski.txt")

avg_r2("../data/saliva.txt")
# 0.2507195
avg_r2("../data/bakulski.txt")
# 0.371119

## Comparison to EpiDISH
library(EpiDISH)
  
LC = epidish(beta.m=meth,ref.m=centEpiFibIC.m,method="RPC")$estF
LC %<>% data.table
LC = LC[,.(Epi,IC)]

frml = paste(c("meth~1",names(LC)),collapse="+")
LC[,meth:=runif(.N)]
m = lm(frml,data=LC)
mm = model.matrix(m)

f = function(meth_i)
{
	m[1:8] = lm.fit(mm,meth_i)
	summary(m)$adj.r.squared
}
f = possibly(f,otherwise=NA)

rr = apply(meth,1,f)
mean(rr,na.rm=TRUE)
# 0.2137797
