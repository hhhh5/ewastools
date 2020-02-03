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

training = list.files("data")
training = stri_sub(training,1,-5)

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

#                      R2    R3    R4    R5; +nRBC    R2    R3    R4    R5
# Bakulski+deGoede  0.284 0.187 0.266 0.146; +nRBC 0.293 0.285 0.449 0.229
# Bakulski+Gervin   0.291 0.177 0.250 0.128; +nRBC 0.308 0.292 0.455 0.225
# Bakulski+Lin      0.288 0.178 0.257 0.139; +nRBC 0.308 0.290 0.450 0.231
# Bakulski+Mill     0.278 0.180 0.279 0.106; +nRBC 0.285 0.304 0.469 0.226
# Bakulski+Reinius  0.303 0.182 0.265 0.143; +nRBC 0.322 0.298 0.463 0.233
# Bakulski+Salas    0.310 0.180 0.281 0.135; +nRBC 0.317 0.303 0.489 0.233
# Bakulski          0.293 0.230 0.299 0.134; +nRBC 0.316 0.286 0.453 0.233
# deGoede+Lin       0.283 0.162 0.265 0.129; +nRBC 0.297 0.285 0.457 0.229
# deGoede+Mill      0.268 0.182 0.320 0.117; +nRBC 0.272 0.286 0.468 0.226
# deGoede+Reinius   0.280 0.166 0.235 0.124; +nRBC 0.297 0.293 0.465 0.227
# deGoede+Salas     0.287 0.173 0.265 0.120; +nRBC 0.269 0.289 0.481 0.233
# deGoede           0.282 0.178 0.254 0.117; +nRBC 0.283 0.282 0.464 0.228
# Gervin+deGoede    0.287 0.168 0.248 0.119; +nRBC 0.287 0.288 0.449 0.226
# Gervin+Lin        0.286 0.168 0.258 0.118;
# Gervin+Mill       0.271 0.178 0.291 0.117;
# Gervin+Reinius    0.288 0.172 0.240 0.129;
# Gervin+Salas      0.286 0.167 0.251 0.121;
# Gervin            0.281 0.183 0.246 0.121;
# Lin+Mill          0.273 0.178 0.290 0.134;
# Lin+Reinius       0.290 0.172 0.245 0.131;
# Lin+Salas         0.286 0.172 0.258 0.118;
# Lin               0.281 0.164 0.258 0.131;
# Mill              0.260 0.149 0.221 0.118;
# Reinius+Mill      0.261 0.191 0.298 0.132;
# Reinius+Salas     0.287 0.168 0.260 0.121;
# Reinius           0.281 0.210 0.312 0.137;
# Salas+Mill        0.268 0.170 0.293 0.143;
# Salas             0.279 0.172 0.255 0.107;
