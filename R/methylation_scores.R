#' Various methylation scores
#' 
#' 
#' @rdname "Methylation scores"
#' @param beta matrix of beta values
#' @param model Name of the model, i.e., one of `maternal_smoking`, `gestational_age`, or `horvath_clock`
#' 
methylation_score = function(beta,model){
	if(model=="cord_blood:maternal_smoking"){

		### model based on Table S5 from Reese et al. (2017) available at doi.org/10.1289/EHP333
		markers = data.table(
			 probe_id = c("cg02256631","cg02482603","cg04103532","cg04180046","cg04506190","cg05549655","cg05575921","cg08698721","cg09743950","cg10799846","cg12186702","cg13834112","cg13893782","cg14179389","cg14351425","cg14633298","cg14743346","cg17397069","cg19381766","cg22154659","cg22802102","cg23304605","cg25189904","cg25949550","cg26764244","cg27291468")
			,coef = c(-0.191,2.706,1.786,14.027,2.318,6.210,-10.909,1.142,-6.330,-4.963,3.847,1.514,-0.963,-6.304,6.361,5.050,2.286,-2.912,5.245,-0.773,-0.254,0.011,-3.903,-46.991,-0.246,0.836)
			)

		i = match(markers$probe_id,rownames(beta))
		y = beta[i,]
		y = y * markers$coef
		y = colSums(y,na.rm=TRUE)
		return(y)

	}else if(model=="cord_blood:gestational_age"){

		### doi.org/10.1186/s13059-016-1068-z
		markers = system.file(paste0("data/13059_2016_1068_MOESM3_ESM.csv"),package="ewastools")
		markers = fread(markers)
		setnames(markers,1:2,c("probe_id","coef"))
		markers = markers[-1]

		i = match(markers$probe_id,rownames(beta))
		y = beta[i,]
		y = y * markers$coef
		y = colSums(y,na.rm=TRUE)
		y = y + 41.72579759
		return(y)

	}else if(model=="placenta:gestational_age"){

		### https://doi.org/10.18632/aging.102049
		markers = system.file(paste0("data/102049-SupFile1.csv"),package="ewastools")
		markers = fread(markers)
		setnames(markers,1:2,c("probe_id","coef"))
		markers = markers[-1]

		i = match(markers$probe_id,rownames(beta))
		y = beta[i,]
		y = y * markers$coef
		y = colSums(y,na.rm=TRUE)
		y = y + 24.99772133
		return(y)

	}else if(model=="horvath_clock"){

		### model based on https://doi.org/10.1186/gb-2013-14-10-r115
		markers = system.file(paste0("data/horvath_clock.csv"),package="ewastools")
		markers = fread(markers,header=TRUE)
		i = match(markers$probe_id,rownames(beta))
		y = beta[i,]
		y = y * markers$coef
		y = colSums(y,na.rm=TRUE)
		y = y + 0.695507258

		inverse_transformation <- function(x,adult.age=20) { ifelse(x<0, (1+adult.age)*exp(x)-1, (1+adult.age)*x+adult.age) }

		return(inverse_transformation(y))

	}else return(NULL)

}
