###
#
# DenverFINRISKHacky (dfh)
# Post-challenge phase
# Analyze importance of variables, and extract risk contours for key variables
#
###

###
#
# Helper functions
#
###

# Scale risk scores between [0,1]
outscale <- \(x){ (x - min(x, na.rm=TRUE))/(max(x, na.rm=TRUE) - min(x, na.rm=TRUE)) }
# Just left shift variables to start at zero
shift <- \(x){ (x-min(x, na.rm=TRUE)) }
# Scale variables with z transformation
zscale <- \(x){
	if(sd(x, na.rm=TRUE) == 0){
		rep(0, times=length(x)) # Singular vector, returning unit values
	}else{
		(x-mean(x, na.rm=TRUE))/sd(x, na.rm=TRUE) # sd != 0, can apply z-score normalization
	}
}
# Left shifted and squared
sqshift <- \(x){ (x-min(x, na.rm=TRUE))^2 }
# Left shifted and n*log(n+1) multiplied
nlognshift <- \(x){ x <- x-min(x, na.rm=TRUE); x*log(x+1) }
## Take all pairwise interactions and create derived variables, and pairwise interactions from two data matrices
#
# From old code: https://github.com/Syksy/ePCR/blob/master/R/int.R
interact.all <- function(input){
	output <- do.call("cbind", lapply(1:ncol(input), FUN=function(z){ 
		do.call("cbind", lapply(z:ncol(input), FUN=function(x){
			tmp <- data.frame(input[,z] * input[,x])
			colnames(tmp)[1] <- paste(colnames(input)[z], "x", colnames(input)[x], sep="") 
			tmp
		}))
	}))
	output
}
interact.part <- function(input, first, second){
	output <- do.call("cbind", lapply(first, FUN=function(z){ 
		do.call("cbind", lapply(second, FUN=function(x){
			tmp <- data.frame(input[,z] * input[,x])
			colnames(tmp)[1] <- paste(z, "x", x, sep="") 
			tmp
		}))
	}))
	output
}
#' Impute median values to missing slots based on columns
imputationNaive <- function(
	x, # Data matrix 
	FUN = \(q){ q[is.na(q)] <- median(q, na.rm=TRUE); q },
	MARGIN = 2,
	...
){
	dimn <- dimnames(x)
	x <- apply(x, MARGIN=MARGIN, FUN=FUN)
	dimnames(x) <- dimn
	x
}
# Create a self-created, relatively arbitrary "disease burden proxy" variable out of multiple risk factors combined
diseaseburden <- \(x){
	burden <- 0
	# Elevated risk for men at younger age; step at 50 years for men, 65 years for women
	burden <- burden + ifelse(x["Sex"] == 1, ifelse(x["Age"] >= 50, 1, 0), ifelse(x["Age"] >= 65, 1, 0))
	# Overweight
	burden <- burden + ifelse(x["BodyMassIndex"] >= 25, 1, 0)
	# Obese adds a second point
	burden <- burden + ifelse(x["BodyMassIndex"] >= 30, 1, 0)
	# Smoking increases overall disease burden
	burden <- burden + ifelse(x["Smoking"] == 1, 1, 0)
	# Diabetes increases risk
	burden <- burden + ifelse(x["PrevalentDiabetes"] == 1, 1, 0)
	# Prevalent conditions increase risk
	burden <- burden + ifelse(x["PrevalentCHD"] == 1, 1, 0)
	burden <- burden + ifelse(x["PrevalentHFAIL"] == 1, 1, 0)
	# Hypertension increases risk of cardiovascular events; setting threshod at 140; could be stratified according to gender as well
	burden <- burden + ifelse(x["SystolicBP"] >= 140, 1, 0)
	# Taking non-HDL cholesterol risk level of 3.37 threshold from https://www.mayoclinic.org/diseases-conditions/high-blood-cholesterol/expert-answers/cholesterol-ratio/faq-20058006
	burden <- burden + ifelse(x["NonHDLcholesterol"] >= 3.37, 1, 0)		
}

## Create data.frames / matrices similar to which were produced by the final submitted model

# clin = "module_agesex"
# clin2 = "module_metamix"

df_clins <- function(
	clin # Raw input of clinical variables
){
	# Naive median imputation to phenodata to remove NAs
	clin <- as.data.frame(imputationNaive(clin))
	
	clin[,"DiseaseBurden"] <- apply(clin, MARGIN=1, FUN=diseaseburden)
	clin[,"AgeShift2"] <- sqshift(clin[,"Age"])
	clin[,"AgeNlogNshift"] <- nlognshift(clin[,"Age"])
	
	# Sticking just to clinical variables, Firmicute to Bacteroidetes ratio as redundant
	clin[,"F2P"] <- 0
	
	# Categorized variables
	clin[,"BMICat"] <- findInterval(clin$BodyMassIndex, c(-Inf, 25, 30, Inf))
	clin[,"AgeCat"] <- findInterval(clin$Age, c(-Inf, 50, 65, 80, Inf))
	clin[,"SysBPCat"] <- findInterval(clin$SystolicBP, c(-Inf, 130, 150, 180, Inf))
	clin[,"NormalMale"] <- as.numeric(clin$BodyMassIndex > 18 & clin$BodyMassIndex < 25 & clin$Sex == 1)
	clin[,"OverweightMale"] <- as.numeric(clin$BodyMassIndex >= 25 & clin$BodyMassIndex < 30 & clin$Sex == 1)
	clin[,"ObeseMale"] <- as.numeric(clin$BodyMassIndex >= 30 & clin$Sex == 1)
	clin[,"NormalFemale"] <- as.numeric(clin$BodyMassIndex > 18 & clin$BodyMassIndex < 25 & clin$Sex == 0)
	clin[,"OverweightFemale"] <- as.numeric(clin$BodyMassIndex >= 25 & clin$BodyMassIndex < 30 & clin$Sex == 0)
	clin[,"ObeseFemale"] <- as.numeric(clin$BodyMassIndex >= 30 & clin$Sex == 0)

	# Create shifted & z-scored clinical variables, with all pairwise interactions incorporated
	clin2a <- apply(clin[,which(!colnames(clin) %in% c("Event", "Event_time"))], MARGIN=2, FUN=shift)
	clin2b <- apply(clin[,which(!colnames(clin) %in% c("Event", "Event_time"))], MARGIN=2, FUN=zscale)
	colnames(clin2a) <- paste0("s_", colnames(clin2a))
	colnames(clin2b) <- paste0("z_", colnames(clin2b))
	clin2 <- cbind(clin2a, clin2b, interact.all(clin2a), interact.all(clin2b))

	# Omit Event_time and Event from raw clinical data, then obtain combinations
	clin <- clin[,which(!colnames(clin) %in% c("Event", "Event_time"))]

	clin <- cbind(clin, interact.all(clin))
	
	# Return the calculated input data for 'module_agesex' and 'module_metamix'
	list(clin, clin2)
}

# Create a grid in which to test model coefficient sums
clin_df <- expand.grid(
	Age = 20:80,
	BodyMassIndex = 15:35,
	Smoking = 0,
	BPTreatment = 0,
	PrevalentDiabetes = 0,
	PrevalentCHD = 0,
	PrevalentHFAIL = 0,
	Event = NA_integer_,
	Event_time = NA_real_,
	SystolicBP = 100,
	NonHDLcholesterol = 7,
	Sex = 0:1
)

clin_df2 <- df_clins(clin = clin_df)

# Aggregate risks over seeds while taking into account which redundant variables were dropped
aggregate <- function(
	clins,
	coefs
){
	# Extract risks across seeds
	seeds <- lapply(1:length(coefs), FUN=\(i){
		coefs1 <- coefs[[i]]$coefs_module_agesex
		coefs2 <- coefs[[i]]$coefs_module_metamix
		
		#clins1 <- clins[[1]][which(colnames(clins[[1]]) %in% names(coefs1))]
		#clins2 <- clins[[2]][which(colnames(clins[[2]]) %in% names(coefs2))]
		
		clins1 <- clins[[1]][,which(colnames(clins[[1]]) %in% names(coefs1))]
		clins2 <- clins[[2]][,which(colnames(clins[[2]]) %in% names(coefs2))]
		
		print(dim(as.matrix(coefs1)))
		print(dim(as.matrix(coefs2)))
		print(dim(as.matrix(clins1)))
		print(dim(as.matrix(clins2)))
		
		# module_agesex-effect
		t(as.matrix(coefs1)) %*% t(as.matrix(clins1)) +
			# module_metamix-effect
			t(as.matrix(coefs2)) %*% t(as.matrix(clins2))
		
	})
	# Predictions were averaged across seeds; take means	
	colMeans(do.call("rbind", seeds))
}

# coefs is loaded from the coefs.RData of interest
tmp <- aggregate(clins = clin_df2, coefs = coefs)
 
contour_m_agebmi <- matrix(NA, nrow=length(unique(clin_df$BodyMassIndex)), ncol=length(unique(clin_df$Age)))
contour_f_agebmi <- matrix(NA, nrow=length(unique(clin_df$BodyMassIndex)), ncol=length(unique(clin_df$Age)))

rownames(contour_m_agebmi) <- rownames(contour_f_agebmi) <- paste0("BMI_", unique(clin_df$BodyMassIndex))
colnames(contour_m_agebmi) <- colnames(contour_f_agebmi) <- paste0("Age_", unique(clin_df$Age))

# Populate the risk score contours for male/female
for(age in 1:length(unique(clin_df$Age))){
	for(bmi in 1:length(unique(clin_df$BodyMassIndex))){
		contour_m_agebmi[bmi, age] <- tmp[which(clin_df$Sex == 1 & clin_df$Age == unique(clin_df$Age)[age] & clin_df$BodyMassIndex == unique(clin_df$BodyMassIndex)[bmi])]
		contour_f_agebmi[bmi, age] <- tmp[which(clin_df$Sex == 0 & clin_df$Age == unique(clin_df$Age)[age] & clin_df$BodyMassIndex == unique(clin_df$BodyMassIndex)[bmi])]
	}
}


# Variant #1, heatmap based plotting

library(hamlet) # My own package available in CRAN

cols <- colorRampPalette(c("orange", "red", "black", "blue", "cyan"), bias=3)(1000)
vals <- seq(from=0, to=2, length.out=1000)

xlabs <- unique(clin_df$Age)
xlabs[-seq(1,length(xlabs),by=5)] <- NA
ylabs <- unique(clin_df$BodyMassIndex)
ylabs[-seq(1,length(ylabs),by=5)] <- NA

par(mfrow=c(1,2), mar=c(4,4,3,1))
h1 <- hmap(contour_m_agebmi, col = cols, valseq = vals, Rowv = NA, Colv = NA, namerows = ylabs, namecols = xlabs)
title(xlab = "\n\nAge", ylab = "\n\nBMI", main = "\n\nMale")
h2 <- hmap(contour_f_agebmi, col = cols, valseq = vals, Rowv = NA, Colv = NA, namerows = ylabs, namecols = xlabs)
title(xlab = "\n\nAge", ylab = "\n\nBMI", main = "\n\nFemale")
hmap.key(h2)

# Variant #2, contour lines in ggplot2

library(ggplot2)

clin_df[,"Risk"] <- tmp

p <- ggplot(clin_df, aes(x = Age, y = BodyMassIndex, z = Risk)) +
	stat_contour(geom = "polygon", aes(fill = ..level..)) +
	geom_tile(aes(fill = Risk)) +
	stat_contour(bins = 15) +
	xlab("Age (Years)") +
	ylab("BMI") +
	guides(fill = guide_colorbar(title = "Predicted risk"))

p + facet_grid(rows = vars(Sex), labeller = labeller(Sex = c("0" = "Female", "1" = "Male"))



