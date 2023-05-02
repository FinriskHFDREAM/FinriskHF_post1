##
#
# Downstream refinement of model(s) submitted to FINRISK DREAM HF Challenge
#
##

#' Reads coefficients' output at 'coefs.RData' which is written to output folder and aggregates this to a human-readable tab-separated table
readCoefsRData <- function(
	path # Path to coefs.RData; e.g. ../DenverFINRISKHacky_Post1/output/ 
){
	# Load the variable 'coefs'
	load(file.path(path, "coefs.RData"))
	
	# Summarize the modules selected across the CV seeds
	modules <- do.call("cbind", lapply(coefs, FUN=\(x) { x$coefs_modules }))
	
	# TODO: Of the modules that performed well, extract information on which variables contributed there-in
	
	# Return modules identified as predictive per each seed
	modules
}
