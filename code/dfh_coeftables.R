##
#
# DenverFINRISKHacky
# Construct Excel-friendly tables for model coefficients from a 'coefs.RData' output
# To be provided as a supplementary table alongside the publication
#
##

# For writing to Excel
library(writexl)

exportcoefs <- function(
	coefs,
	outfile = "coefs.xlsx"
){
	modules <- do.call("cbind", lapply(coefs, FUN=\(x){ x$coefs_modules }))
	colnames(modules) <- paste0("seed_", 1:10)
	modules <- cbind(Variable_Name = rownames(modules), modules)

	module_agesex <- do.call("cbind", lapply(coefs, FUN=\(x){ x$coefs_module_agesex }))
	colnames(module_agesex) <- paste0("seed_", 1:10)
	module_agesex <- cbind(Variable_Name = rownames(module_agesex), module_agesex)

	module_metamix <- do.call("cbind", lapply(coefs, FUN=\(x){ x$coefs_module_metamix }))
	colnames(module_metamix) <- paste0("seed_", 1:10)
	module_metamix <- cbind(Variable_Name = rownames(module_metamix), module_metamix)

	module_alpha <- do.call("cbind", lapply(coefs, FUN=\(x){ x$coefs_module_alpha }))
	colnames(module_alpha) <- paste0("seed_", 1:10)
	module_alpha <- cbind(Variable_Name = rownames(module_alpha), module_alpha)

	module_beta <- do.call("cbind", lapply(coefs, FUN=\(x){ x$coefs_module_beta }))
	colnames(module_beta) <- paste0("seed_", 1:10)
	module_beta <- cbind(Variable_Name = rownames(module_beta), module_beta)
	
	module_relabus <- do.call("cbind", lapply(coefs, FUN=\(x){ x$coefs_module_relabus }))
	colnames(module_relabus) <- paste0("seed_", 1:10)
	module_relabus <- cbind(Variable_Name = rownames(module_relabus), module_relabus)

	module_curated1 <- do.call("cbind", lapply(coefs, FUN=\(x){ x$coefs_module_curated1 }))
	colnames(module_curated1) <- paste0("seed_", 1:10)
	module_curated1 <- cbind(Variable_Name = rownames(module_curated1), module_curated1)

	module_curated2 <- do.call("cbind", lapply(coefs, FUN=\(x){ x$coefs_module_curated2 }))
	colnames(module_curated2) <- paste0("seed_", 1:10)
	module_curated2 <- cbind(Variable_Name = rownames(module_curated2), module_curated2)

	module_trinaryfamilyphylum <- do.call("cbind", lapply(coefs, FUN=\(x){ x$coefs_module_trinaryfamilyphylum }))
	colnames(module_trinaryfamilyphylum) <- paste0("seed_", 1:10)
	module_trinaryfamilyphylum <- cbind(Variable_Name = rownames(module_trinaryfamilyphylum), module_trinaryfamilyphylum)


	dfs <- list(
		"A. Modules" = as.data.frame(modules), 
		"B. module_agesex" = as.data.frame(module_agesex),
		"C. module_metamix" = as.data.frame(module_metamix),
		"D. module_alpha" = as.data.frame(module_alpha),
		"E. module_beta" = as.data.frame(module_beta),
		"F. module_relabus" = as.data.frame(module_relabus),
		"G. module_curated1" = as.data.frame(module_curated1),
		"H. module_curated2" = as.data.frame(module_curated2),
		"I. module_trinaryfamilyphylum" = as.data.frame(module_trinaryfamilyphylum)
	)

	# Write out the list of data frames to an xlsx file	
	writexl::write_xlsx(x = dfs,
	path = outfile)
}

# To be exported as a supplementary table
exportcoefs(coefs = coefs, outfile = "coefs.xlsx")

