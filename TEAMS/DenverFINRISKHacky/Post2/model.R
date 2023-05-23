# DenverFINRISKHacky 
# Post-challenge phase LASSO-testing along with model coefficient output
# Preparing to use the modelling process after the Challenge, extract most important variables, etc

###
#
# Load libraries
#
###
library(survival)
library(glmnet)
library(microbiome)
library(phyloseq)
library(mia)
library(TreeSummarizedExperiment)
library(ecodist)
library(vegan)

# List of coefficients across multiple seeds 
coefs <- list()

###
#
# Define main model function
#
###
model <- function(
	# Clinical data matrix for both training and testing cohort
	train_clin,
	test_clin,
	# Microbiome (raw) read counts for both training and testing cohort
	train_micr,
	test_micr,
	# phyloseq-objects for the train and test cohorts
	train_phyloseq,
	test_phyloseq,
	# Which model submission version to use 
	v,
	# Alpha 'a', which controls the ratio of L1 (LASSO, a=1) to L2 (RR, a=0) regularization
	a = 1,
	# Additional parameters
	# Run across multiple seeds
	seeds = 1:3, # RNG seed for reproducibility
	...
){

	# Small system time print function to help track runtimes; checking that the runtimes stay reasonable in a small-scale Ubuntu 22.04 LTS VM
	catsystime <- \(x){
		cat("\n\nCurrent system time:\n")
		cat(as.character(Sys.time()))
		cat("\n")
		if(!missing(x)) cat(paste("Current step:", x, "\n"))
		cat("\n\n")
	}
	
	catsystime("Start")
	start_time <- Sys.time()

	# Start constructing the output df which will be output as csv
	output_temp <- data.frame(SampleID = rownames(test_clin))
	output_final <- data.frame(SampleID = rownames(test_clin), Score = 0)

	# Discard negative or zero event times
	discard <- which(train_clin$Event_time <= 0)
	train_clin <- train_clin[-discard,]
	otu_table(train_phyloseq) <- otu_table(train_phyloseq)[,-discard]
	train_micr <- train_micr[,-discard]

	# Construct an ensemble training df
	ensemble_temp <- data.frame(SampleID = rownames(train_clin), Event = train_clin[,"Event"], Event_time = train_clin[,"Event_time"])

	# Naive median imputation to phenodata to remove NAs
	train_clin <- as.data.frame(imputationNaive(train_clin))
	test_clin <- as.data.frame(imputationNaive(test_clin))

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
	train_clin[,"DiseaseBurden"] <- apply(train_clin, MARGIN=1, FUN=diseaseburden)
	test_clin[,"DiseaseBurden"] <- apply(test_clin, MARGIN=1, FUN=diseaseburden)
	
	# Squared & shifted age
	train_clin[,"AgeShift2"] <- sqshift(train_clin[,"Age"])
	train_clin[,"AgeNlogNshift"] <- nlognshift(train_clin[,"Age"])

	test_clin[,"AgeShift2"] <- sqshift(test_clin[,"Age"])
	test_clin[,"AgeNlogNshift"] <- nlognshift(test_clin[,"Age"])

	# Basic clinical variables comparable to the baseline model; Age, Sex, and BMI unmodified as numerics
	train_baseclin <- train_clin[,c("Age", "Sex", "BodyMassIndex")]
	test_baseclin <- test_clin[,c("Age", "Sex", "BodyMassIndex")]

	# Literature or other source curated & weighted information on possibly interesting microbiome markers / phenodata
	catsystime("Creating abundance & curated module data...")
	
	# Species
	train_abuspecies <- microbiome::abundances(microbiome::aggregate_taxa(train_phyloseq, level = "Species"), transform="compositional")
	test_abuspecies <- microbiome::abundances(microbiome::aggregate_taxa(test_phyloseq, level = "Species"), transform="compositional")
	# Genus
	train_abugenus <- microbiome::abundances(microbiome::aggregate_taxa(train_phyloseq, level = "Genus"), transform="compositional")
	test_abugenus <- microbiome::abundances(microbiome::aggregate_taxa(test_phyloseq, level = "Genus"), transform="compositional")
	# Family
	train_abufamily <- microbiome::abundances(microbiome::aggregate_taxa(train_phyloseq, level = "Family"), transform="compositional")
	test_abufamily <- microbiome::abundances(microbiome::aggregate_taxa(test_phyloseq, level = "Family"), transform="compositional")
	# Order
	train_abuorder <- microbiome::abundances(microbiome::aggregate_taxa(train_phyloseq, level = "Order"), transform="compositional")
	test_abuorder <- microbiome::abundances(microbiome::aggregate_taxa(test_phyloseq, level = "Order"), transform="compositional")
	# Class
	train_abuclass <- microbiome::abundances(microbiome::aggregate_taxa(train_phyloseq, level = "Class"), transform="compositional")
	test_abuclass <- microbiome::abundances(microbiome::aggregate_taxa(test_phyloseq, level = "Class"), transform="compositional")
	# Phylum
	train_abuphylum <- microbiome::abundances(microbiome::aggregate_taxa(train_phyloseq, level = "Phylum"), transform="compositional")
	test_abuphylum <- microbiome::abundances(microbiome::aggregate_taxa(test_phyloseq, level = "Phylum"), transform="compositional")
	# Domain/Kingdom
	train_abudomain <- microbiome::abundances(microbiome::aggregate_taxa(train_phyloseq, level = "Domain"), transform="compositional")
	test_abudomain <- microbiome::abundances(microbiome::aggregate_taxa(test_phyloseq, level = "Domain"), transform="compositional")
	
	# Certain taxa of interest from various levels taken from literature
	train_curated1 <- data.frame(
		p__Firmicutes = apply(train_abuphylum[grep("p__Firmicutes", rownames(train_abuphylum), value=TRUE),,drop=FALSE], MARGIN=2, FUN=sum),
		p__Bacteroidetes = apply(train_abuphylum[grep("p__Bacteroidetes", rownames(train_abuphylum), value=TRUE),,drop=FALSE], MARGIN=2, FUN=sum),
		g__Lactobacillus = apply(train_abugenus[grep("g__Lactobacillus", rownames(train_abugenus), value=TRUE),,drop=FALSE], MARGIN=2, FUN=sum),
		g__Roseburia = apply(train_abugenus[grep("g__Roseburia", rownames(train_abugenus), value=TRUE),,drop=FALSE], MARGIN=2, FUN=sum),
		s__Escherichia_coli = apply(train_abuspecies[grep("Escherichia_coli", rownames(train_abuspecies), value=TRUE),,drop=FALSE], MARGIN=2, FUN=sum),
		s__Klebsiella_pneumonia = apply(train_abuspecies[grep("Klebsiella_pneumonia", rownames(train_abuspecies), value=TRUE),,drop=FALSE], MARGIN=2, FUN=sum),
		p__Proteobacteria = apply(train_abuphylum[grep("p__Proteobacteria", rownames(train_abuphylum), value=TRUE),,drop=FALSE], MARGIN=2, FUN=sum),
		p__Actinobacteria = apply(train_abuphylum[grep("p__Actinobacteria", rownames(train_abuphylum), value=TRUE),,drop=FALSE], MARGIN=2, FUN=sum),
		f__Enterobacteriaceae = apply(train_abufamily[grep("f__Enterobacteriaceae", rownames(train_abufamily), value=TRUE),,drop=FALSE], MARGIN=2, FUN=sum),
		g__Streptococcus = apply(train_abugenus[grep("g__Streptococcus", rownames(train_abugenus), value=TRUE),,drop=FALSE], MARGIN=2, FUN=sum),
		s__Roseburia_intestinalis = apply(train_abuspecies[grep("s__Roseburia_intestinalis", rownames(train_abuspecies), value=TRUE),,drop=FALSE], MARGIN=2, FUN=sum),
		s__Faecalibacterium_prausnitzii = apply(train_abuspecies[grep("s__Faecalibacterium_prausnitzii", rownames(train_abuspecies), value=TRUE),,drop=FALSE], MARGIN=2, FUN=sum),
		g__Enterococcus = apply(train_abugenus[grep("g__Enterococcus", rownames(train_abugenus), value=TRUE),,drop=FALSE], MARGIN=2, FUN=sum),
		g__Subdoligranulum = apply(train_abugenus[grep("g__Subdoligranulum", rownames(train_abugenus), value=TRUE),,drop=FALSE], MARGIN=2, FUN=sum),
		f__Coriobacteriaceae = apply(train_abufamily[grep("f__Coriobacteriaceae", rownames(train_abufamily), value=TRUE),,drop=FALSE], MARGIN=2, FUN=sum),
		f__Erysipelotrichaceae = apply(train_abufamily[grep("f__Erysipelotrichaceae", rownames(train_abufamily), value=TRUE),,drop=FALSE], MARGIN=2, FUN=sum),
		f__Ruminococcaceae = apply(train_abufamily[grep("f__Ruminococcaceae", rownames(train_abufamily), value=TRUE),,drop=FALSE], MARGIN=2, FUN=sum),
		g__Blautia = apply(train_abugenus[grep("g__Blautia", rownames(train_abugenus), value=TRUE),,drop=FALSE], MARGIN=2, FUN=sum)
	)
	# Add interactions
	train_curated1 <- cbind(train_curated1, interact.all(train_curated1))
		
	# Loosely generated literature-driven abundance combinations		
	train_curated2 <- data.frame(
		F2P = train_curated1$p__Firmicutes / ( train_curated1$p__Bacteroidetes + 0.0001), # Add small systematic epsilon to avoid division by zero
		Atherosclerosis = train_curated1$g__Lactobacillus - train_curated1$g__Roseburia,
		HF = train_curated1$s__Escherichia_coli - train_curated1$s__Klebsiella_pneumonia,
		CKD = train_curated1$p__Firmicutes + train_curated1$p__Proteobacteria + train_curated1$p__Actinobacteria,
		CAD = train_curated1$p__Firmicutes - train_curated1$p__Bacteroidetes + train_curated1$f__Enterobacteriaceae 
			+ train_curated1$g__Streptococcus - train_curated1$s__Roseburia_intestinalis - train_curated1$s__Faecalibacterium_prausnitzii,
		HFP = train_curated1$f__Ruminococcaceae - train_curated1$f__Coriobacteriaceae - train_curated1$f__Erysipelotrichaceae
	)
	# Add interactions
	train_curated2 <- cbind(train_curated2, interact.all(train_curated2))
			
	# Certain taxa of interest from various levels taken from literature
	test_curated1 <- data.frame(
		p__Firmicutes = apply(test_abuphylum[grep("p__Firmicutes", rownames(test_abuphylum), value=TRUE),,drop=FALSE], MARGIN=2, FUN=sum),
		p__Bacteroidetes = apply(test_abuphylum[grep("p__Bacteroidetes", rownames(test_abuphylum), value=TRUE),,drop=FALSE], MARGIN=2, FUN=sum),
		g__Lactobacillus = apply(test_abugenus[grep("g__Lactobacillus", rownames(test_abugenus), value=TRUE),,drop=FALSE], MARGIN=2, FUN=sum),
		g__Roseburia = apply(test_abugenus[grep("g__Roseburia", rownames(test_abugenus), value=TRUE),,drop=FALSE], MARGIN=2, FUN=sum),
		s__Escherichia_coli = apply(test_abuspecies[grep("Escherichia_coli", rownames(test_abuspecies), value=TRUE),,drop=FALSE], MARGIN=2, FUN=sum),
		s__Klebsiella_pneumonia = apply(test_abuspecies[grep("Klebsiella_pneumonia", rownames(test_abuspecies), value=TRUE),,drop=FALSE], MARGIN=2, FUN=sum),
		p__Proteobacteria = apply(test_abuphylum[grep("p__Proteobacteria", rownames(test_abuphylum), value=TRUE),,drop=FALSE], MARGIN=2, FUN=sum),
		p__Actinobacteria = apply(test_abuphylum[grep("p__Actinobacteria", rownames(test_abuphylum), value=TRUE),,drop=FALSE], MARGIN=2, FUN=sum),
		f__Enterobacteriaceae = apply(test_abufamily[grep("f__Enterobacteriaceae", rownames(test_abufamily), value=TRUE),,drop=FALSE], MARGIN=2, FUN=sum),
		g__Streptococcus = apply(test_abugenus[grep("g__Streptococcus", rownames(test_abugenus), value=TRUE),,drop=FALSE], MARGIN=2, FUN=sum),
		s__Roseburia_intestinalis = apply(test_abuspecies[grep("s__Roseburia_intestinalis", rownames(test_abuspecies), value=TRUE),,drop=FALSE], MARGIN=2, FUN=sum),
		s__Faecalibacterium_prausnitzii = apply(test_abuspecies[grep("s__Faecalibacterium_prausnitzii", rownames(test_abuspecies), value=TRUE),,drop=FALSE], MARGIN=2, FUN=sum),
		g__Enterococcus = apply(test_abugenus[grep("g__Enterococcus", rownames(test_abugenus), value=TRUE),,drop=FALSE], MARGIN=2, FUN=sum),
		g__Subdoligranulum = apply(test_abugenus[grep("g__Subdoligranulum", rownames(test_abugenus), value=TRUE),,drop=FALSE], MARGIN=2, FUN=sum),
		f__Coriobacteriaceae = apply(test_abufamily[grep("f__Coriobacteriaceae", rownames(test_abufamily), value=TRUE),,drop=FALSE], MARGIN=2, FUN=sum),
		f__Erysipelotrichaceae = apply(test_abufamily[grep("f__Erysipelotrichaceae", rownames(test_abufamily), value=TRUE),,drop=FALSE], MARGIN=2, FUN=sum),
		f__Ruminococcaceae = apply(test_abufamily[grep("f__Ruminococcaceae", rownames(test_abufamily), value=TRUE),,drop=FALSE], MARGIN=2, FUN=sum),
		g__Blautia = apply(test_abugenus[grep("g__Blautia", rownames(test_abugenus), value=TRUE),,drop=FALSE], MARGIN=2, FUN=sum)
	)
	# Add interactions
	test_curated1 <- cbind(test_curated1, interact.all(test_curated1))

	# Loosely generated literature-driven abundance combinations		
	test_curated2 <- data.frame(
		F2P = test_curated1$p__Firmicutes / ( test_curated1$p__Bacteroidetes + 0.0001), # Add small systematic epsilon to avoid division by zero
		Atherosclerosis = test_curated1$g__Lactobacillus - test_curated1$g__Roseburia,
		HF = test_curated1$s__Escherichia_coli - test_curated1$s__Klebsiella_pneumonia,
		CKD = test_curated1$p__Firmicutes + test_curated1$p__Proteobacteria + test_curated1$p__Actinobacteria,
		CAD = test_curated1$p__Firmicutes - test_curated1$p__Bacteroidetes + test_curated1$f__Enterobacteriaceae 
			+ test_curated1$g__Streptococcus - test_curated1$s__Roseburia_intestinalis - test_curated1$s__Faecalibacterium_prausnitzii,
		HFP = test_curated1$f__Ruminococcaceae - test_curated1$f__Coriobacteriaceae - test_curated1$f__Erysipelotrichaceae
	)
	# Add interactions
	test_curated2 <- cbind(test_curated2, interact.all(test_curated2))
			
	# Add Firmicutes to Bacteroidetes to clinical phenodata due to potential interesting interactions with phenodata
	train_clin[,"F2P"] <- train_curated1[,"p__Firmicutes"] / ( train_curated1[,"p__Bacteroidetes"] + 0.0001) # Add small systematic epsilon to avoid division by zero
	test_clin[,"F2P"] <- test_curated1[,"p__Firmicutes"] / ( test_curated1[,"p__Bacteroidetes"] + 0.0001) # Add small systematic epsilon to avoid division by zero

	# Add categorized versions of the few select continuous phenodata
	# BMI categorized
	train_clin[,"BMICat"] <- findInterval(train_clin$BodyMassIndex, c(-Inf, 25, 30, Inf))
	test_clin[,"BMICat"] <- findInterval(test_clin$BodyMassIndex, c(-Inf, 25, 30, Inf))
	# Age categorized
	train_clin[,"AgeCat"] <- findInterval(train_clin$Age, c(-Inf, 50, 65, 80, Inf))
	test_clin[,"AgeCat"] <- findInterval(test_clin$Age, c(-Inf, 50, 65, 80, Inf))
	# Systolic Blood Pressure categorized
	train_clin[,"SysBPCat"] <- findInterval(train_clin$SystolicBP, c(-Inf, 130, 150, 180, Inf))
	test_clin[,"SysBPCat"] <- findInterval(test_clin$SystolicBP, c(-Inf, 130, 150, 180, Inf))
	# Binary identifiers (particularly for interaction effects that may interact with sex and BMI 
	# Normal weight Male specific identifier
	train_clin[,"NormalMale"] <- as.numeric(train_clin$BodyMassIndex > 18 & train_clin$BodyMassIndex < 25 & train_clin$Sex == 1)
	test_clin[,"NormalMale"] <- as.numeric(test_clin$BodyMassIndex > 18 & test_clin$BodyMassIndex < 25 & test_clin$Sex == 1)
	# Over-weight Male specific identifier
	train_clin[,"OverweightMale"] <- as.numeric(train_clin$BodyMassIndex >= 25 & train_clin$BodyMassIndex < 30 & train_clin$Sex == 1)
	test_clin[,"OverweightMale"] <- as.numeric(test_clin$BodyMassIndex >= 25 & test_clin$BodyMassIndex < 30 & test_clin$Sex == 1)
	# Obese Male specific identifier
	train_clin[,"ObeseMale"] <- as.numeric(train_clin$BodyMassIndex >= 30 & train_clin$Sex == 1)
	test_clin[,"ObeseMale"] <- as.numeric(test_clin$BodyMassIndex >= 30 & test_clin$Sex == 1)
	# Normal weight Female specific identifier
	train_clin[,"NormalFemale"] <- as.numeric(train_clin$BodyMassIndex > 18 & train_clin$BodyMassIndex < 25 & train_clin$Sex == 0)
	test_clin[,"NormalFemale"] <- as.numeric(test_clin$BodyMassIndex > 18 & test_clin$BodyMassIndex < 25 & test_clin$Sex == 0)
	# Over-weight Female specific identifier
	train_clin[,"OverweightFemale"] <- as.numeric(train_clin$BodyMassIndex >= 25 & train_clin$BodyMassIndex < 30 & train_clin$Sex == 0)
	test_clin[,"OverweightFemale"] <- as.numeric(test_clin$BodyMassIndex >= 25 & test_clin$BodyMassIndex < 30 & test_clin$Sex == 0)
	# Obese Female specific identifier
	train_clin[,"ObeseFemale"] <- as.numeric(train_clin$BodyMassIndex >= 30 & train_clin$Sex == 0)
	test_clin[,"ObeseFemale"] <- as.numeric(test_clin$BodyMassIndex >= 30 & test_clin$Sex == 0)

	# Create shifted & z-scored clinical variables, with all pairwise interactions incorporated
	train_clin2a <- apply(train_clin[,which(!colnames(train_clin) %in% c("Event", "Event_time"))], MARGIN=2, FUN=shift)
	train_clin2b <- apply(train_clin[,which(!colnames(train_clin) %in% c("Event", "Event_time"))], MARGIN=2, FUN=zscale)
	colnames(train_clin2a) <- paste0("s_", colnames(train_clin2a))
	colnames(train_clin2b) <- paste0("z_", colnames(train_clin2b))
	train_clin2 <- cbind(train_clin2a, train_clin2b, interact.all(train_clin2a), interact.all(train_clin2b))

	test_clin2a <- apply(test_clin[,which(!colnames(test_clin) %in% c("Event", "Event_time"))], MARGIN=2, FUN=shift)
	test_clin2b <- apply(test_clin[,which(!colnames(test_clin) %in% c("Event", "Event_time"))], MARGIN=2, FUN=zscale)
        colnames(test_clin2a) <- paste0("s_", colnames(test_clin2a))
        colnames(test_clin2b) <- paste0("z_", colnames(test_clin2b))	
	test_clin2 <- cbind(test_clin2a, test_clin2b, interact.all(test_clin2a), interact.all(test_clin2b))

	# Omit combinations that produce NaN
	train_clin2 <- train_clin2[,which(!apply(train_clin2, MARGIN=2, FUN=\(x){ any(!is.finite(x)) }))]
	# Keep same columns
	test_clin2 <- test_clin2[,colnames(train_clin2)]

	# Construct response Surv and remove it from the pheno data
	train_y <- survival::Surv(time = train_clin$Event_time, event = train_clin$Event)	

	# Omit Event_time and Event from raw clinical data, then obtain combinations
	train_clin <- train_clin[,which(!colnames(train_clin) %in% c("Event", "Event_time"))]
	test_clin <- test_clin[,which(!colnames(test_clin) %in% c("Event", "Event_time"))]
	
	# Interaction terms added
	train_clin <- cbind(train_clin, interact.all(train_clin))
	test_clin <- cbind(test_clin, interact.all(test_clin))

	# Microbiome diversity metrics

	catsystime("Alpha diversity metrics")
	# Alpha
	catsystime("Training data...")
	train_alpha <- train_phyloseq |>
		microbiome::alpha() |>
		imputationNaive()

	catsystime("Test data...")
	test_alpha <- test_phyloseq |>
		microbiome::alpha() |>
		imputationNaive()

	# Introduce diversity/entropy metric interactions
	train_alpha <- cbind(train_alpha, interact.all(train_alpha))
	test_alpha <- cbind(test_alpha, interact.all(test_alpha))

	# Beta, examine PCoA vectors up to 10th
	catsystime("Beta diversity metrics")
	catsystime("Training data...")
	train_tse <- train_phyloseq |>
		mia::makeTreeSummarizedExperimentFromPhyloseq() |>
		mia::agglomerateByRank(x = _, rank = "Genus") |>
		mia::transformSamples(x = _, method = "relabundance") |>
		mia::transformSamples(x = _, abund_values = "relabundance", pseudocount = 1, method = "clr")

	catsystime("vegan & ecodist...")
	train_beta <- t(assays(train_tse)$relabundance) |>
		vegan::vegdist(x = _, method = "bray") |>
		ecodist::pco(x = _) |>
		(\(x) { x <- x$vectors[,1:10]; colnames(x) <- paste0("PCoA", 1:10); x })()

	catsystime("Test data...")
        test_tse <- test_phyloseq |>
                mia::makeTreeSummarizedExperimentFromPhyloseq() |>
                mia::agglomerateByRank(x = _, rank = "Genus") |>
                mia::transformSamples(x = _, method = "relabundance") |>
                mia::transformSamples(x = _, abund_values = "relabundance", pseudocount = 1, method = "clr")

	catsystime("vegan & ecodist...")
        test_beta <- t(assays(test_tse)$relabundance) |>
                vegan::vegdist(x = _, method = "bray") |>
                ecodist::pco(x = _) |>
                (\(x) { x <- x$vectors[,1:10]; colnames(x) <- paste0("PCoA", 1:10); x })()

	# Test mia/microbiome packages' suggested approaches
	catsystime("TreeSummarizedExperiment-objects")
	# Little piping function for relative abundances
	pip <- function(phylo, level){
		phylo |>
                mia::makeTreeSummarizedExperimentFromPhyloseq() |>
                mia::agglomerateByRank(x = _, rank = level) |>
                mia::transformSamples(x = _, method = "relabundance") |>
                mia::transformSamples(x = _, abund_values = "relabundance", pseudocount = 1, method = "clr", name = "clr_transformation") |>
                (\(x) { assay(x, "clr_transformation") })() |>
		(\(x) { x[which(!rownames(x) %in% c("s__", "g__", "f__", "o__", "c__", "p__", "k__", "d__")),] })()
	}
	# Training data relative abundances
	catsystime("Training data relative abundances...")
	train_relabus <- t(rbind(
		pip(train_phyloseq, level = "Species"),
		pip(train_phyloseq, level = "Genus"),
		pip(train_phyloseq, level = "Family"),
		pip(train_phyloseq, level = "Order"),
		pip(train_phyloseq, level = "Class")#,
		#pip(train_phyloseq, level = "Phylum")
	))
	# Test data relative abundances
	catsystime("Test data relative abundances...")
	test_relabus <- t(rbind(
		pip(test_phyloseq, level = "Species"),
		pip(test_phyloseq, level = "Genus"),
		pip(test_phyloseq, level = "Family"),
		pip(test_phyloseq, level = "Order"),
		pip(test_phyloseq, level = "Class")#,
		#pip(test_phyloseq, level = "Phylum")
	))	

	# Tri-binarized categorical vars
	# Senior vs. Junior defined at 60 years cutoff
	# Obese vs. Nonobese defined at BMI >= 30
	# Male vs. Female defined by Sex == 1 for male as per challenge
	# Combine these with microbiome relative abundances for HF prediction
	train_tricomb <- data.frame(
		SeniorObeseMale = as.numeric(train_clin$Age >= 60 & train_clin$BodyMassIndex >= 30 & train_clin$Sex == 1), 
		JuniorObeseMale = as.numeric(train_clin$Age < 60 & train_clin$BodyMassIndex >= 30 & train_clin$Sex == 1),
		SeniorObeseFemale = as.numeric(train_clin$Age >= 60 & train_clin$BodyMassIndex >= 30 & train_clin$Sex == 0),
		JuniorObeseFemale = as.numeric(train_clin$Age < 60 & train_clin$BodyMassIndex >= 30 & train_clin$Sex == 0),
		SeniorNonobeseMale = as.numeric(train_clin$Age >= 60 & train_clin$BodyMassIndex < 30 & train_clin$Sex == 1),
		JuniorNonobeseMale = as.numeric(train_clin$Age < 60 & train_clin$BodyMassIndex < 30 & train_clin$Sex == 1),
		SeniorNonobeseFemale = as.numeric(train_clin$Age >= 60 & train_clin$BodyMassIndex < 30 & train_clin$Sex == 0),
		JuniorNonobeseFemale = as.numeric(train_clin$Age < 60 & train_clin$BodyMassIndex < 30 & train_clin$Sex == 0)
	)
        test_tricomb <- data.frame(
                SeniorObeseMale = as.numeric(test_clin$Age >= 60 & test_clin$BodyMassIndex >= 30 & test_clin$Sex == 1),
                JuniorObeseMale = as.numeric(test_clin$Age < 60 & test_clin$BodyMassIndex >= 30 & test_clin$Sex == 1),
                SeniorObeseFemale = as.numeric(test_clin$Age >= 60 & test_clin$BodyMassIndex >= 30 & test_clin$Sex == 0),
                JuniorObeseFemale = as.numeric(test_clin$Age < 60 & test_clin$BodyMassIndex >= 30 & test_clin$Sex == 0),
                SeniorNonobeseMale = as.numeric(test_clin$Age >= 60 & test_clin$BodyMassIndex < 30 & test_clin$Sex == 1),
                JuniorNonobeseMale = as.numeric(test_clin$Age < 60 & test_clin$BodyMassIndex < 30 & test_clin$Sex == 1),
                SeniorNonobeseFemale = as.numeric(test_clin$Age >= 60 & test_clin$BodyMassIndex < 30 & test_clin$Sex == 0),
                JuniorNonobeseFemale = as.numeric(test_clin$Age < 60 & test_clin$BodyMassIndex < 30 & test_clin$Sex == 0)
        )

	catsystime("Training data trinary stratified Family-taxa relative abundances...")
	
	# Use Family and Phylum level taxa, appeared to be prevalent in literature and possibly strong enough signal with respect to trinarized strata
	train_tri <- cbind(train_tricomb, 
		t(pip(train_phyloseq, level = "Family"))#,
		#t(pip(train_phyloseq, level = "Phylum"))
	)
	
	# Interactions between the trinarized indicators and relative abundances
	train_tri <- interact.part(train_tri, first = colnames(train_tri)[1:ncol(train_tricomb)], second = colnames(train_tri)[(ncol(train_tricomb)+1):ncol(train_tri)])

	catsystime("Test data trinary stratified Family-taxa relative abundances...")

	# Use Family and Phylum level taxa, appeared to be prevalent in literature and possibly strong enough signal with respect to trinarized strata
	test_tri <- cbind(test_tricomb, 
		t(pip(test_phyloseq, level = "Family"))#,
		#t(pip(test_phyloseq, level = "Phylum"))
	)
	
	# Interactions between the trinarized indicators and relative abundances
	test_tri <- interact.part(test_tri, first = colnames(test_tri)[1:ncol(test_tricomb)], second = colnames(test_tri)[(ncol(test_tricomb)+1):ncol(test_tri)])

	## Modular modelling of risk
	# Generic use module training with glmnet 10-fold CV
	# ... or other generic use of L1 (LASSO) with possibility to incorporate L2-norm as well via alpha
	module_glmnet <- function(
		trainx,
                trainy,
		testx,
		...
	){
		# Cast data matrices to a matrix to avoid issues with data.frames
		trainx <- as.matrix(trainx)
		testx <- as.matrix(testx)
	
		# Sanitize column and rownames to capped length spanning from the end and removing special symbols such a ':'
		rownames(trainx) <- sanitize(rownames(trainx), charlimit=50)
		rownames(testx) <- sanitize(rownames(testx), charlimit=50)
		colnames(trainx) <- sanitize(colnames(trainx), charlimit=50)
		colnames(testx) <- sanitize(colnames(testx), charlimit=50)
		
		# Omit training samples with NAs for response
		trainx <- trainx[which(!is.na(trainy)),]
		trainy <- trainy[which(!is.na(trainy))]

		# Assign 0 values to non-finite input (i.e. NA/NaN/+-Inf)
		trainx[!is.finite(trainx)] <- 0
		trainy[!is.finite(trainy)] <- 0

		# Debugging info
		print("Debugging; print first 3 rows and 3 first columns from the processed train and test data prior to running glmnet and cv.glmnet")
		print(trainx[1:3,1:3])
		print(testx[1:3,1:3])

		# Fit, CV, predict
		# sub >=6; changing type.measure to C-index, LASSO alpha == 1 to Elastic Net alpha == 0.5; alpha 'a' given as global parameter in main model function
		fit <- glmnet(x = trainx, y = trainy, family = "cox", alpha = a)
		cv <- cv.glmnet(x = trainx, y = trainy, family = "cox", type.measure = "C", alpha = a)

		vars <- colnames(trainx)[predict(fit, s = cv$lambda.1se, type = "nonzero")[[1]]]

		if(!is.null(vars) & length(vars)>0){
			print("lambda.1se, non-zero coefs:")			
			print(data.frame(
				Variable = colnames(trainx)[predict(fit, s = cv$lambda.1se, type = "nonzero")[[1]]],
				Coef = predict(fit, s = cv$lambda.1se, type = "coefficient")[predict(fit, s = cv$lambda.1se, type = "nonzero")[[1]]]
				)
			)
		}else{
			print("All coefficients shrunk to zero.")
		}

		# Bundle both train and test module predictions into a list
		preds <- list()

		# Return two list elements; first is for training data, second is for test data predictions (
		preds[[1]] <- predict(fit, newx = as.matrix(trainx), s = cv$lambda.1se, type = "response")[,1]
		preds[[2]] <- predict(fit, newx = as.matrix(testx), s = cv$lambda.1se, type = "response")[,1]
		# Post-challenge phase, also return the used coefficients
		preds[[3]] <- predict(fit, newx = as.matrix(testx), s = cv$lambda.1se, type = "coefficients")[,1]

		# Debugging; in case any of the predictions are non-finite, replace them with the median
		preds[[1]][!is.finite(preds[[1]])] <- median(preds[[1]], na.rm=TRUE)
		preds[[2]][!is.finite(preds[[2]])] <- median(preds[[2]], na.rm=TRUE)

		preds
	}

	# Run across multiple RNG seeds to alleviate random binning effects
	scores <- do.call("cbind", lapply(seeds, FUN=\(s){
		catsystime(paste("\n\nRunning model CVs, with seed", s, "and alpha", a, "\n\n"))	

		# Part I
		set.seed(s)
		catsystime("\nmodule_baseclin\n")
		module_baseclin <- module_glmnet(trainx = train_baseclin, trainy = train_y, testx = test_baseclin)
		catsystime("\nmodule_agesex\n")
		module_agesex <- module_glmnet(trainx = train_clin, trainy = train_y, testx = test_clin)
		catsystime("\nmodule_metamix\n")
		module_metamix <- module_glmnet(trainx = train_clin2, trainy = train_y, testx = test_clin2)	
		catsystime("\nmodule_alpha\n")
		module_alpha <- module_glmnet(trainx = train_alpha, trainy = train_y, testx = test_alpha)
		catsystime("\nmodule_beta\n")
		module_beta <- module_glmnet(trainx = train_beta, trainy = train_y, testx = test_beta)
		catsystime("\nmodule_relabus\n")
		module_relabus <- module_glmnet(trainx = train_relabus, trainy = train_y, testx = test_relabus)
		catsystime("\nmodule_curated1\n")
		module_curated1 <- module_glmnet(trainx = train_curated1, trainy = train_y, testx = test_curated1)
		catsystime("\nmodule_curated2\n")
		module_curated2 <- module_glmnet(trainx = train_curated2, trainy = train_y, testx = test_curated2)
		catsystime("\nmodule_trinaryfamily\n")
		module_trinaryfamily <- module_glmnet(trainx = train_tri, trainy = train_y, testx = test_tri)
		catsystime("\nmodule_species\n")
		module_species <- module_glmnet(trainx = t(train_abuspecies), trainy = train_y, testx = t(test_abuspecies))

		catsystime("Pt Ia")	
		ensemble_temp[,"module_baseclin"] <- module_baseclin[[1]]
		ensemble_temp[,"module_agesex"] <- module_agesex[[1]]
		ensemble_temp[,"module_metamix"] <- module_metamix[[1]]
		ensemble_temp[,"module_alpha"] <- module_alpha[[1]]
		ensemble_temp[,"module_beta"] <- module_beta[[1]]
		ensemble_temp[,"module_relabus"] <- module_relabus[[1]]
		ensemble_temp[,"module_curated1"] <- module_curated1[[1]]
		ensemble_temp[,"module_curated2"] <- module_curated2[[1]]
		ensemble_temp[,"module_trinaryfamily"] <- module_trinaryfamily[[1]]
		ensemble_temp[,"module_species"] <- module_species[[1]]

		print("Training ensemble modules head")
		print(head(ensemble_temp))

		# Part Ib: Test data, individual modules
		output_temp[,"module_baseclin"] <- module_baseclin[[2]]
		output_temp[,"module_agesex"] <- module_agesex[[2]]
		output_temp[,"module_metamix"] <- module_metamix[[2]]
		output_temp[,"module_alpha"] <- module_alpha[[2]]
		output_temp[,"module_beta"] <- module_beta[[2]]
		output_temp[,"module_relabus"] <- module_relabus[[2]]
		output_temp[,"module_curated1"] <- module_curated1[[2]]
		output_temp[,"module_curated2"] <- module_curated2[[2]]
		output_temp[,"module_trinaryfamily"] <- module_trinaryfamily[[2]]
		output_temp[,"module_species"] <- module_species[[2]]
			
		print("Temp test prediction output head")
		print(head(output_temp))


		catsystime("-- Final modules to include, nested regularization --")
		w1 <- grep("module", colnames(ensemble_temp), value = TRUE)
		w2 <- grep("module", colnames(output_temp), value = TRUE)
		# Return scores from this particular seed; standardize output via z-scores
		tmp <- module_glmnet(trainx = ensemble_temp[,w1], trainy = train_y, test = output_temp[,w2])

		# Store identified coefficients per module, and final element as the identified modules
		# Post-challenge phase; store coefficients inside each module
		coefs[[length(coefs)+1]] <<- list(
			coefs_module_baseclin = module_baseclin[[3]],
			coefs_module_agesex = module_agesex[[3]],
			coefs_module_metamix = module_metamix[[3]],
			coefs_module_alpha = module_alpha[[3]],
			coefs_module_beta = module_beta[[3]],
			coefs_module_relabus = module_relabus[[3]],
			coefs_module_curated1 = module_curated1[[3]],
			coefs_module_curated2 = module_curated2[[3]],
			coefs_module_trinaryfamily = module_trinaryfamily[[3]],
			coefs_module_species = module_species[[3]],
			coefs_modules = tmp[[3]]
		)

		print("Aggregated nested regularized module prediction for test data prior to z-scaling:")
		print(tmp[[2]])

		zscale(tmp[[2]])
	}))

	# Head of scores over different seeds
	print("Head of scores across random seeds (cols=seeds, rows=patients)")
	print(head(scores))		
	output_final[,"Score"] <- apply(scores, MARGIN=1, FUN=mean) # Take mean across predicted scores across random seeds		

	print("Head of the mean of scores across seeds")
	print(head(output_final))

	# Scale within [0,1] as instructed
	# From instructions: "The predictions have to be between 0 and 1, with larger numbers 
	# being associated with higher likelihood of having HF"
	output_final[,"Score"] <- outscale(output_final[,"Score"])

	print("Score head after scaling [0,1]")
	print(head(output_final))

	# Runtime sanity checking
	diff_time <- Sys.time() - start_time
	cat("Total runtime taken:\n")
	cat(paste(as.character(diff_time), attr(diff_time, "units")))
	cat("\n\n")
	# Return scores df
	output_final		
}

###
#
# Helper functions
#
###

# Variable name sanitization; cap character length and remove special problematic symbols
sanitize <- \(x, charlimit = 50){
	# Crop to the end of the character string
	x <- c(unlist(lapply(x, FUN=\(q){
		substr(q, start=max(1, nchar(q)-charlimit+1), stop=nchar(q))
	})))
	# Sanitize special symbols
	gsub(":", "", x)	
}

# Scale risk scores between [0,1]
outscale <- \(x){ (x - min(x, na.rm=TRUE))/(max(x, na.rm=TRUE) - min(x, na.rm=TRUE)) }
# Just left shift variables to start at zero
shift <- \(x){ (x-min(x, na.rm=TRUE)) }
# Scale variables with z transformation
zscale <- \(x){
	if(stats::sd(x, na.rm=TRUE) == 0){
		rep(0, times=length(x)) # Singular vector, returning unit values
	}else{
		(x-mean(x, na.rm=TRUE))/stats::sd(x, na.rm=TRUE) # sd != 0, can apply z-score normalization
	}
}
# Left shifted and squared
sqshift <- \(x){ (x-min(x, na.rm=TRUE))^2 }
# Left shifted and n*log(n+1) multiplied
nlognshift <- \(x){ x <- x-min(x, na.rm=TRUE); x*log(x+1) }

# Modified helper reader function modified from microbiome::read_csv2phyloseq
# package 'microbiome' in bioconductor, also using 'phyloseq'
#
# Taxonomy table provided in 'taxtable.csv' missing row names, so needed to hack the function a bit
#
csv2phylo <- function (otu.file = NULL, taxonomy.file = NULL, metadata.file = NULL, sep = ",") 
{
    s.meta <- read.csv(metadata.file, row.names = 1, check.names = FALSE, sep = sep)
    s.sampledata <- phyloseq::sample_data(s.meta)
    s.otu <- read.csv(otu.file, row.names = 1, check.names = FALSE, sep = sep)
    if (any(rownames(s.otu) %in% rownames(s.meta))) {
        s.otu <- t(s.otu)
    }
    s.otu.table <- phyloseq::otu_table(s.otu, taxa_are_rows = TRUE)
    
    s.tax_table <- phyloseq::tax_table(as.matrix(read.csv(taxonomy.file, sep = sep)))
    rownames(s.tax_table) <- rownames(s.otu)

    if (!all(rownames(s.tax_table) == s.tax_table[, ncol(s.tax_table)])) {
        s.tax_table <- cbind(s.tax_table, OTU = rownames(s.tax_table))
        s.tax_table <- phyloseq::tax_table(s.tax_table)
    }
    pseq <- phyloseq::merge_phyloseq(s.otu.table, s.tax_table, s.sampledata)
    return(pseq)
}

## From: https://gist.github.com/variani/d6a42ac64f05f8ed17e6c9812df5492b
#
#' Inverse normal transformation (INT)
#'
#' https://www.biostars.org/p/80597/
#' See the supplement of Yang et al. Nature 2012. 
#'
#' @example
#' x1 <- 5:1
#' inormal(x1)
#'
#' x2 <- c(NA, 5:1, 5:1, NA)
#' inormal(x2)
inormal <- function(x)
{
  qnorm((rank(x, na.last = "keep") - 0.5) / sum(!is.na(x)))
}

## Take all pairwise interactions and create derived variables, and pairwise interactions from two data matrices
#
# From old code: https://github.com/Syksy/ePCR/blob/master/R/int.R
interact.all <- function(input){
	output <- do.call("cbind", lapply(1:ncol(input), FUN=function(z){ 
		do.call("cbind", lapply(z:ncol(input), FUN=function(x){
			tmp <- data.frame(input[,z] * input[,x])
			colnames(tmp)[1] <- paste(colnames(input)[z], "_X_", colnames(input)[x], sep="") 
			tmp
		}))
	}))
	output
}
interact.part <- function(input, first, second){
	output <- do.call("cbind", lapply(first, FUN=function(z){ 
		do.call("cbind", lapply(second, FUN=function(x){
			tmp <- data.frame(input[,z] * input[,x])
			colnames(tmp)[1] <- paste(z, "_X_", x, sep="") 
			tmp
		}))
	}))
	output
}


## Collapse microbiome data to a desired level
#'
#' Take for example: "k__Bacteria;p__Firmicutes;c__Bacilli;o__Lactobacillales;f__Streptococcaceae;g__Streptococcus;s__Streptococcus_gallolyticus"
#' ... could be collapsed to levels k, p, c, o, f, g binding all rows underneath together
#' Presumably by summing for aligned read counts.
#'
collapseMicrobiome <- function(
	x, # Input matrix
	stratFUN = colnames, # The extractor function for strata, by default column names
	sep = ";", # Separator for the various strata
	strata = c("k", "p", "c", "o", "f", "g"), # Chosen level to collapse to; 's' == 'species' is the lowest level identified, thus omitted. Can also be integer, with k==1, p==2, c==3, ...
	collapseFUN = sum, # Function to apply over rows to collapse, applied across columns within chosen strata
	...
){
	stratas <- c("k", "p", "c", "o", "f", "g")
	if(is.integer(strata) | is.numeric(strata)){
		strata <- stratas[strata]
	}else if(!any(strata %in% stratas)){
		stop("Invalid strata; please indicate letter 'k', 'p', 'c', 'o', 'f', or 'g' or an integer corresponding to these in increasing order")
	}
	
	# Collapse matrix
	mat <- do.call("rbind", lapply(stratFUN(x), FUN=\(q){ strsplit(q, ";")[[1]] }))
	colnames(mat) <- c("k", "p", "c", "o", "f", "g", "s")
	
	# Collapsed names per columns
	nams <- apply(mat[,1:which(colnames(mat) == strata),drop=FALSE], MARGIN=1, FUN=\(x){ paste0(x, collapse=";") })
	# Collapse matrix based on these indices
	cmat <- do.call("cbind", by(t(x), INDICES=nams, FUN=\(q){ apply(q, MARGIN=2, FUN=collapseFUN) }))

	# Return collapsed matrix
	cmat
}

#' Relatively straight-forward filter based on threshold values and proportion of samples that need to qualify
filterMicrobiome <- function(
	x, # Data matrix x
	MARGIN = 2, # MARGIN to filter on (passed to apply etc; 1 = rows, 2 = columns)
	prop = 0.10, # Proportion of samples required to have at least <threshold> value for variable to be included
	threshold = 5, # Threshold value that samples ought to >= above
	...
){
	# 'which' rows/cols to keep
	w <- which(apply(x, MARGIN=MARGIN, FUN=\(q){
		(sum(q>=threshold) / length(q)) >= prop
	}))
	# Return filtered rows or columns
	switch(MARGIN, 
		"1" = x[w,,drop=FALSE],
		"2" = x[,w,drop=FALSE]
	)
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

###
#
# Run script inside the environment
#
###

# Submission number
subv <- 8
subname <- "Post2"

args=(commandArgs(TRUE))
PARAM <- list()
# Path of input folder
PARAM$folder.R <- paste0(args[1]) 

# Just in case for debugging
if(length(args) < 2){
	stop("Two arguments should be passed to the R script; first the output folder location, second either 'test' or 'scoring'")
}else{
	print("args provided as seen by the R script:")
	print(args)
}

# Once model actually runs
# Create team folder
dir.create(file.path(PARAM$folder.R, paste0("DenverFINRISKHacky_", subname)))
# Create team output folder
dir.create(file.path(PARAM$folder.R, paste0("DenverFINRISKHacky_", subname),"output"))
PARAM$folder.data <- paste0(PARAM$folder.R, "/")
PARAM$folder.result <- paste0(PARAM$folder.data, paste0("DenverFINRISKHacky_", subname), "/output/")

## Final submission structure;
# Training data is the old 'training'
# New 'test' data is the 'scoring' data

# Pheno data (both meta as well as response)
# 2nd par determines which test set to use ('test' or 'scoring'):
if(args[2] == "scoring"){
	test_p <- read.csv(file = paste0(PARAM$folder.data, "scoring/pheno_scoring.csv"), row.names=1)
}else if(args[2] == "test"){
	test_p <- read.csv(file = paste0(PARAM$folder.data, "test/pheno_test.csv"), row.names=1)
}else{
        stop("Second input argument to R script should be either 'scoring' or 'test'")
}
train_p <- read.csv(file = paste0(PARAM$folder.data, "train/pheno_training.csv"), row.names=1)

# Read count raw data
# 2nd par determines which test set to use ('test' or 'scoring'):
if(args[2] == "scoring"){
	test_r <- read.csv(file = paste0(PARAM$folder.data, "scoring/readcounts_scoring.csv"), row.names=1)
}else if(args[2] == "test"){
	test_r <- read.csv(file = paste0(PARAM$folder.data, "test/readcounts_test.csv"), row.names=1)
}else{
	stop("Second input argument to R script should be either 'scoring' or 'test'")
}
train_r <- read.csv(file = paste0(PARAM$folder.data, "train/readcounts_training.csv"), row.names=1)

# Read taxa along with raw reads and metadata into a phyloseq object
if(args[2] == "test"){
        test_phylo <- csv2phylo(
                otu.file=paste0(PARAM$folder.data, "test/readcounts_test.csv"), # -> to test
                taxonomy.file=paste0(PARAM$folder.data, "test/taxtable.csv"), # -> to test
                metadata.file=paste0(PARAM$folder.data, "test/pheno_test.csv") # -> to test
        )
}else if(args[2] == "scoring"){
	test_phylo <- csv2phylo(
		otu.file=paste0(PARAM$folder.data, "scoring/readcounts_scoring.csv"), # -> to scoring
		taxonomy.file=paste0(PARAM$folder.data, "scoring/taxtable.csv"), # -> to scoring
		metadata.file=paste0(PARAM$folder.data, "scoring/pheno_scoring.csv") # -> to scoring
	)
}else{
	stop("Second input argument to R script should be either 'scoring' or 'test'")
}
train_phylo <- csv2phylo(
	otu.file=paste0(PARAM$folder.data, "train/readcounts_training.csv"), 
	taxonomy.file=paste0(PARAM$folder.data, "train/taxtable.csv"),
	metadata.file=paste0(PARAM$folder.data, "train/pheno_training.csv")
)

# Run model and obtain scores result
res <- model(
	train_clin = train_p, # Training clinical metadata
	test_clin = test_p, # Test clinical metadata
	train_micr = train_r, # Training microbiome raw read counts
	test_micr = test_r, # Testing microbiome raw read counts
	train_phyloseq = train_phylo, # Training phyloseq object
	test_phyloseq = test_phylo, # Train phyloseq object
	v = subv, # Submission version
	a = 1, # Changing to LASSO
	seeds = 1:3 # Vector of random seeds to alleviate random binning effects, used for multiple runs
)

# Provide debug output
print("Final head of results:")
print(head(res))

# Write the resulting scores.csv
write.csv(res, file=paste0(PARAM$folder.result, "scores.csv"), quote=FALSE, row.names=FALSE)

# Write the list of lists with penalized linear Cox model coefficients as an RData object
save(coefs, file=paste0(PARAM$folder.result, "coefs.RData"))
