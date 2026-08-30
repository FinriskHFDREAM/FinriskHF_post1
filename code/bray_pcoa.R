##########################################################################
# PCoA / Bray-Curtis comparison of Train, Test and Score microbiome datsets
#
Purpose:
  #   Asses whether microbiome community composition is comparable between
  #  training, testing and scoring datsets.
  #
  # Works for both REAL and SYNTHETIC datasets provided tha theinput files
  # have thsame structure and variable names.
  ####################################################################

rm(list = ls())

#####################################################################
# 1. Packages
##################################################################

library(tidyr)
library(ggrepel)
library(pROC)
library(vegan)
library(reshape2)
library(ggplot2)
library(ggpubr)
library(car)
library(dplyr)
library(plyr)
library(glue)
library(tibble)
library(phyloseq)
library(ggExtra)
library(patchwork)

###############################################################
# 2. Input arguments
#############################################################

args <- commandArgs(TRUE)
mainDir <- paste0(args[1])

PARAM <- list()

PARAM$folder.R <- paste0(mainDir, "/")

PARAM$folder.results <- paste0(
  PARAM$folder.R,
  "results/"
)
# Create results directory if itdoes not exist
if (!dir.exists(PARAM$folder.results)) {
  dir.create(
    PARAM$folder.result,
    recursive = TRUE
  )
}

#############################################################
# 3. Load phenotype data
###################################################################

print("Loading phenotype data...")
S.test <- read.csv(
  file = paste0(
    PARAM$folder.R,
    "input_real/test/pheno_test.csv"
  ),
  row.names = 1,
  check.names = FALSE
)

S.train <- read.csv(
  file = paste0(
    PARAM$folder.R,
    "input_real/train/pheno_training.csv"
  ),
  row.names = 1,
  check.names = FALSE
)

S.score <- read.csv(
  file = paste0(
    PARAM$folder.R,
    "input_real/scoring_nohide/pheno_scoring_nohide.csv"
  ),
  row.names = 1,
  check.names = FALSE
)
################################################################
# 4. Add datset labels
#######################################################################

S.test$Group <- "Test"
S.train$Group <- "Train"
S.score$Group <- "Score"

########################################################################
# 6. Combine phenotype data
###################################################################
phe.table <-rbind(
  S.train,
  S.test,
  S.score
)
# Make sure row names are preserved as ample IDs
phe.table$ID <- rownames(phe.table)

###################################################################
# 7.Load microbiome read-count data
#####################################################################
print("Loading microbiome data...")

O.train <- read.csv(
  file = paste0(
    PARAM$folder.R,
    "input_real/train/readcounts_training.csv"
  ),
  row.names = 1,
  check.names = FALSE
)

O.test <- read.csv(
  file = paste0(
    PARAM$folder.R,
    "input_real/test/readcounts_test.csv"
  ),
  row.names = 1,
  check.names = FALSE
)

O.score <-read.csv(
  file = paste0(
    PARAM$folder.R,
    "input_real/scoring_nohide/readcounts_scoring.csv"
  ),
  row.names =1,
  check.names = FALSE
)

############################################################
# 8. Combine microbiome tables
#
# Columns = samples
# Rows    = taxa
##############################################################

otu.total <- cbind(
  O.train,
  O.test,
  O.score
)
###############################################################
# 9. Check that phenotype and microbiome sample IDs match
###################################################################

common_samples <- intersect(
  colnames(otu.total),
  rownames(phe.table)
)

print(
  paste(
    "Number of samples with both phenotype and microbiome data:",
    length(common_samples)
  )
)
# Keep only common samples and putthem in the same orde
otu.total <- otu.total[, common_samples, drop = FALSE]

phe.table <- phe.table[
  common_samples,
  ,
  drop = FALSE
]

###################################################################
# 10.Create taxonomy table
#
# Taxonomic names are assumed to be separated by ";"
####################################################################

print("Creating axonmy table...")

taxa_names_original <- rownames(otu.total)

taxtable <- strsplit(
  taxa_names_original,
  ";"
)
taxtable <- matrix(
  unlist(taxtable),
  nrow = length(taxtable),
  byrow = TRUE
)
# Make sure the taxonomy table has sevn columns
expected_taxa_levels <- c(
  "Domain",
  "Phylum",
  "Class",
  "Order",
  "Family",
  "Genus",
  "Species"
)

if (ncol(taxtable) < length(expected_taxa_levels)) {
  
  taxtable <- cbind(
    taxtable,
    matrix(
      NA,
      nrow = nrow(taxtable),
      ncol = length(expected_taxa_levels) - ncol(taxable)
    )
  )
}

if (ncol(taxtable) > length(expected_taxa_levels)) {
  
  taxtable <- taxable[
    ,
    seq_along(expected_taxa_leves),
    drop = FALSE
  ]
}
colnames(taxtable) <- expected_taxa_levels
rownames(taxtable) <- taxa_names_original

################################################################
# 11. Create phyloseq object
###############################################################

phyOb <- phyloseq(
  otu_table(
    as.matrix(otu.total),
    taxa_are_rows = TRUE
  ),
  tax_table(
    taxtable
  ),
  sample_data(
    phe.table
  ))


##################################################################
# 12. Aggregate to species level
################################################################

print("Aggregating microbiome data to species lvel...")

physec <- tax_glom(
  phyOb,
  taxrank = "Species",
  NArm = FALSE
)

######################################################################
# 13. Remove taxa without species-level assignment
##################################################################

tax_species <- as.data.frame(
  tax_table(physec)
)

keep_species <- !is.na(tax_species$Species) &
  tax_species$Species != "" &
  tax_species$Species != "s__"

physec <- prune_taxa(
  keep_secies,
  physec
)
print(
  paste(
    "Number of species retained:",
    ntaxa(physec)
  )
)
###############################################################
# 14. Extract species abundance matrix
#
# phyloseq convention:
#   rows = taxa
#   columns = samples
####################################################################

otu_spec <- as(
  otu_table(physec),
  "matrix"
)

if (!taxa_are_rows(physec)) {
  otu_spec <- t(otu_spec)
}

################################################################
# 15. Convert counts to relative abundace
#
# Bray-Curtis can be calculated on relative abundances.
#####################################################################

sample_totals <- colSums(
  otu_spec,
  na.rm = TRUE
)

# Remove samples with zero toal abundance
keep_samples <-sample_totals > 0

otu_spec <- otu_spec[
  ,
  keep_samples,
  drop = FALSE
]
sample_totals <-sample_totals[
  keep_samples
]

otu_spec_rel <- sweep(
  otu_spec,
  2,
  sample_totals,
  FUN = "/"
)

##################################################################
# 16. Make sure phenotype data are aligned
####################################################################

comon_samples <- intersect(
  colnames(otu_spec_rel),
  rownames(phe.table)
)
otu_spec_rel <- otu_spec_rel[
  , common_samples,
  drop = FALSE
]

phe.table <- phe.table[
  common_samples,
  ,
  drop = FALSE
]
# Check alignment
stopifnot(
  identical(
    colnames(otu_spec_rel),
    rownames(phe.table)
  )
)
#############################################################
# 17. Bray-Curtis dissmilarity
###############################################################

print("Calculating Bray-Curtis dissimilarity...")

# vegdist expects:
#   rows = samples
#   columns = species

otu_samples <- t(
  otu_spec_rel
)
beta_dist <- vegan::vegdist(
  otu_samples,
  method = "bray"
)
saveRDS(
  beta_dist,
  file = paste0(
    PARAM$folder.results,
    "bray_curtis_species.rds"
  )
)

##############################################################
# 18. PERMANOVA
#
# Tests wheter community composition differs betwen Train/Test/Score.
####################################################################
print("Running PERMANOVA...")

adonis_result <- vegan::adonis2(
  beta_dist ~ Group,
  data = phe.table,
  permutations = 999
)

print(adonis_result)

write.csv(
  as.data.frame(adonis_result),
  file = paste0(
    PARAM$folder.results,
    "permanova_group.csv"
  ),
  quote = FALSE
)

################################################################
# 19. PCoA
#################################################################

print("Running PCoA...")

pcoa <- cmdscale(
  beta_dist,
  eig = TRUE,
  add = TRUE
)

positions <- pcoa$points

colnames(positions) <- c(
  "PCoA1",
  "PCoA2"
)

###################################################################
# 20. Percentage variance explained
##################################################################

eig <- pcoa$eig

# Only positve eigenvalues contribue to explained variance
positive_eig <- eig[eig > 0]

percent_explained <-(
  100 * eig / sum(positive_eig)
)
pretty_pe <- format(
  round(
    percent_explained,
    digits = 1
  ),
  nsmall = 1,
  trim = TRUE
)

x_label <- paste0(
  "PCoA 1 (",
  pretty_pe[1],
  "%)"
)

y_label <- paste0(
  "PCoA 2 (",
  pretty_pe[2],
  "%)"
)

################################################################
# 21. Combine PCoA coordinates with metadat
####################################################################

positions <- as.data.frame(
  positions
)

positions$ID <- rownames(
  positions
)

pcoa_data <- dplyr::inner_join(
  phe.table,
  positions,
  by= "ID"
)

saveRDS(
  pcoa_data,
  file = paste0(
    PARAM$folder.results,
    "pcoa_species.rds"
  )
)

###############################################################
# 22. Set Group as factor
###################################################################

pcoa_data$Group <- factor(
  pcoa_data$Group,
  levels = c(
    "Train",
    "Test",
    "Score"
  )
)

###############################################################
# 23. Density-contour PCoA plot
#
# Points show indivdual samples.
# Filled density conours show where samples are conentrated.
##################################################################

print("Creating density contur PCoA plot..")

#PLOT Density outside    
cb_palette <-c(
  "Train" = "#1b9177",
  "Test"  = "#d95f02",
  "Score" = "#7570b3"
)

p.pcoa <- ggplot(
  pcoa_data,
  aes(
    x = PCoA1,
    y = PCoA2,
    color = Group
  )
) +
  
  geom_point(
    alpha = 0.4,
    size = 1
  ) +
  
  stat_density_2d(
    aes(color = Group),
    linewidth = 0.8
  )+
  
  scale_color_manual(values = cb_palette) +
  
  labs(
    x = x_label,
    y = y_label,
    color = "Dataset"
  ) +
  
  theme_classic() +
  
  theme(
    # Axis labels
    axis.title.x = element_text(size = 20),
    axis.title.y = element_text(size = 20),
    
    # Axis tick labels
    axis.text.x = element_text(size = 16),
    axis.text.y = element_text(size = 16),
    # Legend
    legend.title = element_text(size = 18),
    legend.text = element_text(size = 16),
    
    # Plot title, if you add one
    plot.title = element_text(
      size = 16,
      hjust = 0.5
    )
  )

p.x <- ggplot(
  pcoa_data,
  aes(x = PCoA1, fill = Group, color = Group)
) +
  geom_density(alpha = 0.25) +
  scale_color_manual(values = cb_palette) +
  scale_fill_manual(values = cb_palette) +
  theme_void() +
  theme(
    legend.position = "none"
  )

p.y <- ggplot(
  pcoa_data,
  aes(x = PCoA2, fill = Group, color = Group)
) +
  geom_density(alpha = 0.25) +
  scale_color_manual(values = cb_palette) +
  scale_fill_manual(values = cb_palette) +
  theme_void() +
  theme(
    legend.position = "none"
  )  + coord_flip()
p.final <- p.x + plot_spacer() +
  p.pcoa + p.y +
  plot_layout(
    widths = c(4, 1),
    heights =c(1, 4)
  )

ggsave(
  filename = paste0(
    PARAM$folder.result,
    "pcoa_species_density_real.pdf"
  ),
  plot = p.final,
  width = 7,
  height =6
)
