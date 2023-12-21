#load library
library(knitr)
library(tidyr)
library(ggplot2)
library(ggpubr)
library(car)
library(vegan)
library(glue)
require(Hmisc)
library(dplyr)
library(plyr)
library("tibble")
library(microbiome)
library(phyloseq)

args=(commandArgs(TRUE))
mainDir <- paste0(args[1])
PARAM <- list()
PARAM$folder.R <- paste0(mainDir, "/")
PARAM$folder.results <- paste0(PARAM$folder.R, "results/pcoa_spec_prune2/")

#open the PCoA calc

#phfs.rel.ord.pcoa <- readRDS(paste0(PARAM$folder.results,
#                   "betadis.rds"))

values=c("#1976D2", "#009688", "#FFDB6D", "#C4961A", "#CC79A7", 
         "#D16103", "#C3D7A4", "#52854C", "#4E84C4", "#293352")

#values =c("#e41a1c","#377eb8","#4daf4a", "#984ea3", "#ff7f00", "#a65628"))
# load data
S.test <- read.csv(file = paste0(PARAM$folder.R, 
                                   "input_real/test/pheno_test.csv"), 
                     row.names=1,)
S.train <- read.csv(file = paste0(PARAM$folder.R, 
                                    "input_real/train/pheno_training.csv"),
                      row.names=1,)
S.score <- read.csv(file = paste0(PARAM$folder.R, 
                                    "input_real/scoring_nohide/pheno_scoring_nohide.csv"), 
                      row.names=1,)
S.test$Group <- "Test"
S.train$Group <- "Train"
S.score$Group <- "Score"
#DEFINE AGE_GROUPS
intervals <- c(24, 34, 44, 54, 64, 75)
S.train$age_group = cut(S.train$Age, breaks = intervals)
S.test$age_group = cut(S.test$Age, breaks = intervals)
S.score$age_group = cut(S.score$Age, breaks = intervals)
#S.score$age_group=xxx

O.train <- read.csv(file = paste0(PARAM$folder.R , "input_real/train/readcounts_training.csv"),
               row.names=1,
               check.names=FALSE)

O.test <- read.csv(file = paste0(PARAM$folder.R , "input_real/test/readcounts_test.csv"),
                    row.names=1,
                    check.names=FALSE)

O.score <- read.csv(file = paste0(PARAM$folder.R , "input_real/scoring_nohide/readcounts_scoring.csv"),
                    row.names=1,
                    check.names=FALSE)

# bind all tables
phe.table <- rbind(S.train, S.test, S.score)
otu.total  <- cbind(O.train, O.test, O.score)

#taxtable
taxtable <- strsplit(row.names(as.matrix(otu.total)), ";")
taxtable <- matrix(unlist(taxtable), nrow=length(taxtable), byrow=T)
row.names(otu.total) <- row.names(taxtable)
phyOb <- phyloseq(otu_table(as.matrix(otu.total), taxa_are_rows=T), tax_table(taxtable))
to_be_pruned <- sample_sums(phyOb) > 50000
phyOb <- prune_samples(to_be_pruned, phyOb)
colnames(tax_table(phyOb)) <- c("Domain", "Phylum", "Class", "Order", "Family",  "Genus", "Species")

#Combine the phenotype data with the taxa data
phyOb <- phyloseq(otu_table(phyOb), tax_table(phyOb), sample_data(phe.table))
phyOb <- subset_samples(phyOb, Event %in% c(0, 1))
physpec <- tax_glom(phyOb, taxrank =  "Species")
physpec <- subset_taxa(physpec, Species != "s__")

phf.rel <- microbiome::transform(physpec, 'compositional')

#calculate pcoa
if (file.exists(paste0(PARAM$folder.results,
                  "pcoa_betadis.rds"))) {
  print("file exists")
  phfs.pseq.rel.pcoa <- readRDS(paste0(PARAM$folder.results,
                   "pcoa_betadis.rds"))
} else {
phfs.pseq.rel.pcoa <- ordinate(phf.rel, method="PCoA", distance="bray")
saveRDS(phfs.pseq.rel.pcoa,paste0(PARAM$folder.results,
                   "pcoa_betadis.rds"))
}


meta=phe.table#[1:100,]
meta$ID<-rownames(meta)
otu<- abundances(phf.rel)


#calculate bray
if (file.exists(paste0(PARAM$folder.results,
                   "betadis.rds"))) {
  print("file exists")
  bray.dist.m <- readRDS(paste0(PARAM$folder.results,
                   "betadis.rds"))
} else {
bray.dist.m <- vegdist(t(otu), method="bray")
saveRDS(bray.dist.m,paste0(PARAM$folder.results,               
  "betadis.rds"))
}
# how much the abundance of each phyla correlates with the first two PCoA axis?
# also, we skip Plasmid phyla

## ------------------------------------------------------------------------

tt.m <- as(tax_table(phf.rel), "matrix")
phyla.names <- na.omit(unique(tt.m[, "Phylum"]))

phyla.prop.abd <- lapply(phyla.names,
                         function (x) {
           phyl.taxa <- which(tt.m[, "Phylum"] == x & !grepl("Plasmid", tt.m[, "Domain"])) %>%
       as.integer
           # abundances of taxa
           phyl.abud <- apply(otu_table(phf.rel)[phyl.taxa, ], 2, function (x) sum(x))
                         })
     
phyla.df <- do.call(cbind.data.frame, phyla.prop.abd)
ab.names <- lapply(phyla.names, function (x) paste0(x,'.Abundance') )
colnames(phyla.df) <- ab.names
phyla.df[, 1:length(phyla.names)] <- sapply(phyla.df[, 1:length(phyla.names)], as.numeric)
phyla.df$pcoa1 <- phfs.pseq.rel.pcoa$vectors[,1]
phyla.df$pcoa2 <- phfs.pseq.rel.pcoa$vectors[,2]

phyla.df2 <- merge(meta, phyla.df, by = 0)   
df <- phyla.df2
ord.pcoa <- phfs.pseq.rel.pcoa
# Flip orientation for good layout
head(df)
df$pcoa1 = -df$pcoa1
#df$pcoa2 = -df$pcoa2
#viz
#plot pcoa
p.pcoa <- df %>%
ggplot(aes(x=pcoa1, y=pcoa2, 
           color= as.factor(Group)))+ #replace with Event
  geom_point(alpha=0.8, size=0.5) + 
  labs(x = "PCoA 1", y = "PCoA 2") +
  theme_classic() + 
  xlim(-0.50, 0.50) +
  ylim(-0.50, 0.50) +
  scale_color_manual(labels = c("Without incident HF", "With incident HF", "Score set", "Test set", "Train set"),values = values)+
  scale_fill_manual(values = values) +
  stat_ellipse(inherit.aes = F, #in your example should be able to just use stat_ellipse()
               aes(color=as.factor(Event), x=pcoa1, y=pcoa2))+
  theme(legend.title= element_blank())

# extract relative eigenvalues (proportion of variance explained by each axis)
p.pcoa <- p.pcoa + xlab(paste("PCoA 1 (", round(100*ord.pcoa$values$Relative_eig[1],1), "%", ")", sep = ""))
p.pcoa <- p.pcoa + ylab(paste("PCoA 2 (", round(100*ord.pcoa$values$Relative_eig[2],1), "%", ")", sep = ""))


print("pcoa plot1")
ggsave(p.pcoa, 
       file=paste0(PARAM$folder.results,
                   "betadiv_pcoa_group.pdf"), 
       width = 7, height=5)

p.pcoa <- df %>%
ggplot(aes(x=pcoa1, y=pcoa2, 
           color= as.factor(Group),
           #shape = as.factor(Event)
           ))+ #replace with Event
  geom_point(alpha=0.6, size=1) + 
  labs(x = "PCoA 1", y = "PCoA 2") +
  theme_classic() + 
  xlim(-0.50, 0.50) +
  ylim(-0.50, 0.50) +
  scale_color_manual(labels = c("Score set", "Test set", "Train set"),values = values)+
  #scale_shape_discrete(labels=c("Without incident HF", "With incident HF"))+ 
  scale_fill_manual(values = values)+
  theme(legend.title= element_blank())

# extract relative eigenvalues (proportion of variance explained by each axis)
p.pcoa <- p.pcoa + xlab(paste("PCoA 1 (", round(100*ord.pcoa$values$Relative_eig[1],1), "%", ")", sep = ""))
p.pcoa <- p.pcoa + ylab(paste("PCoA 2 (", round(100*ord.pcoa$values$Relative_eig[2],1), "%", ")", sep = ""))
print("pcoa plot2")
ggsave(p.pcoa, 
       file=paste0(PARAM$folder.results,
                   "betadiv_pcoa_group_noeclipseFlip.pdf"), 
       width = 7, height=5)


p.pcoa <- df %>%
ggplot(aes(x=pcoa1, y=pcoa2, 
           color= as.factor(Event)))+ #replace with Event
  geom_point(alpha=0.6, size=0.5) + 
  labs(x = "PCoA 1", y = "PCoA 2") +
  theme_classic() + 
  xlim(-0.50, 0.50) +
  ylim(-0.50, 0.50) +
  scale_color_manual(labels = c("Without incident HF", "With incident HF"),values = values)+
  scale_fill_manual(values = values) +
  theme(legend.title= element_blank())

# extract relative eigenvalues (proportion of variance explained by each axis)
p.pcoa <- p.pcoa + xlab(paste("PCoA 1 (", round(100*ord.pcoa$values$Relative_eig[1],1), "%", ")", sep = ""))
p.pcoa <- p.pcoa + ylab(paste("PCoA 2 (", round(100*ord.pcoa$values$Relative_eig[2],1), "%", ")", sep = ""))
print("pcoa plot3")
ggsave(p.pcoa, 
       file=paste0(PARAM$folder.results,
                   "betadiv_pcoa_group_event.pdf"), 
       width = 7, height=5)

