cran_packages <- c("tidyr", "BiocManager", "knitr", "ggrepel", "pROC", "vegan",
                   "reshape2", "ggplot2", "ggpubr", "car",  "dplyr", "plyr")
for (i in cran_packages) {
  if (!require(i, character.only = TRUE))
    install.packages(i)
}
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
PARAM$folder.results <- paste0(PARAM$folder.R, "results/pcoa_spec/")

# load data
# synthetic
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
otu.total[1:5,1:5]

#taxtable
taxtable <- strsplit(row.names(as.matrix(otu.total)), ";")
taxtable <- matrix(unlist(taxtable), nrow=length(taxtable), byrow=T)
row.names(otu.total) <- row.names(taxtable)
phyOb <- phyloseq(otu_table(as.matrix(otu.total), taxa_are_rows=T), tax_table(taxtable))
  #Format the tax table.
colnames(tax_table(phyOb)) <- c("Domain", "Phylum", "Class", "Order", "Family",  "Genus", "Species")

values=c("#1976D2", "#009688", "#FFDB6D", "#C4961A", "#CC79A7", 
         "#D16103", "#C3D7A4", "#52854C", "#4E84C4", "#293352")

#create phyloseq object --> extract species only
#otu <- otu_table(otu.total, taxa_are_rows = TRUE)
#tax <- tax_table(taxtable)
#sample <- sample_data(phe.table)
#phyOb <- phyloseq(otu, tax, sample)

#Combine the phenotype data with the taxa data
phyOb <- phyloseq(otu_table(phyOb), tax_table(phyOb), sample_data(phe.table))
physpec <- tax_glom(phyOb, taxrank =  "Species")
physpec <- subset_taxa(physpec, Species != "s__")
phfs.rel <- microbiome::transform(physpec, 'compositional')

#value = c("#b35806","##f1a340","#fee0b6","#d8daeb","#998ec3","#542788")
# motu.abs is the absolue counts of species
# head(motus.abs)
# MMPC35551931ST MMPC41376408ST
#'Candidatus Kapabacteria' thiocyanatum [r_10354]              0              0
#Abiotrophia defectiva [r_04788]                               0              0

# head meta
# environment_material status
#MMPC35551931ST	feces [ENVO:00002003]	PC	1	0	
#MMPC41376408ST	feces [ENVO:00002003]	PC	1	0	
#MMPC59659730ST	feces [ENVO:00002003]

# apply filter and check if we loose so much, we can adjust
#motu.abs=otu.total#[1:100,1:1000]
meta=phe.table#[1:100,]
meta$ID<-rownames(meta)


#motu.abs.fil <- motu.abs[rowSums(motu.abs >= 10^-5) >= 2,
#                         colSums(motu.abs) > 1200 ]

#dim(motu.abs)

motu.abs.fil<- abundances(physpec)

dim(motu.abs.fil)
#[1] 5748 1809
#[1] 5303 1807

# calculate bray curtis
motu.t = t(motu.abs.fil)
min_seq = min(rowSums(motu.t))

print("beta_dist")
beta_dist.rar <-t(motu.abs.fil) %>%
  vegan::vegdist(method = "bray" 
                 )
saveRDS(beta_dist.rar,paste0(PARAM$folder.results,
                   "betadis.rds"))
print("nmds_rar")
nmds.rar <- metaMDS(beta_dist.rar) %>% 
  scores(display=c("sites")) %>% 
  as_tibble(rownames="ID")

# combine metadata and betadiv
meta_nmds.rar <- dplyr::inner_join(meta, nmds.rar)
saveRDS(meta_nmds.rar,paste0(PARAM$folder.results,
                   "nmds.rds"))
###############################################################################
print("plotting betadiv")
p.betadiv.ord.rar <- ggplot(meta_nmds.rar, 
                            aes(x = NMDS1, 
                                y = NMDS2, 
                                #group=Event,
                                color = as.factor(Event))) +
                                #color = as.factor(age_group))) +
                                #color = as.factor(Sex))) + # change color here accordingly such as age groups, gender
  geom_point(aes(shape=as.factor(Sex))) +
  stat_ellipse() +
  scale_color_manual(values = values)+
  scale_fill_manual(values = values) +
  theme_classic()

p.betadiv.ord.rar
ggsave(p.betadiv.ord.rar, 
       file=paste0(PARAM$folder.results,
                   "betadiv_nmds_ordination_Event.pdf"), 
       width = 5, height=4)


p.betadiv.ord.rar <- ggplot(meta_nmds.rar, 
                            aes(x = NMDS1, 
                                y = NMDS2, 
                                #group=Event,
                                #color = as.factor(Event))) +
                                color = as.factor(age_group))) +
                                #color = as.factor(Sex))) + # change color here accordingly such as age groups, gender
  geom_point(aes(shape=as.factor(Sex))) +
  stat_ellipse() +
  scale_color_manual(values = values)+
  scale_fill_manual(values = values) +
  theme_classic()

p.betadiv.ord.rar
ggsave(p.betadiv.ord.rar, 
       file=paste0(PARAM$folder.results,
                   "betadiv_nmds_ordination_agegroup.pdf"), 
       width = 5, height=4)


p.betadiv.ord.rar <- ggplot(meta_nmds.rar, 
                            aes(x = NMDS1, 
                                y = NMDS2, 
                                #group=Event,
                                #color = as.factor(Event))) +
                                #color = as.factor(age_group))) +
                                color = as.factor(Sex))) + # change color here accordingly such as age groups, gender
  geom_point(aes(shape=as.factor(Sex))) +
  stat_ellipse() +
  scale_color_manual(values = values)+
  scale_fill_manual(values = values) +
  theme_classic()

p.betadiv.ord.rar
ggsave(p.betadiv.ord.rar, 
       file=paste0(PARAM$folder.results,
                   "betadiv_nmds_ordination_sex.pdf"), 
       width = 5, height=4)

p.betadiv.ord.rar <- ggplot(meta_nmds.rar, 
                            aes(x = NMDS1, 
                                y = NMDS2, 
                                #group=Event,
                                color = as.factor(Group))) +
                                #color = as.factor(Event))) +
                                #color = as.factor(age_group))) +
                                #color = as.factor(Sex))) + # change color here accordingly such as age groups, gender
  geom_point(aes(shape=as.factor(Sex))) +
  stat_ellipse() +
  scale_color_manual(values = values)+
  scale_fill_manual(values = values) +
  theme_classic()

p.betadiv.ord.rar
ggsave(p.betadiv.ord.rar, 
       file=paste0(PARAM$folder.results,
                   "betadiv_nmds_ordination_group.pdf"), 
       width = 5, height=4)

###############################################################################
# prepare pcoa
print("pcoa ")
dist_tbl <- as_tibble(as.matrix(beta_dist.rar), 
                      rownames="samples")

dist_matrix <- dist_tbl %>%
  pivot_longer(cols=-samples, names_to="b", values_to="distances") %>%
  select(samples, b, distances) %>%
  pivot_wider(names_from="b", values_from="distances") %>%
  select(-samples) %>%
  as.dist()

pcoa <- cmdscale(dist_matrix, eig=TRUE, add=TRUE)
positions <- pcoa$points
colnames(positions) <- c("pcoa1", "pcoa2")
percent_explained <- 100 * pcoa$eig / sum(pcoa$eig)

pretty_pe <- format(round(percent_explained, digits =1), 
                    nsmall=1, 
                    trim=TRUE)

labels <- c(glue("PCoA 1 ({pretty_pe[1]}%)"),
            glue("PCoA 2 ({pretty_pe[2]}%)"))

# plot percent explained
p <- tibble(pe = percent_explained,
       axis = 1:length(percent_explained)) %>%
  ggplot(aes(x=axis, y=pe)) +
  geom_line() +
  coord_cartesian(xlim = c(1, 10), ylim=c(0, 10)) +
  scale_x_continuous(breaks=1:10) +
  theme_classic()
ggsave(p,file=paste0(PARAM$folder.results,
                   "percent_explained.pdf"), 
       width = 5, height=5)
# combine metadata and betadiv
positions <- positions %>% 
  as.data.frame() %>% 
  tibble::rownames_to_column( var = "ID")

pcoa.rar <- dplyr::inner_join(meta, positions)
saveRDS(pcoa.rar,paste0(PARAM$folder.results,
                   "pcoa.rar.rds"))
###############################################################################
#please adapt the code to finrisk meta covariates
#plot pcoa
p.pcoa <- pcoa.rar %>%
ggplot(aes(x=pcoa1, y=pcoa2, 
           color= as.factor(Event), #replace with Event
           size= age_group,#replace with Age group
           shape = as.factor(Sex)))+ 
  geom_point(alpha=0.6, size=2) + 
  labs(x=labels[1], y=labels[2]) +
  theme_classic() + 
  scale_shape_manual(values=c(17, 1), name="Sex") + #keep
  scale_size_continuous(name="age_group") + #change to age group
  scale_color_manual(values = values)+
  scale_fill_manual(values = values) +
  ggtitle("Bray Curtis")  +
  stat_ellipse(inherit.aes = F, #in your example should be able to just use stat_ellipse()
               aes(color=as.factor(Event), x=pcoa1, y=pcoa2)) 

p.pcoa
ggsave(p.pcoa, 
       file=paste0(PARAM$folder.results,
                   "betadiv_pcoa.pdf"), 
       width = 5, height=5)


#plot pcoa
p.pcoa <- pcoa.rar %>%
ggplot(aes(x=pcoa1, y=pcoa2, 
           color= as.factor(Group), #replace with Event
           size= age_group,#replace with Age group
           shape = as.factor(Sex)))+ 
  geom_point(alpha=0.6, size=2) + 
  labs(x=labels[1], y=labels[2]) +
  theme_classic() + 
  scale_shape_manual(values=c(17, 1), name="Sex") + #keep
  scale_size_continuous(name="age_group") + #change to age group
  scale_color_manual(values = values)+
  scale_fill_manual(values = values) +
  ggtitle("Bray Curtis")  +
  stat_ellipse(inherit.aes = F, #in your example should be able to just use stat_ellipse()
               aes(color=as.factor(Event), x=pcoa1, y=pcoa2)) 

p.pcoa
ggsave(p.pcoa, 
       file=paste0(PARAM$folder.results,
                   "betadiv_pcoa_group.pdf"), 
       width = 5, height=5)


#plot pcoa
p.pcoa <- pcoa.rar %>%
ggplot(aes(x=pcoa1, y=pcoa2, 
           color= as.factor(age_group), #replace with Event
           size= Group,#replace with Age group
           shape = as.factor(Sex)))+ 
  geom_point(alpha=0.6, size=2) + 
  labs(x=labels[1], y=labels[2]) +
  theme_classic() + 
  scale_shape_manual(values=c(17, 1), name="Sex") + #keep
  scale_size_continuous(name="age_group") + #change to age group
  scale_color_manual(values = values)+
  scale_fill_manual(values = values) +
  ggtitle("Bray Curtis")  +
  stat_ellipse(inherit.aes = F, #in your example should be able to just use stat_ellipse()
               aes(color=as.factor(Event), x=pcoa1, y=pcoa2)) 

p.pcoa
ggsave(p.pcoa, 
       file=paste0(PARAM$folder.results,
                   "betadiv_pcoa_agegroup.pdf"), 
       width = 5, height=5)