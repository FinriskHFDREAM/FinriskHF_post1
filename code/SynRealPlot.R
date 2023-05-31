# load data
library(readr)
library(reshape2)
library(ggplot2)
library(ggpubr)
library(ggupset)
library(tidyr)
library(patchwork)
library(dplyr)
# set parameters
args=(commandArgs(TRUE))
mainDir <- paste0(args[1])
PARAM <- list()
PARAM$folder.R <- paste0(getwd(), "/")
PARAM$folder.result <- paste0(PARAM$folder.R, "results/output")
PARAM$folder.data <- paste0(PARAM$folder.R, "final_results")
PARAM$folder.synthetic <- paste0(PARAM$folder.R,"input_syn/input_challenge/") #change here accordingly
PARAM$folder.finrisk <- paste0(PARAM$folder.R,"input_real/") #change here accordingly
###############################################################################
# comparing synthetic vs real in different groups
###############################################################################
# finrisk
S.test.F <- read.csv(file = paste0(PARAM$folder.finrisk, 
                                   "test/pheno_test.csv"), 
                     row.names=1,)
S.train.F <- read.csv(file = paste0(PARAM$folder.finrisk, 
                                    "train/pheno_training.csv"),
                      row.names=1,)
S.score.F <- read.csv(file = paste0(PARAM$folder.finrisk, 
                                    "scoring_nohide/pheno_scoring_nohide.csv"), 
                      row.names=1,)
S.test.F$data <- "FINRISK"
S.train.F$data <- "FINRISK"
S.score.F$data <- "FINRISK"

S.test.F$Group <- "Test"
S.train.F$Group <- "Train"
S.score.F$Group <- "Score"
# synthetic
S.test.S <- read.csv(file = paste0(PARAM$folder.synthetic, 
                                   "test/pheno_test.csv"), 
                     row.names=1,)
S.train.S <- read.csv(file = paste0(PARAM$folder.synthetic, 
                                    "train/pheno_training.csv"),
                      row.names=1,)
S.score.S <- read.csv(file = paste0(PARAM$folder.synthetic, 
                                    "scoring/pheno_scoring.csv"), 
                      row.names=1,)

S.test.S$data <- "Synthetic"
S.train.S$data <- "Synthetic"
S.score.S$data <- "Synthetic"

S.test.S$Group <- "Test"
S.train.S$Group <- "Train"
S.score.S$Group <- "Score"

# bind all tables
phe.table <-rbind(S.test.S, S.train.S, S.score.S, S.test.F, S.train.F, S.score.F)
df=melt(phe.table)
###############################################################################
# plot categorical variables
library("viridis")
p.bar <- phe.table %>%
  select(Group, data, Sex, Smoking, BPTreatment, PrevalentCHD,
         PrevalentHFAIL, PrevalentDiabetes, Event) %>%
  melt() %>%
  count(Group, data, variable, value) %>%
  mutate(Groupset = paste0(Group, " ", data)) %>%
  ggplot(aes(x = Group, 
             y = n, 
             fill = Groupset, 
             color=as.factor(value))) +
  geom_bar(stat="identity", 
           position=position_dodge()) +
  facet_wrap(variable~., scales = "free")+
  scale_color_viridis(discrete = TRUE, option = "D")+
  scale_fill_viridis(discrete = TRUE) +
  labs(y = "", 
       x="", 
       color = "Group",
       fill = "Dataset",
       title = "Comparison of FINRISK and Synthetic Daatset in subgroupos")+
  theme_classic() 


p.bar
ggsave(filename=paste0(PARAM$folder.result, "barplot_categorical.pdf"),
       plot=p.bar) 

###############################################################################
# same but density plot
p.density <- phe.table %>%
  select(Group, data, Age, Event_time, BodyMassIndex, SystolicBP, NonHDLcholesterol) %>%
  melt() %>%  
  mutate(Group = paste0(data, " ", Group)) %>%
  ggplot() + 
  #geom_histogram(aes(y=..density..), colour="black", fill="white") +
  geom_density(mapping = aes(x = value, 
                             y = ..count.., 
                             color = as.character(Group)),
               alpha = 0.8,
               position = "identity",
               size=1) + 
  facet_wrap(.~variable, scales = "free", nrow = 4) +
  labs(y = "Total # of Individuals", 
       color = "Dataset", 
       fill = "Dataset",
       title = "Comparing two datasets:0 for FINRISK and 1 for Synthetic") +
  scale_color_viridis(discrete = TRUE, option = "D")+
  scale_fill_viridis(discrete = TRUE) +
  theme_classic() +  theme(legend.position="bottom") 

p.density
ggsave(filename=paste0(PARAM$folder.result, "density_continous.pdf"),
       plot=p.density) 
