# Testing the model
# ============


require(dplyr)
require(survival)
require(ggplot2)
require(ggpubr)
require(Hmisc)

args=(commandArgs(TRUE))
PARAM <- list()
PARAM$folder.R <- paste0(args[1]) 
PARAM$folder.data <- paste0(PARAM$folder.R, "/")

test_path <- paste0(PARAM$folder.data, "input_real/train/pheno_training.csv")        # add path to test data
sb2_path <- paste0(PARAM$folder.data, "TEAMS/SB2/output_ori/scores_train.csv")        # add path to sb2 scores
denver_path <- paste0(PARAM$folder.data, "TEAMS/DenverFINRISKHacky/output_ori/scores_train.csv")      # add path to denver scores

source(paste0(PARAM$folder.data, "code/teststats.r"))  # add complete path to this script

test <- read.csv(test_path, header = TRUE, row.names = 1)
test$SampleID <- rownames(test)

# SB2 scores
SB2 <- read.csv(sb2_path, header = TRUE)

# DENVER scores
DFH <- read.csv(denver_path, header = TRUE)
intervals <- c(0, 34, 44, 54, 64, 75)
# Join datasets and stratify by clinical covariates
data <- 
  test %>% 
  inner_join(., SB2, by = "SampleID") %>% 
  inner_join(., DFH, by = "SampleID", suffix = c("_SB2", "_DFH")) %>% 
  mutate(
    surv = Surv(time = Event_time, event = Event),
    Age_str = cut(Age, breaks = intervals),
    BMI_str = cut2(BodyMassIndex, m = 300),
    SystolicBP_str = cut2(SystolicBP, m = 300),
    NonHDLcholesterol_str = ifelse(NonHDLcholesterol > 4, "high", "normal")
  )

# Make stratified predictions
# =====
# By Age
data_age <- 
  data %>% 
  group_by(Age_str) %>% 
  summarise(
    holmes_DFH = teststats(surv, Score_DFH)$hoslem_p$pval,
    holmes_SB2 = teststats(surv, Score_SB2)$hoslem_p$pval,
    c_DFH = teststats(surv, Score_DFH)$C,
    c_SB2 = teststats(surv, Score_SB2)$C
  ) %>% 
  rename(groups = Age_str) %>% 
  mutate(
    variable = "age"
  )

# By BMI
data_bmi <- 
  data %>% 
  group_by(BMI_str) %>% 
  summarise(
    holmes_DFH = teststats(surv, Score_DFH)$hoslem_p$pval,
    holmes_SB2 = teststats(surv, Score_SB2)$hoslem_p$pval,
    c_DFH = teststats(surv, Score_DFH)$C,
    c_SB2 = teststats(surv, Score_SB2)$C
  ) %>% 
  rename(groups = BMI_str) %>% 
  mutate(
    variable = "bmi"
  )

# By Smoking
data_smoking <- 
  data %>% 
  group_by(Smoking) %>% 
  filter(
    !is.na(Smoking)
  ) %>% 
  summarise(
    holmes_DFH = teststats(surv, Score_DFH)$hoslem_p$pval,
    holmes_SB2 = teststats(surv, Score_SB2)$hoslem_p$pval,
    c_DFH = teststats(surv, Score_DFH)$C,
    c_SB2 = teststats(surv, Score_SB2)$C
  ) %>% 
  rename(groups = Smoking) %>% 
  mutate(
    variable = "smoking",
    groups = ifelse(groups == 1, "smoker", "non-smoker")
  )

# By BPTreatment
data_bptreatment <- 
  data %>% 
  group_by(BPTreatment) %>% 
  filter(
    !is.na(BPTreatment)
  ) %>% 
  summarise(
    holmes_DFH = teststats(surv, Score_DFH)$hoslem_p$pval,
    holmes_SB2 = teststats(surv, Score_SB2)$hoslem_p$pval,
    c_DFH = teststats(surv, Score_DFH)$C,
    c_SB2 = teststats(surv, Score_SB2)$C
  ) %>% 
  rename(groups = BPTreatment) %>% 
  mutate(
    variable = "bptreatment",
    groups = ifelse(groups == 1, "treated", "non-treated")
  )
  

# Join data to plot
toPlot <- rbind.data.frame(
  data_age,
  data_bmi,
  data_bptreatment,
  data_smoking
)

print(head(toPlot))

write.csv(toPlot,paste0(PARAM$folder.data,"results/prelim_test/plot_stratify.csv"))

toPlot <- data.table::melt(toPlot)
names(toPlot) <- c("groups", "variable", "score", "value")
toPlot$score <- as.character(toPlot$score)
toPlot$team <- sapply(strsplit(toPlot$score, "_"), "[", 2)
toPlot$metric <- sapply(strsplit(toPlot$score, "_"), "[", 1)

# Plotting
toPlot_c <- toPlot[which(toPlot$metric == "c"), ]
ggscatter(
  toPlot_c,
  x = "groups",
  y = "value",
  color = "team"
) +
  facet_wrap(~variable , scales = "free") +
  ylab(label = "C-Index")+theme(axis.text.x = element_text(size = 8))  
ggsave( file=(paste0(PARAM$folder.data,"results/prelim_test/Harrel_C_compare_stratify.pdf")),width = 14,height = 7, device="pdf")

toPlot_h <- toPlot[which(toPlot$metric == "holmes"), ]
ggscatter(
  toPlot_h,
  x = "groups",
  y = "value",
  color = "team"
) +
  facet_wrap(~variable , scales = "free") +
  ylab(label = "Hoslem test")+theme(axis.text.x = element_text(size = 8))   
ggsave(file=(paste0(PARAM$folder.data,"results/prelim_test/Hoslem_compare_stratify.pdf")),width = 14,height = 7, device="pdf")