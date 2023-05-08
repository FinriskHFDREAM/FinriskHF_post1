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
sb2 <- read.csv(sb2_path, header = TRUE)

# DENVER scores
denver <- read.csv(denver_path, header = TRUE)

# Join datasets and stratify by clinical covariates
data <- 
  test %>% 
  inner_join(., sb2, by = "SampleID") %>% 
  inner_join(., denver, by = "SampleID", suffix = c("_sb2", "_denver")) %>% 
  mutate(
    surv = Surv(time = Event_time, event = Event),
    Age_str = cut2(Age, m = 300),
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
    holmes_denver = teststats(surv, Score_denver)$hoslem_p$pval,
    holmes_sb2 = teststats(surv, Score_sb2)$hoslem_p$pval,
    c_denver = teststats(surv, Score_denver)$C,
    c_sb2 = teststats(surv, Score_sb2)$C
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
    holmes_denver = teststats(surv, Score_denver)$hoslem_p$pval,
    holmes_sb2 = teststats(surv, Score_sb2)$hoslem_p$pval,
    c_denver = teststats(surv, Score_denver)$C,
    c_sb2 = teststats(surv, Score_sb2)$C
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
    holmes_denver = teststats(surv, Score_denver)$hoslem_p$pval,
    holmes_sb2 = teststats(surv, Score_sb2)$hoslem_p$pval,
    c_denver = teststats(surv, Score_denver)$C,
    c_sb2 = teststats(surv, Score_sb2)$C
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
    holmes_denver = teststats(surv, Score_denver)$hoslem_p$pval,
    holmes_sb2 = teststats(surv, Score_sb2)$hoslem_p$pval,
    c_denver = teststats(surv, Score_denver)$C,
    c_sb2 = teststats(surv, Score_sb2)$C
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
  ylab(label = "C-Index")+theme(axis.text.x = element_text(size = 11))  
ggsave( file=(paste0(PARAM$folder.data,"results/prelim_test/Harrel_C_compare.pdf")),width = 14,height = 7, device="pdf")

toPlot_h <- toPlot[which(toPlot$metric == "holmes"), ]
ggscatter(
  toPlot_h,
  x = "groups",
  y = "value",
  color = "team"
) +
  facet_wrap(~variable , scales = "free") +
  ylab(label = "Hoslem test")+theme(axis.text.x = element_text(size = 11))   
ggsave(file=(paste0(PARAM$folder.data,"results/prelim_test/Hoslem_compare.pdf")),width = 14,height = 7, device="pdf")