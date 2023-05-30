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
intervals <- c(0, 34, 44, 54, 64, 75)
data <- 
  test %>% 
  inner_join(., sb2, by = "SampleID") %>% 
  inner_join(., denver, by = "SampleID", suffix = c("_sb2", "_denver")) %>% 
  mutate(
    surv = Surv(time = Event_time, event = Event),
    Age_str = cut(Age, breaks = intervals),
  ) %>% 
  select(c(surv, Age, Age_str, Score_sb2, Score_denver))

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
  ) %>% 
  rename(
    Age_str = groups
  )

# Join datasets and calculate mean
res <- 
  left_join(data, data_age[, c("Age_str", "holmes_denver", "holmes_sb2", "c_denver", "c_sb2")], by = "Age_str") %>% 
  select(-surv) %>%
  select(-Age)

saveRDS(res, file = paste0(PARAM$folder.data, "results/age_weighed_scores.rds"))
