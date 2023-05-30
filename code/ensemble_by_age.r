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

# Defining weigths
weigths <- 
  res %>%
  group_by(Age_str) %>% 
  slice(1) %>% 
  select(c(Age_str, holmes_denver, holmes_sb2, c_denver, c_sb2)) %>% 
  mutate(
    # holmes
    holmes_sb2_cat = ifelse(holmes_sb2 > 0.05, "calibrated", "non-calibrated"),
    holmes_denver_cat = ifelse(holmes_denver > 0.05, "calibrated", "non-calibrated"),
    
    # c-index
    better_c = ifelse(c_denver - c_sb2 < 0, "sb2", "denver"),
    
    # weigths
    w = -log10(holmes_sb2) - (-log10(holmes_denver)),
    w_sb2 = ifelse(w < 0, (c_sb2 * 10) + abs(w), c_sb2 * 10),
    w_denver = ifelse(w > 0, (c_denver * 10) + abs(w), c_denver * 10),
  )

# based on c-index
ensemble <- 
  left_join(res, weigths[, c("Age_str", "w_sb2", "w_denver")], by = "Age_str") %>% 
  mutate(
    mean = (Score_sb2 + Score_denver) / 2,
    w_mean = (w_sb2 * Score_sb2 + w_denver * Score_denver) / (w_sb2 + w_denver)
  )

saveRDS(ensemble, file = paste0(PARAM$folder.data, "results/age_weighed_scores.rds"))
