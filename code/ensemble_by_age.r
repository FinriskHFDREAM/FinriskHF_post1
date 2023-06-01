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
    hoslem_denver = teststats(surv, Score_denver)$hoslem_p$pval,
    hoslem_sb2 = teststats(surv, Score_sb2)$hoslem_p$pval,
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
  left_join(data, data_age[, c("Age_str", "hoslem_denver", "hoslem_sb2", "c_denver", "c_sb2")], by = "Age_str") %>% 
  select(-surv) %>%
  select(-Age)

# Defining weigths
weigths <- 
  res %>%
  group_by(Age_str) %>% 
  slice(1) %>% 
  select(c(Age_str, hoslem_denver, hoslem_sb2, c_denver, c_sb2)) %>% 
  mutate(
    # hoslem
    hoslem_sb2_cat = ifelse(hoslem_sb2 > 0.05, "calibrated", "non-calibrated"),
    hoslem_denver_cat = ifelse(hoslem_denver > 0.05, "calibrated", "non-calibrated"),
    
    # c-index
    better_c = ifelse(c_denver - c_sb2 < 0, "sb2", "denver"),
    
    # weigths
    w = -log10(hoslem_sb2) - (-log10(hoslem_denver)),
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

##ON TEST SET
#calculating the weighed score at test set and stats
test_path <- paste0(PARAM$folder.data, "input_real/scoring_nohide/pheno_scoring_nohide.csv")        # add path to test data
sb2_path <- paste0(PARAM$folder.data, "TEAMS/SB2/output_ori/scores.csv")        # add path to sb2 scores
denver_path <- paste0(PARAM$folder.data, "TEAMS/DenverFINRISKHacky/output_ori/scores.csv") 

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
  select(c(SampleID,surv, Age, Age_str, Score_sb2, Score_denver))

# Join datasets and calculate mean
res <- 
  left_join(data, data_age[, c("Age_str", "hoslem_denver", "hoslem_sb2", "c_denver", "c_sb2")], by = "Age_str") %>% 
  select(-surv) %>%
  select(-Age)

ensemble <- 
  left_join(res, weigths[, c("Age_str", "w_sb2", "w_denver")], by = "Age_str") %>% 
  mutate(
    mean = (Score_sb2 + Score_denver) / 2,
    w_mean = (w_sb2 * Score_sb2 + w_denver * Score_denver) / (w_sb2 + w_denver)
  )
#SampleID <- data$SampleID
#ensemble<- cbind(ensemble1,SampleID)


#calculate C-index and hoslem based on weighed-mean
# load data
print("Load data")
S.test <- read.csv(file = test_path 
                   ,row.names=1)
endpoints <- c("Event_time", "Event")
df.test <- elim(S.test)

surv.object=Surv(df.test$Event_time,df.test$Event)
df.surv <- data.frame(SampleID=rownames(df.test),True.Score=surv.object) #need to create with Sample ID, to make sure that the order of submitted score is align with true score
range01 <- function(x){(x-min(x, na.rm = TRUE))/(max(x, na.rm = TRUE)-min(x, na.rm = TRUE))}
if (length(unique(ensemble$w_mean))==1){
  ensemble$w_mean=ensemble$w_mean
} else {
  ensemble$w_mean = range01(ensemble$w_mean)
}

#merge true score and predicted score so the order is compatible

data_merge <- merge(ensemble, df.surv, by = c("SampleID"))
evaluation <- teststats(data_merge$True.Score,as.numeric(data_merge[,"w_mean"]))
harrell_c <-evaluation$C
Hosmer_lemeshow<-evaluation$hoslem_p$pval
df.eval_weighed <- data.frame("harrell_c"=evaluation$C,"Hosmer_lemeshow"=evaluation$hoslem_p$pval, row.names = "weighed")

#sanity check
evaluation <- teststats(data_merge$True.Score,as.numeric(data_merge[,"mean"]))
harrell_c <-evaluation$C
Hosmer_lemeshow<-evaluation$hoslem_p$pval
df.eval_mean <- data.frame("harrell_c"=evaluation$C,"Hosmer_lemeshow"=evaluation$hoslem_p$pval, row.names = "means")

output <- rbind(df.eval_weighed,df.eval_mean)
head(output)
write.csv(output, paste0(PARAM$folder.data,"/results/pair_weighedAgeStr_eval.csv"))