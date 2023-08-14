#source:https://codeocean.com/capsule/7326564/tree/v2 (sc1_bootstrap_analysis.R)
#Load Bootstrap samples
# Finally we compute the Bayes factor over the bootstraps for
# consecutively ranked teams
# B:= [sum(score(teamA) < score(teamB))] / [sum(score(teamA) > score(teamB))]
# we consider a Bayes factor larger than 3 significant 

library(survival)
library(tidyr)
library(dplyr) # %>%
library(boot)
library(purrr)
library(readr)

args=(commandArgs(TRUE))
PARAM <- list()
#this argument is to specify the path of input folder
PARAM$folder.R <- paste0(args[1])
PARAM$folder.data <- paste0(PARAM$folder.R, "/")

#functions
#compute bayes factor (https://github.com/Sage-Bionetworks/challengescoring/blob/develop/R/bootstrap.R)
computeBayesFactor <- function(bootstrapMetricMatrix,
                               refPredIndex,
                               invertBayes){
  
  M <- as.data.frame(bootstrapMetricMatrix - bootstrapMetricMatrix[,refPredIndex])
  K <- apply(M ,2, function(x) {
    k <- sum(x >= 0)/sum(x < 0)
    # Logic handles whether reference column is the best set of predictions.
    if(sum(x >= 0) > sum(x < 0)){
      return(k)
    }else{
      return(1/k)
    }
  })
  K[refPredIndex] <- 0
  if(invertBayes == T){K <- 1/K}
  return(K)
}

BF_test <- function(team_a,team_b,bootstrap_table){
  team_a <- as.character(team_a)
  team_b <- as.character(team_b)
  bf = bootstrap_table[c(team_b,team_a)]
  bf_result <-computeBayesFactor(bf,2,invertBayes = FALSE)
  return(bf_result[1])
}


ranked_teams <- read.csv(paste0(PARAM$folder.data, "team_list.csv"))
leaderboard <- read.csv(paste0(PARAM$folder.data, "Final_leaderboard.csv"))

#combine the leaderboard and teamnames to make reference table with ranked harrel'sC
df_ref <- merge(ranked_teams, leaderboard, by = "id", all.x = TRUE)  %>% arrange(-harrell_c)
df_ref_h <- merge(ranked_teams, leaderboard, by = "id", all.x = TRUE)  %>% arrange(-hoslem_test)

#merge bootstrap score from all teams
list_boots <- list.files(PARAM$folder.R, pattern = "\\_boots.csv", recursive = T, full.names = T) 
print(list_boots)
list_team <- list_boots %>%
    strsplit( "/" ) %>%
    sapply( "[", 9 )
boots_df <- list_boots %>% lapply(read.csv)
print(list_team)

#Harell's C
harrell_df <- lapply(boots_df, "[", c("harrell_c")) %>% bind_cols()
colnames(harrell_df) <- c(list_team)
write.csv(harrell_df,file.path(PARAM$folder.data,"bootstrap_Harrell_c_baseline0814.csv"))
harrell_df %>% 
    gather(teams,Harell) %>% group_by(teams) %>%
    summarise(mean_Harell_BS = mean(Harell),
              std_Harell_BS = sd(Harell)) %>% arrange(-mean_Harell_BS) %>%
    write.csv(file.path(PARAM$folder.data,"bootstrap_Harrell_csummary.csv"))


Bayes_factors <- tibble(team_A = df_ref$folder_name[-10], team_B = df_ref$folder_name[-1]) %>%
    rowwise() %>%
    mutate(Bayes_factor = map2_dbl(team_A,team_B,BF_test,harrell_df))

Bayes_factors %>% mutate(Bayes_factor = round(Bayes_factor,2),
                         Bayes_factor = ifelse(
                             is.infinite(Bayes_factor),
                             ">500",as.character(Bayes_factor))) %>%
    write_tsv(file.path(PARAM$folder.data,"harrell_c_bayes_factors.tsv"),col_names = T)

#Hoslem
hoslem_df <- lapply(boots_df, "[", c("hoslem_test")) %>% bind_cols()
colnames(hoslem_df) <- c(list_team)
write.csv(hoslem_df,file.path(PARAM$folder.data,"bootstrap_hoslem_baseline0814.csv"))
hoslem_df %>% 
    gather(teams,Hoslem) %>% group_by(teams) %>%
    summarise(mean_Hoslem_BS = mean(Hoslem),
              std_Hoslem_BS = sd(Hoslem)) %>% arrange(-mean_Hoslem_BS) %>%
    write.csv(file.path(PARAM$folder.data,"bootstrap_hoslem_summary.csv"))
Bayes_factors <- tibble(team_A = df_ref_h$folder_name[-10], team_B = df_ref_h$folder_name[-1]) %>%
    rowwise() %>%
    mutate(Bayes_factor = map2_dbl(team_A,team_B,BF_test,harrell_df))#%>%
    #mutate(t.test = map2_dbl(team_A,team_B,compute_ttest,hoslem_df))

Bayes_factors %>% mutate(Bayes_factor = round(Bayes_factor,2),
                         Bayes_factor = ifelse(
                             is.infinite(Bayes_factor),
                             ">500",as.character(Bayes_factor))) %>%
    write_tsv(file.path(PARAM$folder.data,"hoslem_bayes_factors.tsv"),col_names = T)



##with all possible pairs
pairs_df <- tidyr::expand_grid(df_ref$folder_name,df_ref$folder_name)
names(pairs_df) <- c("team_A","team_B")
Bayes_factors <- pairs_df %>%
  rowwise() %>%
  mutate(Bayes_factor = map2_dbl(team_A,team_B,BF_test,harrell_df))

Bayes_factors %>% mutate(Bayes_factor = round(Bayes_factor,2),
                         Bayes_factor = ifelse(
                           is.infinite(Bayes_factor),
                           ">500",as.character(Bayes_factor))) %>%
  write_tsv(file.path(PARAM$folder.data,"harrell_c_bayes_factors_allpairs_baseline3.tsv"),col_names = T)

pairs_df <- tidyr::expand_grid(df_ref_h$folder_name,df_ref_h$folder_name)
names(pairs_df) <- c("team_A","team_B")
Bayes_factors <- pairs_df %>%
  rowwise() %>%
  mutate(Bayes_factor = map2_dbl(team_A,team_B,BF_test,hoslem_df))

  
Bayes_factors %>% mutate(Bayes_factor = round(Bayes_factor,2),
                         Bayes_factor = ifelse(
                           is.infinite(Bayes_factor),
                           ">500",as.character(Bayes_factor))) %>%
  write_tsv(file.path(PARAM$folder.data,"hoslem_bayes_factors_allpairs_baseline3.tsv"),col_names = T)
