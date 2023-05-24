library(survival)
library(tidyr)
library(dplyr) # %>%
library(survival)
#library(ResourceSelection)
#library(GGally)
#library(StepReg)
library(dplyr) # %>%
library(microbiome)
library(Hmisc) 


args=(commandArgs(TRUE))
PARAM <- list()
#this argument is to specify the path of input folder
PARAM$folder.R <- paste0(args[1])
PARAM$folder.data <- paste0(PARAM$folder.R, "/")

ranked_teams <- read.csv(paste0(PARAM$folder.data, "team_list.csv"))
leaderboard <- read.csv(paste0(PARAM$folder.data, "Final_leaderboard.csv"))
outdir <- paste0(args[2])
#combine the leaderboard and teamnames to make reference table with ranked harrel'sC
df_ref <- merge(ranked_teams, leaderboard, by = "id", all.x = TRUE)  %>% arrange(-harrell_c)
#df_ref_h <- merge(ranked_teams, leaderboard, by = "id", all.x = TRUE)  %>% arrange(-hoslem_test)

#merge bootstrap score from all teams

list_score <- list.files(PARAM$folder.R, pattern = "scores.csv", recursive = T, full.names = T) 
list_team <- list_score %>%
    strsplit( "/" ) %>%
    sapply( "[", 9 )
score_df <-	list_score %>% lapply(read.csv)
Scores_df1 <- lapply(score_df, "[", c("Score")) %>% bind_cols()
colnames(Scores_df1) <- c(list_team)
#read 1 scores file as csv to asign sampleID
scoreEX <- read.csv(file = paste0(PARAM$folder.data,
                                  "SB2_Final_Submission/scores.csv"))
rownames(Scores_df1) <-scoreEX$SampleID

#order columns based on list
col_order <- df_ref$folder_name
#print(col_order)
Scores_df <- Scores_df1#[, col_order]
#write.csv(Scores_df,paste0(outdir,"/results/join_all.score.csv"))

#calculate the average Harrel's C score from combinations 2,3,4,5 combination

#2 combinations
# Get all possible combinations of column names
#col_combinations2 <- combn(names(Scores_df), 2, simplify = FALSE)

# Calculate average of each combination for each row
#averages <- t(sapply(col_combinations2, function(cols) rowMeans(Scores_df[cols])))

# Create output data frame with column names as combinations
#output2 <- data.frame(averages, check.names = FALSE)
#colnames(output2) <- sapply(col_combinations2, paste, collapse = "_")

#calculate average score for combinations of interest
col_combinations <- c("SB2_Final_Submission DenverFINRISKHacky_Final",
"SB2_Final_Submission DenverFINRISKHacky_Final Yuanfang_Guan_and_Hanrui_Zhang_Final_Submission",
"SB2_Final_Submission DenverFINRISKHacky_Final Metformin_121_6",
"SB2_Final_Submission DenverFINRISKHacky_Final TristanF_Final_Submission",
"SB2_Final_Submission DenverFINRISKHacky_Final Yuanfang_Guan_and_Hanrui_Zhang_Final_Submission Metformin_121_6", 
"SB2_Final_Submission DenverFINRISKHacky_Final Yuanfang_Guan_and_Hanrui_Zhang_Final_Submission Metformin_121_6 TristanF_Final_Submission", 
"SB2_Final_Submission DenverFINRISKHacky_Final Yuanfang_Guan_and_Hanrui_Zhang_Final_Submission Metformin_121_6 TristanF_Final_Submission UTK_Bioinformatics_Final_Submission",
"SB2_Final_Submission DenverFINRISKHacky_Final Yuanfang_Guan_and_Hanrui_Zhang_Final_Submission Metformin_121_6 TristanF_Final_Submission UTK_Bioinformatics_Final_Submission Pasolli-team-01"
)


# Calculate the means for each row from different combinations of columns
means <- sapply(col_combinations, function(cols) {
  selected_cols <- unlist(strsplit(cols, " "))
  rowMeans(Scores_df[, selected_cols])
})
SampleID <- scoreEX$SampleID
means2<- cbind(means,SampleID)



range01 <- function(x){(x-min(x, na.rm = TRUE))/(max(x, na.rm = TRUE)-min(x, na.rm = TRUE))}

#functions
#remove unwanted variables
elim <- function(df){ 
  df$Event_time<- as.numeric(df$Event_time)
  # exclude any minus or missing value
  df <- subset(df, Event_time > 0 & Event_time < 17)
  df <- subset(df, !is.na(Event_time))
  df <- subset(df, !is.na(Event)) 
  df = df %>% dplyr::select(!PrevalentHFAIL) 
}

grpkm <- function(grps,surv,weights=NULL) {
  if(!is.Surv(surv)) stop('Please use a survival object')
  res <- NULL
  for(i in sort(unique(grps))) {
    grp <- grps==i
    kpla <- survfit(surv[grp,]~1,weights=weights[grp])
    last <- length(kpla$time)
    n <- sum(grp)
    events <- sum(kpla$n.event)
    KM <- kpla$surv[last]
    std.err <- kpla$stderr
    res <- rbind(res,cbind(n,events,KM,std.err))
  }
  res
}

#in-house hoslem test by Aki Havulinna @ THL
tiles <- c('','','tertiles','quartiles','quintiles','sextiles','septiles','octiles','noniles','deciles')
HosLem.test <- function(surv,pred,plot=FALSE,DF.reduce=2) {
  # Cook-Ridker version test uses DF.reduce=2 (default)
  # Cook NR, Ridker PM. (2009) Ann Intern Med. 150(11): 795-802
  # D'Agostino-Nam version uses DF.reduce=1
  # About the differences: http://hdl.handle.net/1773/22648
  # (Guffey D., Hosmer-Lemeshow goodness-of-fit test: Translations to the Cox Proportional Hazards Model)
  # if plot is a name, a jpg plot is saved with this name (.jpg will be added)
  if(!DF.reduce%in%c(1,2)) stop('Please specify DF.reduce = 1 or 2')
  if(!is.Surv(surv)) stop('Please use a survival object')
  version <- c('D\'Agostino-Nam','Cook-Ridker')[DF.reduce]
  grp <- as.numeric(cut2(pred,g=10))
  pj <- as.numeric(by(pred,grp,mean))
  idx <- TRUE
  while(any(idx)&length(pj)>3) {
    kmj <- 1-grpkm(grp,surv)[,3] # failure, not survival
    pj <- as.numeric(by(pred,grp,mean))
    nj <- table(grp)
    idx <- grp%in%which(nj<5|pj*nj<5|nj*kmj<5)&grp<max(grp)
    if(any(idx)) {
      grp[idx] <- grp[idx]+1
      grp <- match(grp,sort(unique(grp)))
    }
  }
  dan_chi2 <- sum(nj*(kmj-pj)^2/(pj*(1-pj)))
  pval <- pchisq(dan_chi2,length(pj)-DF.reduce,lower.tail=FALSE)
  plotname <- ''
  if(plot!=FALSE) {
    if(is.character(plot)) { plot <- paste0(plot,'.jpg'); jpeg(filename=plot); plotname <- plot }
    barplot(height=rbind(pj*100,kmj*100),beside=T,names.arg=1:length(pj),
            main=paste('Hosmer-Lemeshow (',version,') test\np-value =',format.pval(pval,digits=4),sep=''),
            xlab=xlab <- paste('Risk', tiles[length(pj)]),ylab='Risk, %',col=c('black','lightgray'))
    if(is.character(plot)) dev.off()
  }
  list(version=version,quantiles=tiles[length(pj)],nj=nj,kmj=kmj,pj=pj,chi2=dan_chi2,pval=pval,plotname=plotname)
}


#scoring 
teststats <-  function(S,predicted) {
  # S is the actual survival object that is needed, eg "y" from a fitted model
  C <- rcorr.cens(-predicted,S,outx=FALSE)[1] # in rcorr.cens, when using a Survival object, the sign of the predicted score has to be changed
  print(C)
  hoslem_p <- HosLem.test(S,predicted,plot="post") # fit$y is the actual survival object that is needed 
  list(C=C,hoslem_p=hoslem_p)
}


#create survival object using surv package from the test dataset (or validations dataset during validation phase)
#survival object contains informations of time-to-even data and censoring information
# load data
print("Load data")
S.test <- read.csv(file = paste0(outdir, 
                                 "/input_real/scoring_nohide/pheno_scoring_nohide.csv") 
                   ,row.names=1)
endpoints <- c("Event_time", "Event")
df.test <- elim(S.test)

surv.object=Surv(df.test$Event_time,df.test$Event)
df.surv <- data.frame(SampleID=rownames(df.test),True.Score=surv.object) #need to create with Sample ID, to make sure that the order of submitted score is align with true score


#calculate the Harrel'sC and hoslem for all columns
#merge the true and average score
data_merge <- merge(means2, df.surv, by = c("SampleID"))
#write.csv(data_merge,paste0(outdir,"/results/merge_means.csv"))
#calculate !
# Create a new dataframe to store the results
output <- data.frame(harrell_c = numeric(), Hosmer_lemeshow = numeric())
# Calculate addition and subtraction for each column
for (col in colnames(data_merge)[-c(1,10)]) {
  print(col)
  f (length(unique(as.numeric(data_merge[,col])))==1){
  data_merge[,col]=as.numeric(data_merge[,col])
} else {
  data_merge[,col] = range01(as.numeric(data_merge[,col]))
}
  evaluation <- teststats(data_merge$True.Score,as.numeric(data_merge[,col]))
  harrell_c <-evaluation$C
  Hosmer_lemeshow<-evaluation$hoslem_p$pval
  df.eval <- data.frame("harrell_c"=evaluation$C,"Hosmer_lemeshow"=evaluation$hoslem_p$pval, row.names = col)
  output <- rbind(output,df.eval)
}

write.csv(paste0(outdir,"/results/pair_average_eval.csv"))
