# Copyright: Aki Havulinna @ THL, modified by Pande Putu Erawijantari 

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
#the input folder structure is similar to DreamHF.zip both for synthetic and real dataset
PARAM$folder.R <- paste0(args[1])
PARAM$folder.out <- paste0(args[2])
#please avoid using (.) in your team name
PARAM$folder.data <- paste0(PARAM$folder.R, "/")
PARAM$folder.result <- paste0(PARAM$folder.data, PARAM$folder.out,"/output/")
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
  hoslem_p <- HosLem.test(S,predicted,plot="submission_phase") # fit$y is the actual survival object that is needed 
  list(C=C,hoslem_p=hoslem_p)
}


#create survival object using surv package from the test dataset (or validations dataset during validation phase)
#survival object contains informations of time-to-even data and censoring information
# load data
print("Load data")
S.test <- read.csv(file = paste0(PARAM$folder.data, 
                                 "test/pheno_test.csv") 
                   ,row.names=1)
endpoints <- c("Event_time", "Event")
df.test <- elim(S.test)

surv.object=Surv(df.test$Event_time,df.test$Event)
df.surv <- data.frame(SampleID=rownames(df.test),True.Score=surv.object) #need to create with Sample ID, to make sure that the order of submitted score is align with true score

#read score file from submitted algorithm
scores <- read.csv(file = paste0(PARAM$folder.result, 
                                  "scores.csv"))
if (length(unique(scores$Score))==1){
  scores$Score=scores$Score
} else {
  scores$Score = range01(scores$Score)
}

#merge true score and predicted score so the order is compatible
data_merge <- merge(scores, df.surv, by = c("SampleID"))
evaluation <- teststats(data_merge$True.Score,data_merge$Score)
df.eval <- data.frame("harrell_c"=evaluation$C,"hoslem_test"=evaluation$hoslem_p$pval, row.names = "value")

write.csv(df.eval, file=paste0(PARAM$folder.result, 
                                     "stats.csv"), 
          quote=FALSE, row.names=FALSE)
