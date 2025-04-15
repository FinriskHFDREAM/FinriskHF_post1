# Argument parsing
#.libPaths(c(.libPaths(), "/homes/perawija/opt/miniconda3/lib/R/library"))
args <- commandArgs(trailingOnly = TRUE)

basedir <- args[1]   #/home/joselinares/git/DREAM-FINRISK
inputdir <- args[2]  #input_real
subset <- args[3]    #scoring

#basedir <- "/home/joselinares/git/DREAM-FINRISK"
#inputdir <- "input_real"
#subset <- "scoring"
#code_path <- "/home/joselinares/git/DREAM-FINRISK/src/refining_v2/"

# Inner arguments
team <- "SB2_ref"
model_version <- "ref_lassonophylum_CLR"
code_path <- "/homes/perawija/data_pande/DREAM/FINRISK2022/post_challenge/FinriskHF_post1/TEAMS/SB2_ref/code_ref"

# specify input and output directories
inputdir <- file.path(basedir, inputdir)
outputdir <- file.path(basedir, "TEAMS", team, model_version)

print(inputdir)
print(outputdir)

# Creating output folder
dir.create(outputdir, recursive = TRUE)

# Load packages and functions
start <- Sys.time()
source(file.path(code_path, "utils_seeds/requirements.r"))
source(file.path(code_path, "utils_seeds/importPseq.r"))
source(file.path(code_path, "utils_seeds/prepro_functions.r"))
source(file.path(code_path, "utils_seeds/co-abundances.r"))
source(file.path(code_path, "utils_seeds/get_scores.r"))
source(file.path(code_path, "utils_seeds/fit_model.r"))
source(file.path(code_path, "utils_seeds/clr_transform.r"))

print("Start preprocessing...")

# Preprocess
# ======
source(file.path(code_path, "utils_seeds/preprocessing_clr.r"))

# Fit model
# ========
method = c("lasso","enet")
cvrts <- c("Age", "BodyMassIndex",
           "SystolicBP", "NonHDLcholesterol",
           "Sex")

print("Fitting the model ...")

models <- fit_biospear(data = train,
                       biomarkers = setdiff(colnames(train),
                                            c(cvrts, colnames(y_train))),
                       surv = c("Event_time", "Event"),
                       cvrts = cvrts,
                       inter = FALSE,
                       methods = method)
colnames(train)
saveRDS(models, file = file.path(outputdir, "models.rds"))

# Build cox with all train data
# ====
risk <- function(model, newdata, time) {
  as.numeric(1-summary(
    survfit(model, newdata = newdata, se.fit = F, conf.int = F), 
    times = time, extend = TRUE)$surv)
}

for (i in seq_along(models)) {
  selected_betas <- models[[i]]
  model_name <- names(models)[i]
  print("Selected betas:")
  print(names(selected_betas))
  print("Train columns:")
  print(colnames(train))
    train_ <- train[, match(
    c("Event", "Event_time", names(selected_betas)),
    colnames(train)
  )]
  surv <- "Surv(Event_time, Event) ~"
  f <- as.formula(paste0(surv, paste(names(selected_betas),collapse='+')))
  model <- coxph(f, train_)

  # Predict risk
  # ===
  pred <- risk(model, test, time = 15)
  print(paste0("Predicting with ", model_name))
  print(coef(model))
  saveRDS(model, file = file.path(outputdir, paste0("cox_",model_name,".rds")))
  # Save q
  scores <- data.frame(
      SampleID = rownames(test),
      Score = pred
  )
  write.csv(scores, quote = FALSE, row.names = FALSE,
            file = file.path(
              outputdir, 
              paste0("scores_", model_name, ".csv")))

  
}

end <- Sys.time()
time <- difftime(end, start, units = "mins")
print(time)