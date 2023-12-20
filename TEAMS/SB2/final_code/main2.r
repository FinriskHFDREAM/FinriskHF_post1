# Argument parsing
args <- commandArgs(trailingOnly = TRUE)

basedir <- args[1]   #/home/joselinares/git/DREAM-FINRISK
inputdir <- args[2]  #input_real
subset <- args[3]    #scoring

#basedir <- "/home/joselinares/git/DREAM-FINRISK"
#inputdir <- "input_real"
#subset <- "scoring"
#code_path <- "/home/joselinares/git/DREAM-FINRISK/src/refining/"

set.seed(20233001)
# Inner arguments
team <- "SB2"
model_version <- "final_0814"
code_path <- "/mnt"

# specify input and output directories
inputdir <- file.path(basedir, inputdir)
outputdir <- file.path(basedir, "TEAMS", team, model_version)

print(inputdir)
print(outputdir)

# Creating output folder
dir.create(outputdir, recursive = TRUE)

# Load packages and functions
start <- Sys.time()
source(file.path(code_path, "utils/requirements.r"))
source(file.path(code_path, "utils/importPseq.r"))
source(file.path(code_path, "utils/prepro_functions.r"))
source(file.path(code_path, "utils/co-abundances.r"))
source(file.path(code_path, "utils/get_scores.r"))
source(file.path(code_path, "utils/fit_model.r"))

print("Start preprocessing...")

# Preprocess
# ======
source(file.path(code_path, "utils/preprocessing.r"))

# Fit model
# ========
method = c("lasso")
cvrts <- c("Age", "BodyMassIndex",
           "SystolicBP", "NonHDLcholesterol",
           "Sex", "disbiosis")

print("Fitting the model ...")
models <- fit_biospear(data = train,
                       biomarkers = setdiff(colnames(train),
                                            c(cvrts, colnames(y_train))),
                       surv = c("Event_time", "Event"),
                       cvrts = cvrts,
                       inter = FALSE,
                       methods = method)

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

  train_ <- train[, match(
    c("Event", "Event_time", names(selected_betas)),
    colnames(train)
  )]

  surv <- "Surv(Event_time, Event) ~"
  f <- as.formula(paste0(surv, paste(names(selected_betas), collapse='+')))
  model <- coxph(f, train_)
  print(coef(model))
  saveRDS(model, file = file.path(outputdir, "cox_final.rds"))

  # Predict risk
  # ===
  pred <- risk(model, test, time = 10)
  print(coef(model))
  print(summary(model))
  saveRDS(summary(model), file = file.path(outputdir, paste0("summarycox_",model_name,".rds")))
  saveRDS(model, file = file.path(outputdir, paste0("cox_",model_name,".rds")))
  # Save scores
  scores <- data.frame(
      SampleID = rownames(test),
      Score = pred
  )
  write.csv(scores, quote = FALSE, row.names = FALSE,
            file = file.path(
              outputdir, 
              paste0("scores_", model_name, ".csv")))

  print(paste0("Predicting with ", model_name))
}

end <- Sys.time()
time <- difftime(end, start, units = "mins")
print(time)
