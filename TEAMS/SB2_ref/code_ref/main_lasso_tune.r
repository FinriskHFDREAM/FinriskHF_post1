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
model_version <- "ref_lassonophylum_tune"
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

print("Start preprocessing...")

# Preprocess
# ======
source(file.path(code_path, "utils_seeds/preprocessing.r"))

# Fit model
# ========
method = c("lasso", "enet")
cvrts <- c("Age", "BodyMassIndex",
           "SystolicBP", "NonHDLcholesterol",
           "Sex")

# Define tuning grid
param_grid <- expand.grid(
  method = c("lasso", "enet"),
  ss.nsub = c(100, 200),
  ss.thr = c(0.6, 0.5),
  dfmax = c(70, 50),
  stringsAsFactors = FALSE
)

print("Fitting the model ...")

# Function: Run BMsel with given params
fit_model <- function(data, biomarkers, surv, cvrts, params) {
  fit <- BMsel(
    data = data,
    x = biomarkers,
    y = surv,
    z = cvrts,
    inter = FALSE,
    tt = NULL,
    std.x = TRUE,
    std.i = FALSE,
    std.tt = FALSE,
    method = params$method,
    folds = 10,
    uni.fdr = 0.05,
    uni.test = 1,
    ss.rando = FALSE,
    ss.nsub = params$ss.nsub,
    ss.fsub = 0.5,
    ss.fwer = 1,
    ss.thr = params$ss.thr,
    dfmax = params$dfmax,
    pct.rep = 1,
    pct.qtl = 0.95,
    trace = FALSE
  )
  return(fit)
}


# Function: Evaluate Cox model using C-index
evaluate_model <- function(cox_model, test_data) {
  surv_obj <- Surv(test_data$Event_time, test_data$Event)
  risk_scores <- predict(cox_model, newdata = test_data, type = "risk")
  cindex <- concordance.index(x = risk_scores, surv.time = test_data$Event_time,
                              surv.event = test_data$Event, method = "noether")
  return(cindex$c.index)
}

# Store results
tuned_models <- list()
evaluation_results <- data.frame()

# === Run tuning loop ===
for (i in seq_len(nrow(param_grid))) {
  cat("Running model config", i, "\n")
  params <- param_grid[i, ]
  model_id <- paste0("model_", i)

  fit <- fit_model(
  data = train,
  biomarkers = setdiff(colnames(train), c(cvrts, colnames(y_train))),
  surv = c("Event_time", "Event"),
  cvrts = cvrts,
  params = params
)
  if (length(selected_features) == 0) {
    cat("No features selected for", model_id, "\n")
    next
  }

  # Fit Cox model with selected features
  train_subset <- train[, c("Event_time", "Event", selected_features)]
  formula_str <- paste0("Surv(Event_time, Event) ~ ", paste(selected_features, collapse = "+"))
  cox_model <- coxph(as.formula(formula_str), data = train_subset)

  # Evaluate
  cidx <- evaluate_model(cox_model, test)

  # Store results
  tuned_models[[model_id]] <- list(config = params, model = fit, cox = cox_model)
  evaluation_results <- rbind(evaluation_results, cbind(model_id, params, cindex = cidx))

  # Save model and predictions
  saveRDS(cox_model, file = file.path(outputdir, paste0("cox_", model_id, ".rds")))

  # Save risk scores
  risk_score <- predict(cox_model, newdata = test, type = "risk")
  scores <- data.frame(SampleID = rownames(test), Score = risk_score)
  write.csv(scores, file = file.path(outputdir, paste0("scores_", model_id, ".csv")),
            row.names = FALSE, quote = FALSE)
}

# === Save everything ===
saveRDS(tuned_models, file = file.path(outputdir, "all_tuned_models.rds"))
write.csv(evaluation_results, file = file.path(outputdir, "model_evaluation_summary.csv"),
          row.names = FALSE)

cat("Tuning & evaluation complete. Results saved to:", outputdir, "\n")


end <- Sys.time()
time <- difftime(end, start, units = "mins")
print(time)