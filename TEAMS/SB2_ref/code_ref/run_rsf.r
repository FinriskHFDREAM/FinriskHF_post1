require(survival)
require(survcomp)
#require(microbiome)
require(randomForestSRC)

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
model_version <- "ref_lassonophylum_RSF"
code_path <- "/homes/perawija/data_pande/DREAM/FINRISK2022/post_challenge/FinriskHF_post1/TEAMS/SB2_ref/code_ref"

# specify input and output directories
inputdir <- file.path(basedir, inputdir)
outputdir <- file.path(basedir, "TEAMS", team, model_version)


models <- readRDS(file.path(outputdir, "models.rds"))
train <- readRDS(file.path(outputdir, "train.rds"))
test <- readRDS(file.path(outputdir, "test.rds"))

# Loop through each model in models list
for (i in seq_along(models)) {
  
  # Get the selected betas from the model and the model name
  selected_betas <- models[[i]]
  model_name <- names(models)[i]
  
  # Print selected betas for tracking
  print("Selected betas:")
  print(names(selected_betas))
  
  # Print columns available in train data
  print("Train columns:")
  print(colnames(train))
  
  # Subset the train data to include selected features and target variables
  train_ <- train[, match(c("Event", "Event_time", names(selected_betas)), colnames(train))]
  
# Features selected from the model
features_used <- names(selected_betas)

# Build training data that includes event info and selected features
train_model_data <- train[, c("Event_time", "Event", features_used)]

# Fit the RSF model
rsf_model <- rfsrc(
  Surv(Event_time, Event) ~ ., 
  data = train_model_data,
  ntree = 1000,
  importance = TRUE,
  time.interest = seq(0, max(train$Event_time), length.out = 100)
)

# Prepare test data with the same features
test_ <- test[, features_used, drop = FALSE]

# Predict survival
prediction <- predict(rsf_model, newdata = test_, proximity = TRUE)

# Extract predicted risk at 15 years
closest_time <- which.min(abs(rsf_model$time.interest - 15))
pred <- 1 - prediction$survival[, closest_time]

print("Prediction object:")
print(str(prediction))

# Check available time points
print("RSF model time points:")
print(rsf_model$time.interest)

# Then extract survival at t=15
time_index <- which.min(abs(rsf_model$time.interest - 15))
print(paste("Using time index:", time_index, "which corresponds to time", rsf_model$time.interest[time_index]))
  # Print feature information for reference
  print("Train features used:")
  print(features_used)

  print("Test features available:")
  print(colnames(test_))

  print(paste0("Predicting with ", model_name))
  
  # Save the trained Random Survival Forest model
  saveRDS(rsf_model, file = file.path(outputdir, paste0("rsf_", model_name, ".rds")))
  
  # Save predicted risk scores as a CSV file
  scores <- data.frame(
    SampleID = rownames(test),
    Score = pred
  )
  write.csv(scores, quote = FALSE, row.names = FALSE,
            file = file.path(outputdir, paste0("scores_", model_name, ".csv")))
}


end <- Sys.time()
time <- difftime(end, start, units = "mins")
print(time)