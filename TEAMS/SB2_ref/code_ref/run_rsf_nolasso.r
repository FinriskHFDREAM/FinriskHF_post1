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

# === Define the features to use ===
# Exclude Event and Event_time from predictors
features_used <- setdiff(colnames(train), c("Event", "Event_time"))

# === Build training and test datasets ===
train_model_data <- train[, c("Event_time", "Event", features_used)]
test_ <- test[, features_used, drop = FALSE]

# === Fit the RSF model ===
rsf_model <- rfsrc(
  Surv(Event_time, Event) ~ ., 
  data = train_model_data,
  ntree = 1000,
  importance = TRUE,
  time.interest = seq(0, max(train$Event_time), length.out = 100)
)

# === Predict survival on test set ===
prediction <- predict(rsf_model, newdata = test_, proximity = TRUE)

# === Extract predicted risk at 15 years ===
closest_time <- which.min(abs(rsf_model$time.interest - 15))
pred <- 1 - prediction$survival[, closest_time]

# === Diagnostics ===
print("Prediction object:")
print(str(prediction))

print("RSF model time points:")
print(rsf_model$time.interest)

print(paste("Using time index:", closest_time, 
            "which corresponds to time", rsf_model$time.interest[closest_time]))

print("Train features used:")
print(features_used)

print("Test features available:")
print(colnames(test_))

# === Save model and risk scores ===
model_name <- "rsf_full_features"

saveRDS(rsf_model, file = file.path(outputdir, paste0("rsf_", model_name, ".rds")))

scores <- data.frame(
  SampleID = rownames(test),
  Score = pred
)

write.csv(scores, quote = FALSE, row.names = FALSE,
          file = file.path(outputdir, paste0("scores_", model_name, ".csv")))



end <- Sys.time()
time <- difftime(end, start, units = "mins")
print(time)