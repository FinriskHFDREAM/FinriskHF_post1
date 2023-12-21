# The Heart Failure Prediction: Microbiome-FINRISK DREAM Challenge

The repository contains informations related to The Heart Failure Prediction: Microbiome-FINRISK DREAM Challenge, including:
* The submissions retrieval from synapse throughout the challenge
* The code for evaluations
* Post-challenge record on model ensembles and refinement
* Post-challenge supporting analysis

## A. Overview directory information

* input_syn : placement of sythetic dataset (the dataset are available through the THL biobank upon submission of a research plan and signing a daata transfer agreement, https://thl.fi/en/web/thl-biobank/for-researchers/ ).

* TEAMS: the directory to coordinate  code for testing in real dataset during model refinement process with 2 top performing teams

* code: directory for general usage code and reproducing figures for manuscript

* challenge_infrastructure: submissions retrieval, evaluations code, baseline model

* Results: placement for figures, tables and general outputs for papers, the output of testing refine model will be placed in TEAMS/(name_of_team)/(model_version)

* input_real: the real dataset (the dataset are available through the THL biobank upon submission of a research plan and signing a daata transfer agreement, https://thl.fi/en/web/thl-biobank/for-researchers/ )

The challenge participants model were available as the singularity image. Due to size limitations, the link of singularity image for reproducible analysis is provided in Supplementary Table 3 and can be accessed through the synapse platform (step is given [here](https://help.synapse.org/docs/Downloading-Data-Programmatically.2003796248.html). We are aware that **Singularity** is now known as **Apptainer**. The use of singularity and apptainer is interchengable and therefore, here we refered it as **Singularity**.

IMPORTANT: in the code please allow arguments that will point to the input directory (real or synthetic)
The structure of the input directory will be as described as below. The pheno_scoring_nohide.csv and  pheno_test_nohide.csv are for evaluating the individual risk score from the true value (event and event time information available).

```
├── scoring
│   ├── pheno_scoring.csv
│   ├── readcounts_scoring.csv
│   ├── pheno_scoring_nohide.csv
│   └── taxtable.csv
└── train
│   ├── pheno_train.csv
│   ├── readcounts_train.csv
│   └── taxtable.csv
└── test
│   ├── pheno_test.csv
│   ├── readcounts_test.csvc
│   ├── pheno_test_nohide.csv
│   └── taxtable.csv
```

## B. The submissions retrieval, evaluation, and status update from synapse throughout the challenge period

Notes on retriving, evaluating, and status updating of the challenge submission is available in the `challenge_infrastructure` directory with detail instructions in `submission_retrival.md`.
All singularity submission from participants is available in synapse platform and listed in the Supplementary Table 3 of the manuscript.

To ensure a robust ranking of participants, we additionally performed 1000 bootstrap iterations of random sampling on the individual’s risk scores calculated by each model. The evaluation metrics, Harrell’s C-index and Hosmer-Lemeshow p-value, were then re-calculated to generate a distribution of evaluation scores for each submission. We used these metrics to calculate the Bayes factor, using the computeBayesFactor functions from the challenge scoring R package (https://github.com/Sage-Bionetworks/challengescoring/blob/develop/R/bootstrap.R) and comparing them to the top-performing model as well as to the baseline models (Supp. Table 8). 

The code to create the 1000 bootstrap of evaluation metrics from each submission available in `challenge_infrastructure/boot_score.R`. 
The Bayes Factor reported in Supp. Table 8 were  calculated using the code available in `code/BayesFactor_calc.R`. 


## C. Post challenge analysis

Due to participants’ limited access to the real dataset, we also worked with the top two performing teams to further improve the model performance and calibrations for the final submitted models after the challenge formally concluded. Top 2 performing models has submitted different version of model refinement and was evaluated in scoring dataset. All models were available in `TEAMS` directory. Among them we primarly reported are 2 modifications of SB2 final model and the ensemble model among all received model. 

### 1. SB2 model refinement

The SB2 final model R code already containerized in the submitted singularity image listed in Supplementary Table 3.
Command to run the submission during the final phase when evaluations were carried out on scoring dataset is as follow:

```
singularity run SB2_final.sif input_real scoring
```
During the post-challenge phase, 2 more modifications were tested, which are:

a. **SB2 refined model** where the phylum informations were removed from the features for LASSO feature selection due to it representing redundant information with species information. The model can be found in: `TEAMS/SB2/code_v2/main_lasso.r`

b. **SB2 age fixed** almost similar to SB2 refined model above, however in this model, only Age is selected as unpenalized features. The model can be found in: `TEAMS/SB2/code_v2/main_Agefix.r`

Similar singularity image that has been submitted during the challenge can be used to execute the refine model as follow:

```
singularity exec SB2_final.sif Rscript <each script listed above> ./ input_real scoring
```

### 2. Ensemble model

We constructed ensemble models  as shown in Supp. Figure 5, which were generated by taking the mean of prediction scores obtained from the combined individual models. 



## D. Reproducing the Figures in mansucript




