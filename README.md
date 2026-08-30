# The Heart Failure Prediction: Microbiome-FINRISK DREAM Challenge

The repository contains informations related to The Heart Failure Prediction: Microbiome-FINRISK DREAM Challenge, including:
* The submissions retrieval from synapse throughout the challenge
* The code for evaluations
* Post-challenge record on model ensembles and refinement
* Post-challenge supporting analysis

## A. Overview directory information

* input_syn : placement of sythetic dataset (the dataset are available through the THL biobank upon submission of a research plan and signing a data transfer agreement, https://thl.fi/en/web/thl-biobank/for-researchers/ ).

* TEAMS: the directory to coordinate  code for testing in real dataset during model refinement process with 2 top performing teams

* code: directory for general usage code and reproducing figures for manuscript

* challenge_infrastructure: submissions retrieval, evaluations code, baseline model

* Results: placement for figures, tables and general outputs for papers, the output of testing refine model will be placed in TEAMS/(name_of_team)/(model_version)

* input_real: the real dataset (the dataset are available through the THL biobank upon submission of a research plan and signing a data transfer agreement, https://thl.fi/en/web/thl-biobank/for-researchers/ )

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

Notes on retriving, evaluating, and status updating of the challenge submission is available in the **[challenge_infrastructure](https://github.com/FinriskHFDREAM/FinriskHF_post1/tree/main/challenge_infrastructure)** directory with detail instructions in **[submission_retrival.md](https://github.com/FinriskHFDREAM/FinriskHF_post1/blob/main/challenge_infrastructure/submission_retrival.md)**.
All singularity submission from participants is available in synapse platform and listed in the **Supplementary Table 3** of the manuscript.

To ensure a robust ranking of participants, we additionally performed 1000 bootstrap iterations of random sampling on the individual’s risk scores calculated by each model. The evaluation metrics, Harrell’s C-index and Hosmer-Lemeshow p-value, were then re-calculated to generate a distribution of evaluation scores for each submission. We used these metrics to calculate the Bayes factor, using the computeBayesFactor functions from the challenge scoring R package (https://github.com/Sage-Bionetworks/challengescoring/blob/develop/R/bootstrap.R) and comparing them to the top-performing model as well as to the baseline models (**Supp. Table 4**). 

The code to create the 1000 bootstrap of evaluation metrics from each submission available in **[challenge_infrastructure/boot_score.R](https://github.com/FinriskHFDREAM/FinriskHF_post1/blob/main/challenge_infrastructure/boot_score.R)**. 
The Bayes Factor reported in **Supp. Table 4** were  calculated using the code available in **[code/BayesFactor_calc.R](https://github.com/FinriskHFDREAM/FinriskHF_post1/blob/main/code/BayesFactor_calc.R)**. 


## C. Post challenge analysis

Due to participants’ limited access to the real dataset, we also worked with the top two performing teams to further improve the model performance and calibrations for the final submitted models after the challenge formally concluded. Top 2 performing models has submitted different version of model refinement and was evaluated in scoring dataset. All models were available in [TEAMS](https://github.com/FinriskHFDREAM/FinriskHF_post1/tree/main/TEAMS) directory. Among them we primarly reported are 2 modifications of SB2 final model and the ensemble model among all received model. 

### 1. SB2 model refinement

The SB2 final model R code already containerized in the submitted singularity image listed in **Supplementary Table 3**.
Command to run the submission during the final phase when evaluations were carried out on scoring dataset is as follow:

```
singularity run SB2_final.sif input_real scoring
```

During the post-challenge phase, 2 more modifications were tested which include the **SB2_v2** and **SB2_v3**. For the *SB2_v2*, we revised the taxonomic filtering criteria by introducing more stringent and interpretable thresholds for detection and prevalence. Instead of filtering taxa based on raw count abundance (present in >80% of samples with counts >2), we applied a relative abundance threshold of 0.1% and retained taxa present in at least 10% of the samples. This ensured that rare and potentially noisy taxa were excluded more robustly. Additionally,  we modified the downstream Gene Set Variation Analysis (GSVA) score calculation to operate on relative abundance data instead of raw counts. This adjustment allowed GSVA to better capture variation across biologically coherent microbial clusters. To reduce redundancy and improve model efficiency, we excluded phylum-level information. Additionally, the Firmicutes-Bacteroidetes ratio was subjected to penalization, rather than being treated as an unpenalized covariate as in the original model. Finally, we simplified the alpha diversity estimation step by focusing on Shannon diversity and observed richness. **SB2_v3** using the similar approach as the **SB2_v2**, however in this model, only Age is selected as unpenalized features. We also created one alternative model based on the Random Survival Forest (RSF). The unpenalized features in all of three model were subjected to features selection via LASSO or EN. We also performed RSF model without feature selection. All of the code for each model variation will be described in details bellow:

a. The code for **SB2_v2_++CoxPh** can be found in the directory below with code implemented LASSO or EN feature selection in parallel : **[TEAMS/SB2/code_ref/main_lasso.r](https://github.com/FinriskHFDREAM/FinriskHF_post1/blob/main/TEAMS/SB2_ref/code_ref/main_lasso_clr.r)**

b. **SB2_v3_++CoxPh** almost similar to SB2_v2, however in this model, only Age is selected as unpenalized features. The code can be found in: **[TEAMS/SB2/code_ref/main_lasso_clrAgeFix.r](https://github.com/FinriskHFDREAM/FinriskHF_post1/blob/main/TEAMS/SB2_ref/code_ref/main_lasso_clrAgeFix.r)**

c. **SB2_v2_++RSF** is RSF implementation based on feature engineering applied in SB2 modelwith LASSO or EN feature selection. The code can be found in: **[TEAMS/SB2/code_ref/main_rsf.r](https://github.com/FinriskHFDREAM/FinriskHF_post1/blob/main/TEAMS/SB2_ref/code_ref/main_rsf.r)**

d. **SB2_v3_++RSF** is similar as model c, with Age as the only penalized feature. The code can be found in: **[TEAMS/SB2/code_ref/main_rsfAgeFix.r](https://github.com/FinriskHFDREAM/FinriskHF_post1/blob/main/TEAMS/SB2_ref/code_ref/main_rsfAgeFix.r)**

e. **SB2_v2_AllRSF** is similar to model c without LASSo or EN penalization. The code can be found in: **[TEAMS/SB2/code_ref/run_rsf_nolasso.r](https://github.com/FinriskHFDREAM/FinriskHF_post1/blob/main/TEAMS/SB2_ref/code_ref/run_rsf_nolasso.r)**

The ++ indicated whether the LASSO or EN feature selection was applied which was implemented in the same code. Similar singularity image that has been submitted during the challenge can be used to execute the refine model as follow:

```
singularity exec SB2_final.sif Rscript <each script listed above> ./ input_real scoring
```

### 2. Ensemble model

We constructed ensemble models and the evaluation metrices was shown in Extended Data Figure 5. The ensembles model were generated by taking the mean of prediction scores obtained from the combined individual models as indicated in the figure.

The code to analyze available in [code/AveScore_allteam.R](https://github.com/FinriskHFDREAM/FinriskHF_post1/blob/main/code/AveScore_allteam.R) and the output is available in [results/pair_average_eval.csv](https://github.com/FinriskHFDREAM/FinriskHF_post1/blob/main/results/pair_average_eval.csv) as also provided in Supp. Table 7. 


## D. Reproducing the Figures in mansucript

- Figure 1A. Geographical distribution across Finland for the individuals within the national FINRISK 2002 cohort was adapted from https://gitlab.com/openfinrisk2002/article-2021-finrisk02-taxonomic-signatures-mortality-risk/-/blob/master/code/create_figure1_map_and_pcoa/finland.R
- Figure 1B. Principal Coordinate Analysis (PCoA) using Bray-Curtis dissimilarity metrics between randomly selected subsets of the data (training, testing, scoring sets) : [code/pcoa_viz.R](https://github.com/FinriskHFDREAM/FinriskHF_post1/blob/main/code/pcoa_viz.R)
- Figure 2B were created according to the Supplementary Table 4 and 6 that was extracted from the output of each model.
- Figure 4. Cumulative incidence of HF over 15 years of follow up, stratified by quintiles of predicted risk score: [code/km_curves.html](https://github.com/FinriskHFDREAM/FinriskHF_post1/blob/main/code/km_curves.html)
- The rest of the figure were created by [code/plotting_final.R](https://github.com/FinriskHFDREAM/FinriskHF_post1/blob/main/code/plotting_final.R).




