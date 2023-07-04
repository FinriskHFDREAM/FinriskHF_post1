# The Heart Failure Prediction: Microbiome-FINRISK DREAM Challenge

The repository contains informations related to The Heart Failure Prediction: Microbiome-FINRISK DREAM Challenge, including:
* The submissions retrieval from synapse throughout the challenge
* The code for evaluations
* Post-challenge record on model ensembles and refinement
* Post-challenge supporting analysis

## A. Overview directory information

* input_syn : placement of sythetic dataset, only training and test set will be shared. Inside the directory, you can find different versions of synthetic data (if updated, for example with additional metadata).

* TEAMS: the directory to coordinate  code for testing in real dataset during model refinement process with 2 top performing teams

* code: directory for general usage code and reproducing figures for manuscript

*challenge_infrastructure: submissions retrieval, evaluations code, baseline model

* Results: placement for figures, tables and general outputs for papers, the output of testing refine model will be placed in TEAMS/(name_of_team)/(model_version)

* input_real: the real dataset (will not be shared in git)

Singularity image for reproducible analysis can be accessed through synapse platform (example given [here](https://help.synapse.org/docs/Downloading-Data-Programmatically.2003796248.html) 

IMPORTANT: in the code please allow arguments that will point to the input directory (real or synthetic)
The structure of the input directory will be as described as below. The pheno_scoring_nohide.csv is for evaluating the individual risk score from the true value (event and event time information available).

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
│   ├── readcounts_test.csv
│   └── taxtable.csv
```

## B. The submissions retrieval from synapse throughout the challenge period
