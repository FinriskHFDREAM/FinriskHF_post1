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

The challenge participants model were available as the singularity image. Due to size limitations, the link of singularity image for reproducible analysis is provided in Supplementary Table 3 and can be accessed through the synapse platform (step is given [here](https://help.synapse.org/docs/Downloading-Data-Programmatically.2003796248.html).

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

## C. Post challenge analysis

## D. Reproducing the Figures in mansucript




