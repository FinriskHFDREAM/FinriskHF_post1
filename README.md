# Repository for sharing and testing code on model refinement

## Directory information

* input_syn : placement of sythetic dataset, only training and test set will be shared. Inside the directory you can find different versions of synthetic data (if updated, for example with additional metadata).

* TEAMS: the directory to coodinate your code for testing in real dataset in model refinement

* code: directory for general usage code and best final model results from refinement

* results: placement for figures, tables and general outputs for papers, the output of testing refine model will be placed in TEAMS/(name_of_team)/output_(date)

* input_real: the real dataset (will not be shared in git)

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

