import os
import synapseclient
syn = synapseclient.Synapse()
syn = synapseclient.login()
from challengeutils import annotations
import pandas as pd
from challengeutils import utils
import argparse

#record arguments
parser = argparse.ArgumentParser()
parser.add_argument("-s", "--subID", required=True, help="Submission ID")
parser.add_argument("-n", "--subsname", required=True, help="Submission title")
#parser.add_argument("-c", "--synapseConfig", required=True, help="credentials file")
parser.add_argument("-d", "--workingdir", required=True, help="working directory path")
args = parser.parse_args()
#get the score real dataset
scores = pd.read_csv(os.path.join(args.workingdir,"real_dataset",args.subsname, "output", "stats.csv"))
res = {'harrell_c': scores.loc[0, 'harrell_c'],'hoslem_test': scores.loc[0, 'hoslem_test']}
annotations.annotate_submission(syn, args.subID, res, is_private=False)

#get the score training dataset
scores = pd.read_csv(os.path.join(args.workingdir,"real_dataset_train",args.subsname, "output", "stats.csv"))
res = {'harrell_c_train': scores.loc[0, 'harrell_c'],'hoslem_train': scores.loc[0, 'hoslem_test']}
annotations.annotate_submission(syn, args.subID, res, is_private=False)
utils.change_submission_status(syn, args.subID, "SCORED")

